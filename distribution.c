#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <nlopt.h>

#include "common.h"
#include "shortest_path.h"
#include "convolution.h"
#include "simulation.h"

#include "distribution.h"

static const double INF = DBL_MAX/2 - 1;
static const unsigned int _seed = INT_MAX - 1;

static NODE *node;
static unsigned int NODE_NUM = 0;
static FILE *fp;
static unsigned int *delay_t = NULL;
static int dnum_t = 0;
static double *pdf = NULL;
static int pdf_len = 0;
static int pdf_cur = 0;
static peerlist *p_ccdf = NULL;
static double sim_rev = 0;
static long sim_delay = 0;
static int sim_delivery = 0;
static int _candidate = -1;
static double best_rev = 0;
static char *rev_test;

#define FIXED_ROUTE
#define SINGLE_SELECT
//#define _DEBUG

#define TSLOT	60
#define KTHRESH	6
#define PRICE 	50
#define COST	20

#define remove_data(data, node)	\
{	\
	FDATA *__last = &(node->buffer[node->buff_cur - 1]);	\
	FDATA *__tmp = data;	\
	data = __last;	\
	__last = __tmp;	\
	node->buff_cur--;	\
}

#ifdef _DEBUG
#define debug(num) \
	do {	\
		int i;	\
		for(i=0; i<num; i++) {	\
			M_NODE *_debug = &(node[i]);	\
			for(i=0; i<_debug->buff_cur; i++) {	\
				FDATA *_tmp = &(_debug->buffer[i]);	\
				if(_tmp->type == FILE_DISTRIBUTION)	\
					printf("(%d) - src:%d, dest: %d\n", __LINE__, _tmp->src, _tmp->dest);	\
			}	\
		}	\
	} while(0)
#define _dprintf	printf
#else
#define debug(num)
#define _dprintf	
#endif

static int int_cmp(const void *n1, const void *n2)
{
	return (*(unsigned int *)n1 - *(unsigned int *)n2);
}

static int nei_cmp(const void *n1, const void *n2)
{
	return (((NEIGHBOR *)n1)->id - ((NEIGHBOR *)n2)->id);
}

static double cal_pdf(unsigned int *array, int num)
{
	if(expect_false(pdf_len < num))
		array_needsize(double, pdf, pdf_len, num, array_zero_init);
	
	int i, time = TSLOT, sum = 0;
	double total_delay = 0;
	pdf_cur = 0;
	for(i=0; i<num; i++) {
		if(array[i] > time) {
			pdf[pdf_cur++] = (double)sum/(double)num;
			time += TSLOT;
			sum = 0;
		}
		sum++;
		total_delay += array[i];
	}
	if(sum > 0)
		pdf[pdf_cur++] = (double)sum/(double)num;

	return (total_delay/(double)num);
}

static void neighbor_wb(unsigned int *delay/*neighbor delay distribution*/, int num/*num of inter contact time record*/, int nei_id/*neighbor id*/, NODE *n/*node*/)
{
	//sort for future use at first...
	qsort(delay, num, sizeof(unsigned int), int_cmp);

	if(expect_false(n->cur == n->num))
		array_needsize(NEIGHBOR, n->nei, n->num, n->num + 1, array_zero_init);

	NEIGHBOR *p = &(n->nei[n->cur++]);
	p->id = nei_id;
	p->delay_average = cal_pdf(delay, num);
	p->num = pdf_cur;
	p->delay_pdf = (double *)calloc(pdf_cur, sizeof(double));
	memcpy(p->delay_pdf, pdf, pdf_cur * sizeof(double));
}

static void get_node_info(unsigned int id)
{
#define NEIGHBOR_THRESHOLD	30
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	unsigned int p_neighbor = 0;
	int cur_t = 0;

	NODE *n = &(node[id-1]);
	array_needsize(NEIGHBOR, n->nei, n->num, 2, array_zero_init);

	while((read = getline(&line, &len, fp)) != -1) {
		unsigned int time, node, neighbor;
		sscanf(line, "%d,%d,%d", &node, &neighbor, &time);

		//the infor of node id has been all handled
		if(node > id) {
			long fpos = ftell(fp) - read;	
			fseek(fp, fpos, SEEK_SET);	//roll back

			//check if we need to write back the neighbor infor before return
			if(cur_t > NEIGHBOR_THRESHOLD)
				neighbor_wb(delay_t, cur_t - 2, p_neighbor - 1, n);

			break;
		}

		//the infor of neighbor p_neighbor has been all handled
		if(p_neighbor != neighbor) {
			if(cur_t > NEIGHBOR_THRESHOLD)
				neighbor_wb(delay_t/*neighbor delay distribution*/, cur_t - 2/*num of inter contact time record*/, p_neighbor - 1/*neighbor id*/, n/*node*/);

			//reset for new neighbor infor
			cur_t = 0;
			p_neighbor = neighbor;
		}

		//record the neighbor info...
		if(expect_false(cur_t == dnum_t))
			array_needsize(unsigned int, delay_t, dnum_t, dnum_t + 1, array_zero_init);
		if(cur_t == 0)
			delay_t[cur_t++] = time;
		else {
			delay_t[cur_t - 1] = time - delay_t[cur_t - 1];
			delay_t[cur_t++] = time;
		}
	}
	free(line);
#undef NEIGHBOR_THRESHOLD
}

static void exit_clean(void)
{
	if(delay_t) {
		free(delay_t);
		delay_t = NULL;
	}

	if(pdf) {
		free(pdf);
		pdf = NULL;
	}

	fclose(fp);
}

static void node_free(void)
{
	if(node) {
		int i;
		for(i=0; i<NODE_NUM; i++) {
			if(node[i].nei) {
				int j;
				for(j=0; j<node[i].cur; j++)
					if(node[i].nei[j].delay_pdf) {
						free(node[i].nei[j].delay_pdf);
						node[i].nei[j].delay_pdf = NULL;
					}
				free(node[i].nei);
				node[i].nei = NULL;
			}
		}
		free(node);
		node = NULL;
	}
}

static bool write_distribution(const char *filename)
{
	char *buff;
	size_t buff_len = 1024;
	buff = (char *)calloc(buff_len, sizeof(char));

	FILE *f = fopen(filename, "w");
	if(f == NULL) {
		perror("fopen");
		return false;
	}

	FILE *fw = fopen("./weight.csv", "w");
	if(fw == NULL) {
		perror("fopen");
		return false;
	}

	int i;
	for(i=0; i<NODE_NUM; i++) {
		if(node[i].cur) {
			int j;
			for(j=0; j<node[i].cur; j++) {
				NEIGHBOR *nei = &(node[i].nei[j]);
				char *tmp;
				if(expect_false(buff_len < nei->num * 10))
					array_needsize(char, buff, buff_len, nei->num * 10, array_zero_init);
				tmp = buff;
				tmp[0] = 0;
				tmp += sprintf(tmp, "%d,%d,", i, nei->id);
				double_to_string(tmp, nei->delay_pdf, nei->num);
				fwrite(buff, sizeof(char), strlen(buff), f);
				fprintf(fw, "%d,%d,%lf\r\n", i, nei->id, nei->delay_average);
			}
		}
	}

	free(buff);
	buff = NULL;
	fclose(f);
	fclose(fw);

	return true;
}

bool cal_distribution(const char *inputF, const char *outputF)
{
	char cmd[512] = {0};

	fp = fopen(inputF, "r");
	if(fp == NULL) {
		perror("fopen");
		return false;
	}

	sprintf(cmd, "cut -d , -f 2 %s | sort | uniq | wc -l", inputF);
	NODE_NUM = atoi(cmd_system(cmd));	//num of nodes in the network
	node = (NODE *)calloc(NODE_NUM, sizeof(NODE));

	int i;
	for(i=1; i<=NODE_NUM; i++)
		get_node_info(i);

	bool res = write_distribution(outputF);

	exit_clean();
	return res;
}

MATRIX *distribution_to_matrix(bool *dense)
{
	MATRIX *m = (MATRIX *)calloc(NODE_NUM * NODE_NUM, sizeof(MATRIX));
	int i, j;
	double degree = 0;
	for(i=0; i<NODE_NUM; i++) {
		for(j=0; j<NODE_NUM; j++) {
			if(i == j)
				m_weight(m, i, j, NODE_NUM) = 0;
			else
				m_weight(m, i, j, NODE_NUM) = INF;

			m_parent(m, i, j, NODE_NUM) = -1;
		}

		for(j=0; j<node[i].cur; j++) {
			NEIGHBOR *n = &(node[i].nei[j]);
			m_weight(m, i, n->id, NODE_NUM) = n->delay_average;
			degree++;
		}
	}
	
	*dense = (degree/(double)NODE_NUM >= (double)NODE_NUM/2.) ? true : false;
	return m;
}

double *update_convolution(PATH *p, int *len)
{
	double *D[p->cur - 2];
	double *a1, *a2;
	int n1, n2;
	int c = 0;
	int s = p->path[0];
	NEIGHBOR key, *res;
	key.id = p->path[1];
	res = bsearch(&key, node[s].nei, node[s].cur, sizeof(NEIGHBOR), nei_cmp);
	a1 = res->delay_pdf;
	n1 = res->num;

	int i;
	for(i=1; i<p->cur - 1; i++) {
		s = p->path[i];
		key.id = p->path[i+1];
		res = bsearch(&key, node[s].nei, node[s].cur, sizeof(NEIGHBOR), nei_cmp);
		a2 = res->delay_pdf;
		n2 = res->num;

		D[c] = convolution(a1, n1, a2, n2);
		a1 = D[c];
		n1 = n1 + n2 - 1;
		c++;
	}
	
	//make cdf
	for(i=1; i<n1; i++)
		a1[i] += a1[i-1];

	//free memory
	for(i=0; i<c-1; i++)
		free(D[i]);

	if(len)
		*len = n1;
	return a1;
}

//double *convolution(double *a1, int n1, double *a2, int n2)
void do_convolution(FILE *f, PATH *p, PEER *n)
{
	static char conv_buff[10240] = {0};
	char *str = conv_buff;
	int len;

	str[0] = 0;
	str += sprintf(str, "%d,%d,%d,", p->path[0], p->path[p->cur - 1], p->cur - 1);

	n->stime = p->cur - 1;
	n->cdf = update_convolution(p, &len);

	double_to_string(str, n->cdf, len);
	fwrite(conv_buff, sizeof(char), strlen(conv_buff), f);
}

void write_record(MATRIX *G, bool type)
{
	FILE *fc = fopen("./ccdf.csv", "w");
	FILE *fpa = fopen("./path.csv", "w");
	int i, j;
	
	char *buff;
	size_t buff_len = 1024;
	buff = (char *)calloc(buff_len, sizeof(char));
	char pbuff[1024] = {0};

	int id, cur;
	for(i=0; i<NODE_NUM; i++) {
		id = i;
		cur = 0;
		for(j=0; j<NODE_NUM; j++) {
			m_path(G, i, j, NODE_NUM) = path(type, G, i, j, NODE_NUM);
			//do the convolution
			PATH *tmp =  m_path(G, i, j, NODE_NUM);
			if(tmp == NULL)
				continue;

			//write path
			char *str;
			str = pbuff;
			str[0] = 0;
			int_to_string(str, tmp->path, tmp->cur);
			fwrite(pbuff, sizeof(char), strlen(pbuff), fpa);

			cur++;
			p_ccdf[i] = (PEER *)realloc(p_ccdf[i], cur * sizeof(PEER));
			PEER *pt = p_ccdf[i] + cur - 1;
			pt->id = j;

			if(tmp->cur == 2) {
				NEIGHBOR key, *res;
				key.id = j;
				res = bsearch(&key, node[i].nei, node[i].cur, sizeof(NEIGHBOR), nei_cmp);

				if(expect_false(buff_len < res->num * 10))
					array_needsize(char, buff, buff_len, res->num * 10, array_zero_init);
				str = buff;
				str[0] = 0;
				str += sprintf(str, "%d,%d,%d,", i, res->id, 1);

				//make cdf
				int k;
				double *tmp_res = (double *)calloc(res->num, sizeof(double));
				tmp_res[0] = res->delay_pdf[0];
				for(k=1; k<res->num; k++)
					tmp_res[k] = res->delay_pdf[k] + tmp_res[k - 1];
				
				pt->stime = 1;
				pt->cdf = tmp_res;

				double_to_string(str, tmp_res, res->num);
				fwrite(buff, sizeof(char), strlen(buff), fc);
				continue;
			}
			//need to recal the delay_cdf
			do_convolution(fc, tmp, pt);
		}
		p_ccdf[id] = (PEER *)realloc(p_ccdf[id], (cur + 1) * sizeof(PEER));
		memset(p_ccdf[id] + cur, 0, sizeof(PEER));
	}

	free(buff);
	fclose(fc);
	fclose(fpa);
}

static inline double r1(void)
{
	return (double)rand() / (double)RAND_MAX;
}

PATH *path_merge(PATH *a, PATH *b)
{
	PATH *res = (PATH *)calloc(1, sizeof(PATH));
	res->num = a->cur + b->cur - 1;
	res->cur = res->num;
	res->path = (int *)calloc(res->cur, sizeof(int));

	int *tmp = res->path;
	memcpy(tmp, a->path, (a->cur - 1) * sizeof(int));
	tmp += a->cur - 1;
	memcpy(tmp, b->path, b->cur * sizeof(int));

	return res;
}

double cal_probability(double *cdf, int stime, int time)
{
	int rtime = stime * TSLOT;
	double *res = cdf;
	for(;;) {
		if(*res == 1)
			return *res;

		if(rtime >= time)
			return *res;

		rtime += TSLOT;
		res++;
	}
}

double get_probability(MATRIX *G, int s, int i, int j, int time)
{
	PATH *si;
	if(s == -1)
		si = NULL;
	else
		si = m_path(G, s, i, NODE_NUM);	

	PATH *ij = m_path(G, i, j, NODE_NUM);
	PATH *ji = m_path(G, j, i, NODE_NUM);
	if((s != -1 && si == NULL) || ij == NULL || ji == NULL) {
//		printf("no path %d-%d-%d\n", s, i, j);
		return 0;
	}

	PATH *sj;
	if(si == NULL)
		sj = ij;
	else
		sj = path_merge(si, ij);

	PATH *sji = path_merge(sj, ji);
	PATH *f = path_merge(sji, ij);

	double *new_cdf = update_convolution(f, NULL);
	double res = cal_probability(new_cdf, f->cur - 1, time);

	free(new_cdf);
	if(si) {
		free(sj->path);
		free(sj);
	}
	free(sji->path);
	free(sji);
	free(f->path);
	free(f);

	return res;
}

void write_cdf(MATRIX *G, int s, int time)
{
	FILE *f = fopen("probability.gams", "w");
	fprintf(f, "TABLE CP(I, J)	the cdf for path from source node to requestor J by candidate I\n\t");

	int i, j;
	for(i=0; i<NODE_NUM; i++)
		fprintf(f, "%d\t", i);
	fprintf(f, "\n");

	for(i=0; i<NODE_NUM; i++) {
		fprintf(f, "%d\t", i);
		for(j=0; j<NODE_NUM; j++) {
	//		printf("%d %d\n", i, j);
			fprintf(f, "%.3lf\t", get_probability(G, s, i, j, time));
		}
		fprintf(f, "\n");
	}

	fclose(f);
}

double cal_mrev(MATRIX *G, PINFO *n, int s, const char *x, int time)
{
	int i, j;
	double m, sum = 0, c = 0;
	for(j=0; j<NODE_NUM; j++) {
		if(n[j].interest == 0)
			continue;

		m = 1;
		for(i=0; i<NODE_NUM; i++) {
			if(x[i] == 0)
				continue;

			m = m * (1 - get_probability(G, s, i, j, time));
		}
		sum += (1 - m) * n[j].interest * PRICE;
	}

	for(i=0; i<NODE_NUM; i++)
		if(x[i] != 0)
			c += COST;

	m = sum - c;
	return (m < 0 ? 0 : m);
}

PEER *peer_search(peerlist p, int id)
{
	PEER *res = p;
	while(res->stime) {
		if(res->id == id)
			return res;
		res++;
	}
	return NULL;
}

PINFO *build_node_info(peerlist *p, int s, int time)
{
	PINFO *res = (PINFO *)calloc(NODE_NUM, sizeof(PINFO));

	int i;
	for(i=0; i<NODE_NUM; i++) {
		PEER *tmp = peer_search(p[s], i);
		if(tmp == NULL)
			continue;

		res[i].probability = cal_probability(tmp->cdf, tmp->stime, time);
		res[i].interest = r1();
	}

	return res;
}

void free_peerlist(peerlist *p, int num)
{
	int i;
	for(i=0; i<num; i++) {
		PEER *tmp = p[i];
		if(tmp == NULL)
			continue;

		while(tmp->stime) {
			free(tmp->cdf);
			tmp++;
		}
		free(p[i]);
	}
}

double obj_func(unsigned n, const double *x, double *grad, void *func_data)
{
	FUNC_DATA *d = (FUNC_DATA *)func_data;
	int s = d->snode, time = d->wtime;
	MATRIX *G = d->graph;
	PINFO *node = d->node;

	return cal_mrev(G, node, s, x, time);
}

double constraint_func1(unsigned n, const double *x, double *grad, void *data)
{
	int i;
	double sum = 0.;
	for(i=0; i<n; i++)
		sum += x[i];

	return sum - KTHRESH;
}

//should be 0
double constraint_func2(unsigned n, const double *x, double *grad, void *data)
{
	int i;
	double sum = 0.;
	for(i=0; i<n; i++)
		sum += x[i] * (x[i] - 1.);

	return sum;
}

void write_node_interest(PINFO *n)
{
	FILE *f = fopen("interest.gams", "w");
	int i = 0;
	fprintf(f, "\tR(J)	the probability of requestor needs the file\n\t/");
	fprintf(f, "\t%d\t%.3lf\n", i++, n[i].interest);
	for(i; i<NODE_NUM-1; i++) {
		fprintf(f, "\t\t%d\t%.3lf\n", i, n[i].interest);
	}
	fprintf(f, "\t\t%d\t%.3lf\t/;", i, n[i].interest);
	fclose(f);
}

void knapsack(dp_item **t, int num, int k, int *w, double *v, MATRIX *G, PINFO *n, int s, int time)
{
	dp_item *table = *t;
	int i, j;
	for(i=0; i<num; i++) {
		for(j=0; j<k; j++) {
		//	printf("%d:%d\n", i, j);
			dp_item *cur = matrix(table, i, j, k);

			if(i == 0 || j == 0) {
				cur->value = 0;
				continue;
			}

			dp_item *prev = matrix(table, i - 1, j, k);
			if(w[i - 1] <= j) {
				char tmp[NODE_NUM];
				memcpy(tmp, ((dp_item *)matrix(table, i - 1, j - w[i - 1], k))->selection, NODE_NUM * sizeof(char));
				tmp[cur->id] = 1;

				double rev = cal_mrev(G, n, s, tmp, time);
				if(rev > prev->value) {
					cur->value = rev;
					memcpy(cur->selection, tmp, NODE_NUM * sizeof(char));
				}
				else {
					cur->value = prev->value;
					memcpy(cur->selection, prev->selection, NODE_NUM * sizeof(char));
				}
			}
			else {
				cur->value = prev->value;
				memcpy(cur->selection, prev->selection, NODE_NUM * sizeof(char));
			}

		}
	}
}

void item_init(int **weight, double **value, int num, MATRIX *G, PINFO *n, int s, int time)
{
	int i;
	int *w = *weight;
	double *v = *value;
	char x[NODE_NUM];
	for(i=0; i<num; i++) {
		memset(x, 0, NODE_NUM * sizeof(char));
		w[i] = 1;
		
		x[i] = 1;
		v[i] = cal_mrev(G, n, s, x, time);
	}
}

void generate_adv(M_NODE *n, int time)
{
	int left = n->buff_len - n->buff_cur - NODE_NUM;
	if(left < 0)
		array_needsize(FDATA, n->buffer, n->buff_len, n->buff_len - left, array_zero_init);

	int i;
	for(i=0; i<NODE_NUM; i++) {
		FDATA *b = &(n->buffer[n->buff_cur++]);
		b->src = n->id;
		b->dest = i;
		b->type = FILE_ADV;
		b->stime = time;
	}
}

int find_next_hop(MATRIX *G, M_NODE *n, int dest)
{
	int res = -1, i;
	PATH *p = m_path(G, n->id, dest, NODE_NUM);
	if(p == NULL)
		return -1;
#ifdef FIXED_ROUTE
	//must follow the path we compute...sucks...
	if(n->neighbor[p->path[1]])
		return p->path[1];
	else
		return -1;
#else
	for(i=p->cur - 1; i>=0; i--) {
		int j = p->path[i];
		if(n->neighbor[j]) {
			res = j;
			break;
		}
	}
	return res;
#endif
}

static int check_data_duplicate(M_NODE *n, FDATA data)
{
	int i, res = 0;
	FDATA *tmp;
	for(i=0; i<n->buff_cur; i++) {
		tmp = &(n->buffer[i]);
		if(tmp->type == data.type && tmp->src == data.src && tmp->dest == data.dest) {
			res = 1;
			break;
		}
	}
	return res;	/*0 - no duplicate data; 1 - duplicate data*/
}

static void send_data(M_NODE *n, FDATA *b)
{
	if(check_data_duplicate(n, *b) == 0) {
		if(expect_false(n->buff_cur == n->buff_len))
			array_needsize(FDATA, n->buffer, n->buff_len, n->buff_len + 1, array_zero_init);

		memcpy(&(n->buffer[n->buff_cur++]), b, sizeof(FDATA));
	}
}

void handle_node(M_NODE *list[], int num, M_NODE *node, int stime, int rtime, MATRIX *G)
{
	int i, j;
	M_NODE *n;
	FDATA *b;
/*
 * 1. handle possible direct transfer
 * 2. check buffer for multiple hop data receive
 */
	//step 1. handle all possible FILE_TRANS
	for(i=0; i<num; i++) {
		n = list[i];
		//handle source node at first
		if(n->source) {
			for(j=0; j<NODE_NUM; j++) {
				if(n->neighbor[j]) {
#ifdef FIXED_ROUTE
					if(j != find_next_hop(G, n, j))
						continue;
#endif
					//if it's candidate
					if(node[j].candidate) {
						//if the candidate hasn't received the shared file yet, distribute the file to the candidate
						if(node[j].have_file == false) {
							node[j].have_file = true;
							sim_delay += rtime - stime;
							sim_delivery++;
							sim_rev += node[j].interest * PRICE - COST;
							_dprintf("SRC_DISTRIBUTION@%d: #%d - #%d\n", rtime, n->id, j);
							//generate adv
							generate_adv(&(node[j]), rtime);
						}	
					}
					//if it's requestor and does not have the file
					else if(node[j].have_file == false) {
						node[j].have_file = true;
						sim_delivery++;
						sim_delay += rtime - stime;
						_dprintf("SRC_TRANS@%d: #%d - #%d\n", rtime, n->id, j);
						sim_rev += node[j].interest * PRICE;
					}
				}
			}
		}
	}
//	debug(num);
	//handle if candidate has file to transfer at second
	for(i=0; i<num; i++) {
		n = list[i];
		if(!n->source && n->candidate && n->have_file) {
			//search neighbor
			for(j=0; j<NODE_NUM; j++) {
				if(n->neighbor[j] && 
					node[j].candidate == false &&
					node[j].have_file == false) {
#ifdef FIXED_ROUTE
					if(j != find_next_hop(G, n, j))
						continue;
#endif
					node[j].have_file = true;
					sim_delivery++;
					sim_delay += rtime - stime;
					_dprintf("CAN_TRANS@%d: #%d - #%d\n", rtime, n->id, j);
					sim_rev += node[j].interest * PRICE;
				}
			}
		}
	}
//	debug(num);
	//check if any node has file in buffer to transfer at last
	for(i=0; i<num; i++) {
		n = list[i];
		for(j=0; j<n->buff_cur; j++) {
			b = &(n->buffer[j]);
			if(b->type == FILE_TRANS && 
				n->neighbor[b->dest] &&
				node[b->dest].have_file == false) {
#ifdef FIXED_ROUTE
				if(b->dest != find_next_hop(G, n, b->dest))
					continue;
#endif
				node[b->dest].have_file = true;
				sim_delivery++;
				sim_delay += rtime - stime;
				_dprintf("RELAY_TRANS@%d: #%d - #%d\n", rtime, b->src, b->dest);
				sim_rev += node[b->dest].interest * PRICE;
				
				remove_data(b, n);
			}
		}
	}
//	debug(num);
	//step 2: handle FILE_DISTRIBUTION to generate ADV
	for(i=0; i<num; i++) {
		n = list[i];
		for(j=0; j<n->buff_cur; j++) {
			b = &(n->buffer[j]);
			if(b->type == FILE_DISTRIBUTION &&
				n->neighbor[b->dest] &&
				node[b->dest].candidate &&
				node[b->dest].have_file == false) {
#ifdef FIXED_ROUTE
				if(b->dest != find_next_hop(G, n, b->dest))
					continue;
#endif
				node[b->dest].have_file = true;
				sim_delivery++;
				sim_delay += rtime - stime;
				_dprintf("RELAY_DISTRIBUTION@%d: #%d - #%d\n", rtime, b->src, b->dest);
				sim_rev += node[b->dest].interest * PRICE - COST;
				generate_adv(&(node[b->dest]), rtime);

				remove_data(b, n);
			}
		}
	}

//	debug(num);
	//step 3: handle ADV to generate REQ
	for(i=0; i<num; i++) {
		n = list[i];
		for(j=0; j<n->buff_cur; j++) {
			b = &(n->buffer[j]);
			if(b->type == FILE_ADV &&
				n->neighbor[b->dest] &&
				node[b->dest].candidate == false &&
				node[b->dest].source == false &&
				node[b->dest].have_file == false) {
				M_NODE *cn = &(node[b->dest]);
				cn->candidate_list[b->src] = 1;
				int c, k, flag = 0;
#ifdef FIXED_ROUTE
				if(b->dest != find_next_hop(G, n, b->dest))
					continue;
#endif
				//!!!!! could always generate REQ when there are neighbors around...
				node[b->dest].candidate_list[b->src] = 1;
				_dprintf("recv FILE_ADV@%d: #%d - #%d by: %d\n", rtime, b->src, b->dest, n->id);
				//check if the REQ we are going to creat is duplicated
				FDATA *dbuff;
#ifdef IMPROVE_SCHEME
//!!!!! could always generate REQ when there are neighbors around...
				for(c=0; c<NODE_NUM; c++) {
					flag = 0;
					if(cn->candidate_list[c] == 0)
						continue;
#else
					c = b->src;
					_dprintf("recv FILE_ADV@%d: #%d - #%d by: %d\n", rtime, b->src, b->dest, n->id);
#endif
					//check if the REQ we are going to creat is duplicated
					FDATA tmp_data = {.src = cn->id, .dest = c, .type = FILE_REQ};
					flag = check_data_duplicate(cn, tmp_data);
					if(flag == 0) {
						if(expect_false(cn->buff_cur == cn->buff_len))
							array_needsize(FDATA, cn->buffer, cn->buff_len, cn->buff_len + 1, array_zero_init);
						dbuff = &(cn->buffer[cn->buff_cur++]);
						dbuff->src = cn->id;
						dbuff->dest = c;
						dbuff->type = FILE_REQ;
						dbuff->stime = rtime;
					}
#ifdef IMPROVE_SCHEME
				}
#endif
#if 0
				c = b->src;
				for(k=0; k<cn->buff_cur; k++) {
					dbuff = &(cn->buffer[k]);
					if(dbuff->type == FILE_REQ && 
						dbuff->src == cn->id && 
						dbuff->dest == c) {
						flag = 1;
						break;
					}
				}
				if(flag == 0) {
					if(expect_false(cn->buff_cur == cn->buff_len))
						array_needsize(FDATA, cn->buffer, cn->buff_len, cn->buff_len + 1, array_zero_init);
					dbuff = &(cn->buffer[cn->buff_cur++]);
					dbuff->src = b->dest;
					dbuff->dest = c;
					dbuff->type = FILE_REQ;
					dbuff->stime = rtime;
				}
#endif
				remove_data(b, n);
			}
		}
	}

//	debug(num);
	//step 4: handle REQ to generate TRANS
	for(i=0; i<num; i++) {
		n = list[i];
		for(j=0; j<n->buff_cur; j++) {
			b = &(n->buffer[j]);
#ifdef IMPROVE_SCHEME
			if(b->type == FILE_REQ) {
				if(n->have_file && b->src == n->id) {
					remove_data(b, n);
					continue;
				}

				//!!!!! could be recv by different candidate...
				int c;
				for(c=0; c<NODE_NUM; c++) {
					if(n->neighbor[c] &&
						node[c].candidate &&
						node[c].have_file == true) {
							FDATA tmp_data = {.src = c, .dest = b->src, .type = FILE_TRANS};
							if(check_data_duplicate(&(node[c]), tmp_data))
								continue;

							FDATA *dbuff;
							if(expect_false(node[c].buff_cur == node[c].buff_len))
								array_needsize(FDATA, node[c].buffer, node[c].buff_len, node[c].buff_len + 1, array_zero_init);
							dbuff = &(node[c].buffer[node[c].buff_cur++]);
							dbuff->src = c;
							dbuff->dest = b->src;
							//must be sent back by multiple hops
							dbuff->type = FILE_TRANS;
							dbuff->stime = rtime;
						}
				}
				remove_data(b, n);
			}
#endif
			if(b->type == FILE_REQ &&
				n->neighbor[b->dest] &&
				node[b->dest].candidate &&
				node[b->dest].have_file == true) {
#ifdef FIXED_ROUTE
				if(b->dest != find_next_hop(G, n, b->dest))
					continue;
#endif
				if(n->have_file && b->src == n->id) {
					remove_data(b, n);
				}
				else {
					FDATA tmp_data = {.src = b->dest, .dest = b->src, .type = FILE_TRANS};
					if(check_data_duplicate(&node[b->dest], tmp_data) == 0) {
						_dprintf("recv FILE_REQ@%d: #%d - #%d by: %d\n", rtime, b->src, b->dest, n->id);
						FDATA *dbuff;
						if(expect_false(node[b->dest].buff_cur == node[b->dest].buff_len))
							array_needsize(FDATA, node[b->dest].buffer, node[b->dest].buff_len, node[b->dest].buff_len + 1, array_zero_init);
						dbuff = &(node[b->dest].buffer[node[b->dest].buff_cur++]);
						dbuff->src = b->dest;
						dbuff->dest = b->src;
						//must be sent back by multiple hops
						dbuff->type = FILE_TRANS;
						dbuff->stime = rtime;
					}
					remove_data(b, n);
				}
			}
		}
	}

//	debug(num);
	//step 5: handle all the data left in buffer for multiple-hop transfer
	int next_hop;
	for(i=0; i<num; i++) {
		n = list[i];	
		for(j=0; j<n->buff_cur; j++) {
			b = &(n->buffer[j]);
			if(b->type == FILE_DISTRIBUTION) {
				if(node[b->dest].have_file == false) {
					//find possible next hop
					next_hop = find_next_hop(G, n, b->dest);
					if(next_hop != -1) {
						send_data(&node[next_hop], b);
						_dprintf("send FILE_DISTRIBUTION@%d - src:%d, dest:%d, from:%d, to:%d\n", rtime, b->src, b->dest, n->id, next_hop);
						remove_data(b, n);
					}
				}
				else {
					remove_data(b, n);
				}
			}
			else if(b->type == FILE_ADV) {
				if(node[b->dest].have_file == false) {
					next_hop = find_next_hop(G, n, b->dest);
					if(next_hop != -1) {
						_dprintf("send FILE_ADV@%d - src:%d, dest:%d, from:%d, to:%d\n", rtime, b->src, b->dest, n->id, next_hop);
						send_data(&node[next_hop], b);
						remove_data(b, n);
					}
				}
				else {
					remove_data(b, n);
				}
			}
			else if(b->type == FILE_REQ) {
				if(n->have_file == true && b->src == n->id) {
					remove_data(b, n);
				}
				else {
					next_hop = find_next_hop(G, n, b->dest);
					if(next_hop != -1) {
						_dprintf("send FILE_REQ@%d - src:%d, dest:%d, from:%d, to:%d\n", rtime, b->src, b->dest, n->id, next_hop);
						send_data(&node[next_hop], b);
						remove_data(b, n);
					}
				}
			}
			else if(b->type == FILE_TRANS) {
				if(node[b->dest].have_file == false) {
					next_hop = find_next_hop(G, n, b->dest);
					if(next_hop != -1) {
						_dprintf("send FILE_TRANS@%d - src:%d, dest:%d, from:%d, to:%d\n", rtime, b->src, b->dest, n->id, next_hop);
						send_data(&node[next_hop], b);
						remove_data(b, n);
					}
				}
				else {
					remove_data(b, n);
				}
			}
		}
	}
//	debug(num);
}

void node_communication(char *neighbor, M_NODE *node, int stime, int rtime, MATRIX *G)
{
	M_NODE *list[NODE_NUM];
	int i, cur = 0;
	for(i=0; i<NODE_NUM; i++) {
		if(neighbor[i]) {
			memcpy(node[i].neighbor, neighbor, sizeof(char) * NODE_NUM);
			node[i].neighbor[i] = 0;	//clear self bit
			list[cur++] = &(node[i]);
		}
	}

	handle_node(list, cur, node, stime, rtime, G);
}

void simulation_loop(int source_node, int stime, long wtime, char *candidate, PINFO *info, MATRIX *G)
{
/*
	start loop...
*/
	//init mobilie node at first...
	M_NODE *node = (M_NODE *)calloc(NODE_NUM, sizeof(M_NODE));
	int i;
	for(i=0; i<NODE_NUM; i++) {
		node[i].id = i;
		node[i].neighbor = (char *)calloc(NODE_NUM, sizeof(char));
		node[i].candidate_list = (char *)calloc(NODE_NUM, sizeof(char));

		if(i == source_node) {	//init source node
			node[i].candidate = true;
			node[i].have_file = true;
			node[i].source = true;

			array_needsize(FDATA, node[i].buffer, node[i].buff_len, NODE_NUM, array_zero_init);
			//generate data
			int j;
			for(j=0; j<NODE_NUM; j++) {
				if(expect_false(node[i].buff_cur == node[i].buff_len))
					array_needsize(FDATA, node[i].buffer, node[i].buff_len, node[i].buff_len + 1, array_zero_init);

				FDATA *t = &(node[i].buffer[node[i].buff_cur++]);
				//generate file distribution data
				if(candidate[j]) {
					t->src = i;
					t->dest = j;
					t->stime = stime;
					t->type = FILE_DISTRIBUTION;
				}
				//generate the adv data
				else {
					t->src = i;
					t->dest = j;
					t->stime = stime;
					t->type = FILE_ADV;
				}
			}
		}
		else {
			if(candidate[i])
				node[i].candidate = true;
			node[i].interest = info[i].interest;
		}
	}

	//init done, start loop
	char neighbor[NODE_NUM];
	memset(neighbor, 0, sizeof(char) * NODE_NUM);
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	FILE *f = fopen("./mobility.csv", "r");
	int rtime = stime;

	while((read = getline(&line, &len, f)) != -1) {
		if(wtime != -1 && rtime - stime > wtime)
			break;

		int n1, n2, time;
		sscanf(line, "%d,%d,%d", &n1, &n2, &time);
		n1--; n2--;
		if(time == rtime) {
			neighbor[n1] = 1;
			neighbor[n2] = 1;
			continue;
		}

		if(time > rtime) {
			//handle previous time slot record
			node_communication(neighbor, node, stime, rtime, G);

			//reset record
			rtime = time;
			memset(neighbor, 0, sizeof(char) * NODE_NUM);
			neighbor[n1] = 1;
			neighbor[n2] = 1;
		}
	}

	//free M_NODE
	for(i=0; i<NODE_NUM; i++) {
		if(node[i].neighbor)
			free(node[i].neighbor);
		if(node[i].candidate_list)
			free(node[i].candidate_list);
		if(node[i].buffer)
			free(node[i].buffer);
	}
	free(node);
	free(line);
	fclose(f);
}

int get_meetingEvent(int node, int stime, long wtime)
{
	FILE *f = fopen("./mobility.csv", "r");
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	int res = 0;

	while((read = getline(&line, &len, f)) != -1) {
		int i, j, time;
		sscanf(line, "%d,%d,%d", &i, &j, &time);
		i--; j--;
		if(time > stime + wtime) 
			break;
		if(i == node || j == node)
			res++;
	}

	fclose(f);
	free(line);
	return res;
}

//return the candidate node id
int get_max_obRev(int node, int stime, int wtime, int k, int *best_candidate, PINFO *n, MATRIX *G)
{
	FILE *f = fopen("./mobility.csv", "r");
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	int cnt = 0;
	char x[NODE_NUM];
	double ob_max = 0, res = 0;
	int s_can = -1;
	char *ob_test, *can_test;
	ob_test = (char *)calloc(NODE_NUM, sizeof(char));
	can_test = (char *)calloc(NODE_NUM, sizeof(char));
	long wtime_s = wtime * 60;

	while((read = getline(&line, &len, f)) != -1) {
		int i, j, time;
		sscanf(line, "%d,%d,%d", &i, &j, &time);
		i--; j--;
		if(time > stime + wtime_s)
			break;
		if(i == node || j == node) {
			cnt++;

			int neighbor = (node - i == 0) ? j : i;

			if(cnt <= k && ob_test[neighbor])
				continue;

			if(cnt > k && can_test[neighbor])
				continue;

			memset(x, 0, sizeof(char) * NODE_NUM);
			x[neighbor] = 1;
			double rev = cal_mrev(G, n, -1, x, wtime * 3/*actual file transfer time window*/);	//we do not need to consider the path from souce node to candidate
			if(cnt <= k) {
			//within the observing time window
				ob_test[neighbor] = 1;
				ob_max = max(ob_max, rev);
			}
			else {	//start the real candidate selection
				can_test[neighbor] = 1;
				if(rev >= ob_max) {
					if(s_can == -1) {	//first node we meet with higher value than the ob_max after the first k meeting event
						s_can = neighbor;
						printf("rev: %lf\n", rev);
					}
					else {
						*best_candidate = neighbor;
						res = max(res, rev);
					}
				}
			}
		}
	}

	printf("max rev: %lf\n", res);

	fclose(f);
	free(line);
	free(ob_test);
	free(can_test);
	return s_can;
}

void direct_communication(char *neighbor, M_NODE *node, int stime, int rtime, MATRIX *G, PINFO *info, int source_node, int wtime, double threshold)
{
	static bool get_can = false;
	bool trans = false;
	int i;
	char x[NODE_NUM];
	for(i=0; i<NODE_NUM; i++) {
		if(neighbor[i]) {
			//if it is the candidate and it has the file, then we could transfer the file to neighbors
			if(node[i].candidate && node[i].have_file)
				trans = true;

			if(rev_test[i])
				continue;

			rev_test[i] = 1;
			memset(x, 0, sizeof(char) * NODE_NUM);
			x[i] = 1;
			double tmp = cal_mrev(G, info, source_node, x, wtime);
			//if it's equal or larger than threshold, we do not select the candidate yet, and it is not candidate yet, then we select the node as the candidate.
			if(tmp >= threshold && 
				get_can == false &&
				node[i].candidate == false) {
				node[i].candidate = true;
				get_can = true;
			}
			//record the best choice...
			if(tmp > best_rev) {
				best_rev = tmp;
				_candidate = i;
			}
		}
	}
	if(trans) {
		for(i=0; i<NODE_NUM; i++) {
			if(neighbor[i] && node[i].have_file == false) {
				node[i].have_file = true;
				sim_rev += node[i].interest * PRICE;
				if(!node[i].source && node[i].candidate)
					sim_rev -= COST;
				sim_delivery++;
				sim_delay += rtime - stime;
			}
		}
	}

}

void simulation_start(int source_node, int stime, int ob_time, long wtime, PINFO* n, MATRIX *G, double threshold)
{
	M_NODE *node = (M_NODE *)calloc(NODE_NUM, sizeof(M_NODE));
	int i;
	for(i=0; i<NODE_NUM; i++) {
		node[i].id = i;
		if(i == source_node) {
			node[i].source = true;
			node[i].candidate = true;
			node[i].have_file = true;
			node[i].interest = 0;
		}
		else
			node[i].interest = n[i].interest;
	}

	FILE *f = fopen("./mobility.csv", "r");
	char *line = NULL;
	size_t len = 0;
	ssize_t read;

	int ori_wtime = wtime;
	int ori_stime = stime;
	wtime = wtime - (ob_time - stime);
	stime = ob_time;
	int rtime = stime;
	char neighbor[NODE_NUM];
	memset(neighbor, 0, sizeof(char) * NODE_NUM);
	rev_test = (char *)calloc(NODE_NUM, sizeof(char));

	while((read = getline(&line, &len, f)) != -1) {
		if(rtime - stime > wtime)
			break;

		int n1, n2, time;
		sscanf(line, "%d,%d,%d", &n1, &n2, &time);
		if(time < stime)
			continue;	//not start yet

		n1--; n2--;
		if(time == rtime) {
			neighbor[n1] = 1;
			neighbor[n2] = 1;
			continue;
		}

		if(time > rtime) {
			//handle previous time slot record
			//node_communication(neighbor, node, stime, rtime, G);
			direct_communication(neighbor, node, ori_stime, rtime, G, n, source_node, ori_wtime, threshold);

			//reset record
			rtime = time;
			memset(neighbor, 0, sizeof(char) * NODE_NUM);
			neighbor[n1] = 1;
			neighbor[n2] = 1;
		}
	}

	for(i=0; i<NODE_NUM; i++) {
		if(!node[i].source && node[i].candidate) {
			printf("the selected candidate: %d\n", i);
			break;
		}
	}

	free(node);
	free(line);
	free(rev_test);
	fclose(f);
}

void distributed_simulation(int source_node, int stime, int wtime, PINFO *n, MATRIX * G)
{
#define DRATIO	0.4
	int best_candidate = -1;
	sim_delivery = 0;
	sim_delay = 0;
	sim_rev = 0;
#ifdef SINGLE_SELECT
	int total_events = get_meetingEvent(source_node, stime, (wtime/3)*60/*the candidate selection time window*/);
	int ob_events = total_events * DRATIO;
	int candidate = get_max_obRev(source_node, stime, wtime/3, ob_events, &best_candidate, n, G);
	printf("selected candidate: %d, best candidate: %d\n", candidate, best_candidate);
	//start real file distributioin. Start time: stime+wtime; file validate time: wtime
//	simulation_start(source_node, stime, ob_time, wtime, n, G, ob_rev);
#else
#endif
#undef DRATIO
}


#define DP_OPT
int main(int argc, char *argv[])
{
	cal_distribution(argv[1], "./pdf.csv");
	int i;

	bool flag;
	MATRIX *G = distribution_to_matrix(&flag);
	if(flag)
		folyd_warshall(G, NODE_NUM);
	else {
		double dist[NODE_NUM];
		for(i=0; i<NODE_NUM; i++)
			dijkstra(G, dist, i, NODE_NUM);
		
	}

	p_ccdf = (peerlist *)calloc(NODE_NUM, sizeof(peerlist));

	write_record(G, flag);

	srand(_seed);
	int source_node = 0, wtime = 1500;	//time window is xx min
	PINFO *ni = build_node_info(p_ccdf, source_node, wtime);

	write_node_interest(ni);
	write_cdf(G, source_node, wtime);

	char x[NODE_NUM];
	memset(x, 0, NODE_NUM * sizeof(char));
/*
	x[8] = 1;
	x[9] = 1;
	x[12] = 1;
	x[33] = 1;
	x[58] = 1;
*/
	x[31] = 1;
	x[2] = 1;
	x[12] = 1;
	x[25] = 1;
	x[43] = 1;
#ifdef USE_NLOPT
	FUNC_DATA fdata;
	fdata.snode = 0;
	fdata.wtime = wtime;
	fdata.graph = G;
	fdata.node = ni;

	double maxf;
	nlopt_opt opt;
	opt = nlopt_create(NLOPT_LN_COBYLA, NODE_NUM);
	nlopt_set_lower_bounds1(opt, 0);
	nlopt_set_upper_bounds1(opt, 1);
	nlopt_set_max_objective(opt, obj_func, &fdata);
	nlopt_add_inequality_constraint(opt, constraint_func1, NULL, 0.1);
	nlopt_add_equality_constraint(opt, constraint_func2, NULL, 0.1);

	nlopt_set_xtol_rel(opt, 0.1);
	if(nlopt_optimize(opt, x, &maxf) < 0)
		printf("nlopt failed\n");
	else {
		printf("maxv: %lf\n", maxf);
		for(i=0; i<NODE_NUM; i++)
			if(x[i] == 1)
				printf("#%d ", i);
		printf("\n");
	}
	nlopt_destroy(opt);
#else
	//cal_mrev() func test...
	double rev = cal_mrev(G, ni, source_node, x, wtime);
	printf("\n=============OPT RESULTS===================\n\n");
	printf("random rev: %lf\n", rev);
	printf("random candidates:\t");
	for(i=0; i<NODE_NUM; i++) {
		if(x[i] != 0)
			printf("#%d\t", i);
	}
	printf("\n\n");
#endif	

#ifdef DP_OPT
	int max_weight = 6;	//the total num of nodes we could choose is (max_weight - 1)
	int *i_weight = (int *)calloc(NODE_NUM, sizeof(int));
	double *i_value = (double *)calloc(NODE_NUM, sizeof(double));
	item_init(&i_weight, &i_value, NODE_NUM, G, ni, source_node, wtime);

	dp_item *dp = (dp_item *)calloc(max_weight * NODE_NUM, sizeof(dp_item));
	for(i=0; i<max_weight * NODE_NUM; i++) {
		dp[i].selection = (char *)calloc(NODE_NUM, sizeof(char));
		dp[i].id = i/max_weight;
	}

	knapsack(&dp, NODE_NUM, max_weight, i_weight, i_value, G, ni, source_node, wtime);
	dp_item final = dp[i-1];
	printf("max rev: %lf\n", final.value);
	printf("candidates:\t");
	for(i=0; i<NODE_NUM; i++) {
		if(final.selection[i] != 0)
			printf("#%d\t", i);
	}
	printf("\n");
#endif
/////////////////////////SIMULATION///////////////////////////////////////////
	printf("\n============SIMULATION RESULTS===================\n\n");

	char s_cmd[] = "head -1 ./graph.csv | cut -d , -f 3";
	int stime = atoi(cmd_system(s_cmd));
	simulation_loop(source_node, stime, wtime * 60, final.selection, ni, G);
	printf("sim revenue: %lf\n", sim_rev);
	printf("total sharing: %d\n", sim_delivery);
	printf("average delay: %lf\n", (double)sim_delay/(double)sim_delivery);

	printf("\n=====================================================\n\n");

	distributed_simulation(source_node, stime, wtime, ni, G);
//	printf("the best candidate: %d\n", best_candidate);
	printf("sim revenue: %lf\n", sim_rev);
	printf("total sharing: %d\n", sim_delivery);
	printf("average delay: %lf\n", (double)sim_delay/(double)sim_delivery);

	printf("\n===============DONE===========================\n\n");
/////////////////////////CLEAN UP//////////////////////////////////////////////

#ifdef DP_OPT
	free(i_weight);
	free(i_value);
	for(i=0; i<max_weight * NODE_NUM; i++) {
		if(dp[i].selection)
			free(dp[i].selection);
	}
	free(dp);
#endif
	for(i=0; i<NODE_NUM * NODE_NUM; i++) {
		PATH *tmp = G[i].path;
		if(tmp && tmp->path) {
			free(tmp->path);
			free(tmp);
		}
	}
	free(G);

	node_free();
	free(ni);
	free_peerlist(p_ccdf, NODE_NUM);
	return 0;
}

