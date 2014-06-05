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

#include "distribution.h"

static const double INF = DBL_MAX/2 - 1;
static const unsigned int _seed = INT_MAX - 1;

static char *line = NULL;
static NODE *node;
static unsigned int NODE_NUM = 0;
static FILE *fp;
static unsigned int *delay_t = NULL;
static int dnum_t = 0;
static double *pdf = NULL;
static int pdf_len = 0;
static int pdf_cur = 0;
static peerlist *p_ccdf = NULL;

#define TSLOT	60
#define KTHRESH	6

typedef struct func_data_t {
	int snode, wtime;
	MATRIX *graph;
	PINFO *node;
} FUNC_DATA;

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
#undef NEIGHBOR_THRESHOLD
}

static void exit_clean(void)
{
	if(delay_t) {
		free(delay_t);
		delay_t = NULL;
	}

	if(line) {
		free(line);
		line = NULL;
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
	PATH *si = m_path(G, s, i, NODE_NUM);	
	PATH *ij = m_path(G, i, j, NODE_NUM);
	if(si == NULL || ij == NULL) {
//		printf("no path %d-%d-%d\n", s, i, j);
		return 0;
	}

	PATH *sj = path_merge(si, ij);

	double *new_cdf = update_convolution(sj, NULL);
	double res = cal_probability(new_cdf, sj->cur - 1, time);

	free(new_cdf);
	free(sj->path);
	free(sj);
	return res;
}

double cal_mrev(MATRIX *G, PINFO *n, int s, const double *x, int time)
{
#define PRICE 50
	int i, j;
	double m, sum = 0;
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

	return sum;
#undef PRICE
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
	int i, sum = 0;
	for(i=0; i<n; i++)
		sum += x[i];

	return sum - KTHRESH;
}

//should be 0
double constraint_func2(unsigned n, const double *x, double *grad, void *data)
{
	int i, sum = 0;
	for(i=0; i<n; i++)
		sum += x[i] * (x[i] - 1);

	return sum;
}

int main(int argc, char *argv[])
{
	cal_distribution(argv[1], "./pdf.csv");

	bool flag;
	MATRIX *G = distribution_to_matrix(&flag);
	if(flag)
		folyd_warshall(G, NODE_NUM);
	else {
		double dist[NODE_NUM];
		int i;
		for(i=0; i<NODE_NUM; i++)
			dijkstra(G, dist, i, NODE_NUM);
		
	}

	p_ccdf = (peerlist *)calloc(NODE_NUM, sizeof(peerlist));

	write_record(G, flag);

	srand(_seed);
	int source_node = 0, wtime = 500;
	PINFO *ni = build_node_info(p_ccdf, source_node, wtime);

	double x[NODE_NUM];
	memset(x, 0, NODE_NUM * sizeof(double));
	int i;
	for(i=0; i<KTHRESH; i++)
		x[i] = 1;

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
	nlopt_add_inequality_constraint(opt, constraint_func1, NULL, 1e-8);
	nlopt_add_equality_constraint(opt, constraint_func2, NULL, 1e-8);

	nlopt_set_xtol_rel(opt, 1e-4);
#if 1
	if(nlopt_optimize(opt, x, &maxf) < 0)
		printf("nlopt failed\n");
	else
		printf("maxv: %lf\n", maxf);
#else
	double rev = cal_mrev(G, ni, source_node, x, wtime);
	printf("%lf\n", rev);
#endif	
	nlopt_destroy(opt);
	node_free();
	free(ni);
	free_peerlist(p_ccdf, NODE_NUM);
	return 0;
}

