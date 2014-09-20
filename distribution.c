#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include "common.h"
#include "shortest_path.h"
#include "convolution.h"
#include "simulation.h"
#include "info_collect.h"

#include "distribution.h"

//#define FIXED_ROUTE
//#define SINGLE_SELECT
#define DP_OPT
#define DISTRI_SIM
#define CENTRA_SIM
//#define INTR_TEST
//#define _DEBUG

#ifdef INTR_TEST
#include <gsl/gsl_rng.h>
#endif
#ifdef USE_NLOPT
#include <nlopt.h>
#endif

static const double INF = DBL_MAX/2 - 1;
static const unsigned int _seed = INT_MAX - 1;

NODE *node;
unsigned int NODE_NUM = 0;
static FILE *fp;
static unsigned int *delay_t = NULL;
static size_t dnum_t = 0;
static double *pdf = NULL;
static size_t pdf_len = 0;
static size_t pdf_cur = 0;
static peerlist *p_ccdf = NULL;
static double sim_rev = 0;
static long sim_delay = 0;
static int sim_delivery = 0;
static int _candidate = -1;
static double best_rev = 0;
static char *rev_test;
static double *g_interest = NULL;
static FILE *f_src;
static FILE *f_can;
static FILE *f_dst;
static FILE *f_log;
static FILE *f_scan;
static CAN_NODE c_can_log;
static DST_NODE c_dst_log;
static CAN_NODE d_can_log;
static DST_NODE d_dst_log;
static int c_runtime;
static int d_runtime;
static int sim_type;
static int sim_no;

//#define NEIGHBOR_THRESHOLD		10000
//#define TSLOT	1800	//seconds
#define NEIGHBOR_THRESHOLD		100
#define TSLOT	300
#define SIM_ROUND	10

//#define KTHRESH	6

//default variables value
//int TIME_WINDOW = 10080;		//mins
int TIME_WINDOW = 600;
int CAN_NUM = 4;
int PRICE = 50;
int COST = 20;
int OB_WINDOW = 7;
double DRATIO = 0.4;

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

static int unsigned_cmp(const void *n1, const void *n2)
{
	return (*(unsigned int *)n1 - *(unsigned int *)n2);
}

static int int_cmp(const void *n1, const void *n2)
{
	return (*(int *)n1 - *(int *)n2);
}

static int nei_cmp(const void *n1, const void *n2)
{
	return (((NEIGHBOR *)n1)->id - ((NEIGHBOR *)n2)->id);
}

static double cal_pdf(unsigned int *array, size_t num)
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

static void neighbor_wb(unsigned int *delay/*neighbor delay distribution*/, size_t num/*num of inter contact time record*/, int nei_id/*neighbor id*/, NODE *n/*node*/)
{
	//sort for future use at first...
	qsort(delay, num, sizeof(unsigned int), unsigned_cmp);

	if(expect_false(n->cur == n->num))
		array_needsize(NEIGHBOR, n->nei, n->num, n->num + 1, array_zero_init);

	NEIGHBOR *p = &(n->nei[n->cur++]);
	p->id = nei_id;
	p->cn = num;
	p->delay_average = cal_pdf(delay, num);
	p->num = pdf_cur;
	p->delay_pdf = (double *)calloc(pdf_cur, sizeof(double));
	memcpy(p->delay_pdf, pdf, pdf_cur * sizeof(double));
}

static void get_node_info(unsigned int id)
{

	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	unsigned int p_neighbor = 0;
	size_t cur_t = 0;

	NODE *n = &(node[id]);
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
				neighbor_wb(delay_t, cur_t - 2, p_neighbor, n);

			break;
		}

		//the infor of neighbor p_neighbor has been all handled
		if(p_neighbor != neighbor) {
			if(cur_t > NEIGHBOR_THRESHOLD)
				neighbor_wb(delay_t/*neighbor delay distribution*/, cur_t - 2/*num of inter contact time record*/, p_neighbor/*neighbor id*/, n/*node*/);

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
		pdf_len = 0;
		pdf_cur = 0;
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

	FILE *fc = fopen("./contact.csv", "w");
	if(fc == NULL) {
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
				if(nei->cn)
					fprintf(fc, "%d,%d,%ld\r\n", i, nei->id, nei->cn);
			}
		}
	}

	free(buff);
	buff = NULL;
	fclose(f);
	fclose(fw);
	fclose(fc);

	return true;
}

bool cal_distribution(const char *inputF)
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
	for(i=0; i<NODE_NUM; i++)
		get_node_info(i);

	//bool res = write_distribution(outputF);

	exit_clean();
	return true;
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

double *update_convolution(PATH *p, int time, int *len)
{
	double *D[p->cur - 2];
	double *a1, *a2;
	int n1, n2, max_n;
	int c = 0;
	int s = p->path[0];
	NEIGHBOR key, *res;
	key.id = p->path[1];
	res = bsearch(&key, node[s].nei, node[s].cur, sizeof(NEIGHBOR), nei_cmp);
	a1 = res->delay_pdf;
	n1 = res->num;
	max_n = time /TSLOT / 2;

	int i;
	for(i=1; i<p->cur - 1; i++) {
		s = p->path[i];
		key.id = p->path[i+1];
		res = bsearch(&key, node[s].nei, node[s].cur, sizeof(NEIGHBOR), nei_cmp);
		a2 = res->delay_pdf;
		n2 = res->num;

		if(time != -1) {
			n1 = n1 > max_n ? max_n : n1;
			n2 = n2 > max_n ? max_n : n2;
		}
		
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
	static char conv_buff[102400] = {0};	//100k
	char *str = conv_buff;
	int len;

	str[0] = 0;
	str += sprintf(str, "%d,%d,%d,", p->path[0], p->path[p->cur - 1], p->cur - 1);

	n->stime = p->cur - 1;
	n->cdf = update_convolution(p, -1, &len);

	double_to_string(str, n->cdf, len);
	fwrite(conv_buff, sizeof(char), strlen(conv_buff), f);
}

double cal_peerProb(double *cdf, int len, int time)
{
	double *ptr = cdf;
	int i, tmp = TSLOT;
	for(i=0; i<len; i++) {
		tmp += i * TSLOT;
		if(tmp >= time)
			return cdf[i];
	}
	return 1.;
}

void write_record(MATRIX *G)
{
	FILE *fc = fopen("./ccdf.csv", "w");
	FILE *fpa = fopen("./path.csv", "w");
	FILE *fpb = fopen("./peer_pro.csv", "w");
	int i, j;
	
	char *buff;
	size_t buff_len = 10240;
	buff = (char *)calloc(buff_len, sizeof(char));
	char pbuff[10240] = {0};

	int id, cur;
	for(i=0; i<NODE_NUM; i++) {
		id = i;
		cur = 0;
		for(j=0; j<NODE_NUM; j++) {
			//m_path(G, i, j, NODE_NUM) = path(type, G, i, j, NODE_NUM);
			//do the convolution
			PATH *tmp =  m_path(G, i, j, NODE_NUM);
			if(tmp == NULL) {
				if(i != j)
					fprintf(fpb, "%d,%d,%lf\n", i, j, 0.);
				continue;
			}

			//write path
			char *str;
			str = pbuff;
			str[0] = 0;
			int_to_string(str, tmp->path, tmp->cur);
			fwrite(pbuff, sizeof(char), strlen(pbuff), fpa);

			cur++;
			p_ccdf[i] = (PEER *)realloc(p_ccdf[i], cur * sizeof(PEER));
			PEER *pt = p_ccdf[i] + cur - 1;

			//init it
			pt->id = j;
			pt->cdf = NULL;
			pt->stime = 0;

			//write cdf
#ifndef FIXED_ROUTE
			{
					NEIGHBOR key, *res;
					key.id = j;
					res = bsearch(&key, node[i].nei, node[i].cur, sizeof(NEIGHBOR), nei_cmp);
					if(!res) {
						pt->stime = -1;	//mark 0 as end...
						fprintf(fpb, "%d,%d,%lf\n", i, j, 0.);
						continue;
					}

					if(expect_false(buff_len < res->num * 10))
						array_needsize(char, buff, buff_len, res->num * 10, array_zero_init);
					str = buff;
					str[0] = 0;
					str += sprintf(str, "%d,%d,", i, res->id);

					int k;
					double *tmp_res = (double *)calloc(res->num, sizeof(double));
					tmp_res[0] = res->delay_pdf[0];
					for(k=1; k<res->num; k++)
						tmp_res[k] = res->delay_pdf[k] + tmp_res[k - 1];
				
					pt->stime = 1;
					pt->cdf = tmp_res;

					double_to_string(str, tmp_res, res->num);
					fwrite(buff, sizeof(char), strlen(buff), fc);
					fprintf(fpb, "%d,%d,%lf\n", i, j, cal_peerProb(pt->cdf, res->num, TIME_WINDOW * 60));
			}
#else
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
#endif
		}
		p_ccdf[id] = (PEER *)realloc(p_ccdf[id], (cur + 1) * sizeof(PEER));
		memset(p_ccdf[id] + cur, 0, sizeof(PEER));
	}

	free(buff);
	fclose(fc);
	fclose(fpa);
	fclose(fpb);
}

static double r1(void)
{
	int r;
	while((r = rand()) == 0);
	return (double)r / (double)RAND_MAX;
}

#ifdef INTR_TEST
static double r2(int l, gsl_rng *r)
{
	double res;
	double high = 0.1 * (double)(l + 1);
	double low = 0.1 * (double)l;
	
	for(;;) {
		res = gsl_rng_uniform(r);
		if(res >= low && res <= high)
			return res;
	}
}
#endif

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

double cal_probability(double *cdf, int len, int stime, int time)
{
	int rtime = stime * TSLOT, i;
	double *res = cdf;
	for(i=0; i<len; i++) {
		if(*res == 1)
			return *res;

		if(rtime >= time)
			return *res;

		rtime += TSLOT;
		res++;
	}

	return *(--res) > 1 ? 1 : *res;
}

PATH *build_direct_path(int s, int d, int num, bool symmetry)
{
	if(s == d)
		return NULL;

//check if node s could reach node d before building the path
	NEIGHBOR _d, *_s;
	_d.id = d;
	_s = bsearch(&_d, node[s].nei, node[s].cur, sizeof(NEIGHBOR), nei_cmp);
	if(!_s)
		return NULL;

	if(symmetry) {
		_d.id = s;
		_s = bsearch(&_d, node[d].nei, node[d].cur, sizeof(NEIGHBOR), nei_cmp);
		if(!_s)
			return NULL;
	}
	
	PATH *res = (PATH *)calloc(1, sizeof(PATH));
	res->path = (int *)calloc(num, sizeof(int));
	res->num = num;
	res->path[res->cur++] = s;
	res->path[res->cur++] = d;

	return res;
}


double get_probability(const MATRIX *G, int s, int i, int j, int time)
{
	PATH *si = NULL, *sj = NULL, *sji = NULL, *ij = NULL, *ji = NULL, *f = NULL;
	double *new_cdf;
	double res = 0;
	int cdf_len;

#ifndef FIXED_ROUTE
	if(s != -1 && (s == i || i == j)) {
		f = build_direct_path(s, j, 4, true);
		if(f == NULL)
			return 0;

		f->path[f->cur++] = s;
		f->path[f->cur++] = j;
		goto PRO_CAL;
	}
#endif

	if(s == -1)
		si = NULL;
	else {
#ifdef FIXED_ROUTE
		si = m_path(G, s, i, NODE_NUM);	
#else	
		si = build_direct_path(s, i, 2, false);
#endif
	}

#ifdef FIXED_ROUTE
	ij = m_path(G, i, j, NODE_NUM);
	ji = m_path(G, j, i, NODE_NUM);
#else
	ij = build_direct_path(i, j, 2, false);
	ji = build_direct_path(j, i, 2, false);
#endif

	if((s != -1 && si == NULL) || 
		(ij == NULL) || 
		(ji == NULL)) {
//		printf("no path %d-%d-%d\n", s, i, j);
		goto GP_EXIT;
	}
	
	if(si == NULL)
		sj = ij;
	else
		sj = path_merge(si, ij);
		
	sji = path_merge(sj, ji);
	f = path_merge(sji, ij);	
	
PRO_CAL:
	new_cdf = update_convolution(f, time, &cdf_len);
	//res = new_cdf[cdf_len - 1];
	res = cal_probability(new_cdf, cdf_len, f->cur - 1, time);
	free(new_cdf);

	if(si && sj) {
		free(sj->path);
		free(sj);
		sj = NULL;
	}
	if(sji) {
		free(sji->path);
		free(sji);
		sji = NULL;
	}	
	if(f) {
		free(f->path);
		free(f);
		f = NULL;
	}

GP_EXIT:
#ifndef FIXED_ROUTE
	if(si) {
		free(si->path);
		free(si);
		si = NULL;
	}
	if(ij) {
		free(ij->path);
		free(ij);
		ij = NULL;
	}
	if(ji) {
		free(ji->path);
		free(ji);
		ji = NULL;
	}
#endif

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

double cal_mrev(const MATRIX *G, PINFO *n, int s, const char *x, int time)
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
		if(res->id == id && res->stime != -1)
			return res;
		res++;
	}
	return NULL;
}

PINFO *build_node_info(const peerlist *p, int s, int time)
{
	PINFO *res = (PINFO *)calloc(NODE_NUM, sizeof(PINFO));

	int i;
	for(i=0; i<NODE_NUM; i++) {
		
		PEER *tmp = peer_search(p[s], i);
		if(tmp == NULL)
			continue;
		
//		res[i].probability = cal_probability(tmp->cdf, tmp->stime, time);
		res[i].interest = g_interest[i];
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

#ifdef USE_NLOPT
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
#endif

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

void knapsack(dp_item **t, int num, int k, int *w, double *v, const MATRIX *G, PINFO *n, int s, int time)
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

void item_init(int **weight, double **value, int num, const MATRIX *G, PINFO *n, int s, int time)
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

void generate_adv(M_NODE *node, int id, int time)
{
	M_NODE *n = &(node[id]);
	int left = n->buff_len - n->buff_cur - NODE_NUM;
	if(left < 0)
		array_needsize(FDATA, n->buffer, n->buff_len, n->buff_len - left, array_zero_init);

	int i;
	for(i=0; i<NODE_NUM; i++) {
		if(node[i].candidate)
			continue;
		
		FDATA *b = &(n->buffer[n->buff_cur++]);
		b->src = n->id;
		b->dest = i;
		b->type = FILE_ADV;
		b->stime = time;
	}
}

int find_next_hop(const MATRIX *G, M_NODE *n, int dest)
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

#if 0
double rev_probability(const MATRIX *G, M_NODE *node, int dst)
{
	int i;
	double res = 1;
	for(i=0; i<NODE_NUM; i++) {
		if(node[i].candidate == false || node[i].source)
			continue;

		res = res * (1 - get_probability(G, sn, i, dst, TIME_WINDOW * 60));
	}
	
	return res;
}
#endif

//new trans method
void can_trans(M_NODE *n, M_NODE *node, int stime, int rtime, const MATRIX *G)
{
	int i;
	for(i=0; i<NODE_NUM; i++) {
		if(n->neighbor[i] && node[i].recv_file == false) {
			if(node[i].candidate) {
				node[i].recv_file = true;
				node[i].have_file = true;
				generate_adv(node, i, rtime);

				if(sim_type == 0) {
					c_can_log.storage_load[i]++;
					c_can_log.comm_load[i]++;
					c_can_log.comm_load[n->id]++;
				}
				else {
					d_can_log.storage_load[i]++;
					d_can_log.comm_load[i]++;
					d_can_log.comm_load[n->id]++;
				}
			}
			else{
				if(node[i].interest) {
					node[i].recv_file = true;
					
					if(sim_type == 0) {
						c_can_log.storage_load[i]++;
						c_can_log.comm_load[i]++;
						c_can_log.comm_load[n->id]++;
					}
					else {
						d_can_log.storage_load[i]++;
						d_can_log.comm_load[i]++;
						d_can_log.comm_load[n->id]++;
					}
				}
			}

			if(node[i].interest) {
				sim_delay += rtime - stime;
				sim_delivery++;
				
				sim_rev += node[i].interest * PRICE;
				//fprintf(f_log, "#%d (%lf) \t", i, node[i].interest);

				if(sim_type == 0) {
					c_dst_log.delay[i] += rtime - stime;
					c_dst_log.receivings[i]++;
				}
				else {
					d_dst_log.delay[i] += rtime - stime;
					d_dst_log.receivings[i]++;
				}
			}
		}
	}
}

void relay_trans(M_NODE *n, M_NODE *node, int stime, int rtime, const MATRIX *G)
{
	int i;
	FDATA *b;
	
	for(i=0; i<n->buff_cur; i++) {
		b = &(n->buffer[i]);
		if((b->type == FILE_TRANS || b->type == FILE_DISTRIBUTION) && 
			n->neighbor[b->dest] &&
			node[b->dest].recv_file == false) {
			n->have_file = false;
			node[b->dest].recv_file = true;
			if(node[b->dest].candidate) {
				node[b->dest].have_file = true;
				generate_adv(node, b->dest, rtime);
			}
				
			if(node[b->dest].interest) {
				sim_delivery++;
				sim_delay += rtime - stime;
				sim_rev += node[b->dest].interest * PRICE;

				//fprintf(f_log, "#%d (%lf) \t", i, node[i].interest);

				if(sim_type == 0) {
					c_dst_log.delay[b->dest] += rtime - stime;
					c_dst_log.receivings[b->dest]++;
				}
				else {
					d_dst_log.delay[b->dest] += rtime - stime;
					d_dst_log.receivings[b->dest]++;
				}
			}

			if(sim_type == 0) {
				c_can_log.storage_load[b->dest]++;
				//do this later on buffer handle...
				//can_log.storage_load[i]--;
				c_can_log.comm_load[b->dest]++;
				c_can_log.comm_load[i]++;
			}
			else {
				d_can_log.storage_load[b->dest]++;
				d_can_log.comm_load[b->dest]++;
				d_can_log.comm_load[i]++;
			}
				
			remove_data(b, n, i);
		}
	}
}

void file_trans(M_NODE *n, M_NODE *node, int stime, int rtime, const MATRIX *G)
{
	if(n->candidate)
		can_trans(n, node, stime, rtime, G);
	else
		relay_trans(n, node, stime, rtime, G);
}

void Adv2Req(M_NODE *node, FDATA *b, int time, int id)
{
	M_NODE *cn = &(node[b->dest]);
	cn->candidate_list[b->src] = 1;
	int c, k, flag = 0;
	FDATA *dbuff;
	
	//!!!!! could always generate REQ when there are neighbors around...
	node[b->dest].candidate_list[b->src] = 1;
#ifdef IMPROVE_SCHEME
//!!!!! could always generate REQ when there are neighbors around...
	for(c=0; c<NODE_NUM; c++) {
		flag = 0;
		if(cn->candidate_list[c] == 0)
			continue;
#else
		c = b->src;
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
			dbuff->stime = time;

			if(sim_type == 0) {
				c_can_log.comm_load[cn->id]++;
				c_can_log.comm_load[id]++;
			}
			else {
				d_can_log.comm_load[cn->id]++;
				d_can_log.comm_load[id]++;
			}
		}
#ifdef IMPROVE_SCHEME
	}
#endif
}

#ifdef IMPROVE_SCHEME
void Req2Trans(M_NODE *node, FDATA *b, int time, int id, char *neighbor)
{
	//!!!!! could be recv by different candidates...
	int c;
	for(c=0; c<NODE_NUM; c++) {
		if(neighbor[c] &&
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
				dbuff->stime = time;

				if(sim_type == 0) {
					c_can_log.comm_load[id]++;
					c_can_log.comm_load[c]++;
				}
				else {
					d_can_log.comm_load[id]++;
					d_can_log.comm_load[c]++;
				}
			}
	}
}
#else
void Req2Trans(M_NODE *node, FDATA *b, int time, int id)
{
	FDATA tmp_data = {.src = b->dest, .dest = b->src, .type = FILE_TRANS};
	if(check_data_duplicate(&node[b->dest], tmp_data) == 0) {
		FDATA *dbuff;
		if(expect_false(node[b->dest].buff_cur == node[b->dest].buff_len))
			array_needsize(FDATA, node[b->dest].buffer, node[b->dest].buff_len, node[b->dest].buff_len + 1, array_zero_init);
		dbuff = &(node[b->dest].buffer[node[b->dest].buff_cur++]);
		dbuff->src = b->dest;
		dbuff->dest = b->src;
		//must be sent back by multiple hops
		dbuff->type = FILE_TRANS;
		dbuff->stime = time;

		if(sim_type == 0) {
			c_can_log.comm_load[id]++;
			c_can_log.comm_load[b->dest]++;
		}
		else {
			d_can_log.comm_load[id]++;
			d_can_log.comm_load[b->dest]++;
		}
	}
}
#endif

void handle_buffer(M_NODE *n, M_NODE *node, int stime, int rtime, const MATRIX *G)
{
	int i, next_hop;
	FDATA *b;
	
	for(i=0; i<n->buff_cur; i++) {
		b = &(n->buffer[i]);
		if(b->type == FILE_DISTRIBUTION || b->type == FILE_TRANS) {
			if(node[b->dest].recv_file == false) {
				//find possible next hop
				next_hop = find_next_hop(G, n, b->dest);
				if(next_hop != -1) {
					send_data(&node[next_hop], b);
					node[next_hop].have_file = true;
					if(sim_type == 0)
						c_can_log.storage_load[next_hop]++;
					else
						d_can_log.storage_load[next_hop]++;
					
					remove_data(b, n, i);
					if(n->candidate == false) {
						if(sim_type == 0)
							c_can_log.storage_load[n->id]--;
						else
							d_can_log.storage_load[n->id]--;
						n->have_file = false;
					}

					if(sim_type == 0) {
						c_can_log.comm_load[next_hop]++;
						c_can_log.comm_load[n->id]++;
					}
					else {
						d_can_log.comm_load[next_hop]++;
						d_can_log.comm_load[n->id]++;
					}
				}
			}
			else {
				remove_data(b, n, i);
				if(n->candidate == false) {
					if(sim_type == 0)
						c_can_log.storage_load[n->id]--;
					else
						d_can_log.storage_load[n->id]--;
					n->have_file = false;
				}
			}
		}
		else if(b->type == FILE_ADV) {
			if(node[b->dest].recv_file == false) {
				next_hop = find_next_hop(G, n, b->dest);
				if(next_hop != -1) {
					if(next_hop == b->dest)
						Adv2Req(node, b, rtime, n->id);
					else {
						send_data(&node[next_hop], b);
						
						if(sim_type == 0) {
							c_can_log.comm_load[next_hop]++;
							c_can_log.comm_load[n->id]++;
						}
						else {
							d_can_log.comm_load[next_hop]++;
							d_can_log.comm_load[n->id]++;
						}
					}
					remove_data(b, n, i);
				}
			}
			else {
				remove_data(b, n, i);
			}
		}
		else if(b->type == FILE_REQ) {
			if(n->recv_file == true && b->src == n->id) {
				remove_data(b, n, i);
			}
			else {
				next_hop = find_next_hop(G, n, b->dest);
				if(next_hop != -1) {
					if(next_hop == b->dest)
						Req2Trans(node, b, rtime, n->id);
					else {
						send_data(&node[next_hop], b);

						if(sim_type == 0) {
							c_can_log.comm_load[next_hop]++;
							c_can_log.comm_load[n->id]++;
						}
						else {
							d_can_log.comm_load[next_hop]++;
							d_can_log.comm_load[n->id]++;
						}
					}
					remove_data(b, n, i);
				}
			}
		}
	}
}

void handle_node(M_NODE *list[], int num, M_NODE *node, int stime, int rtime, const MATRIX *G)
{
	int i;
	M_NODE *n;
	
	for(i=0; i<num; i++) {
		n = list[i];
		if(n->have_file)
			file_trans(n, node, stime, rtime, G);
			
		handle_buffer(n, node, stime, rtime, G);
	}
}

decl_shuffle(nptr);

void node_communication(char *neighbor, M_NODE *node, int stime, int rtime, const MATRIX *G)
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
	if(cur == 0)
		return;
	
	shuffle_nptr(list, cur);	//shuffle the neighbors before starting data trans
	handle_node(list, cur, node, stime, rtime, G);
}

//not a generic function... just for simulation_loop to pre-generate the adv data for source node/candidate nodes
void generate_initData(M_NODE *n, char *candidate, int stime, int source_node, int type)
{ 
	int j;
	for(j=0; j<NODE_NUM; j++) {
		if(j == n->id || j == source_node)
			continue;

		if(expect_false(n->buff_cur == n->buff_len))
			array_needsize(FDATA, n->buffer, n->buff_len, n->buff_len + 1, array_zero_init);

		FDATA *t;
		if(!candidate[j]) {
		//generate the adv data
			t = &(n->buffer[n->buff_cur++]);
			t->src = n->id;
			t->dest = j;
			t->stime = stime;
			t->type = FILE_ADV;
		}
		else if(!type){
		//generate file distribution data
			t = &(n->buffer[n->buff_cur++]);
			t->src = n->id;
			t->dest = j;
			t->stime = stime;
			t->type = FILE_DISTRIBUTION;
		}
	}
}

void write_can_log(void)
{
#define FILE_SIZE	1000
	int i;
	fprintf(f_can, "###centralized comm_load###\n");
	for(i=0; i<NODE_NUM; i++)
		fprintf(f_can, "%.5lf\n", c_can_log.comm_load[i]/(double)c_runtime);
	fprintf(f_can, "\n###centralized storage_load###\n");
	for(i=0; i<NODE_NUM; i++)
		fprintf(f_can, "%.5lf\n", c_can_log.storage_load[i]/(double)c_runtime*FILE_SIZE);

	fprintf(f_can, "\n###distributed comm_load###\n");
	for(i=0; i<NODE_NUM; i++)
		fprintf(f_can, "%.5lf\n", d_can_log.comm_load[i]/(double)d_runtime);
	fprintf(f_can, "\n###distributed storage_load###\n");
	for(i=0; i<NODE_NUM; i++)
		fprintf(f_can, "%.5lf\n", d_can_log.storage_load[i]/(double)d_runtime*FILE_SIZE);
#undef FILE_SIZE
}

void write_dst_log(void)
{
	int i;
	fprintf(f_dst, "###centralized delay###\n");
	for(i=0; i<NODE_NUM; i++)
		fprintf(f_dst, "%.5lf\n", c_dst_log.delay[i]/(double)c_runtime);
	fprintf(f_dst, "\n###centralized receivings###\n");
	for(i=0; i<NODE_NUM; i++)
		fprintf(f_dst, "%.5lf\n", c_dst_log.receivings[i]/(double)c_runtime);

	fprintf(f_dst, "\n###distributed delay###\n");
	for(i=0; i<NODE_NUM; i++)
		fprintf(f_dst, "%.5lf\n", d_dst_log.delay[i]/(double)d_runtime);
	fprintf(f_dst, "\n###distributed receivings###\n");
	for(i=0; i<NODE_NUM; i++)
		fprintf(f_dst, "%.5lf\n", d_dst_log.receivings[i]/(double)d_runtime);
}


void simulation_loop(int source_node, int stime, long wtime, char *candidate, PINFO *info, const MATRIX *G, int type/*0 - centrilized; 1 - distributed*/)
{
	if(type * (type - 1) != 0) {
		printf("unknown simulation type\n");
		return;
	}

	sim_delivery = 0;
	sim_delay = 0;
	sim_rev = 0;
	
/*
	start loop...
*/
	//init mobilie node at first...
	M_NODE *node = (M_NODE *)calloc(NODE_NUM, sizeof(M_NODE));
	M_NODE *tn;
	int i, j;
	for(i=0; i<NODE_NUM; i++) {
		node[i].id = i;
		node[i].neighbor = (char *)calloc(NODE_NUM, sizeof(char));
		node[i].candidate_list = (char *)calloc(NODE_NUM, sizeof(char));

		if(i == source_node) {	//init source node
			node[i].candidate = true;
			node[i].have_file = true;
			node[i].recv_file = true;
			if(sim_type == 0)
				c_can_log.storage_load[i] = 1;
			else
				d_can_log.storage_load[i] = 1;
			node[i].source = true;

			array_needsize(FDATA, node[i].buffer, node[i].buff_len, NODE_NUM, array_zero_init);
			//generate data
			tn = &(node[i]);
			generate_initData(tn, candidate, stime, source_node, type);
		}
		else {
			if(candidate[i]) {
				node[i].candidate = true;
				if(type) {
					tn = &(node[i]);
					tn->have_file = true;
					tn->recv_file = true;
					if(sim_type == 0)
						c_can_log.storage_load[i] = 1;
					else
						d_can_log.storage_load[i] = 1;
					generate_initData(tn, candidate, stime, source_node, type);
				}
			}
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
		if(time < stime)	//not start yet
			continue;

//		n1--; n2--;
		if(time == rtime) {
			neighbor[n1] = 1;
			neighbor[n2] = 1;
			continue;
		}

		if(time > rtime) {
			//handle previous time slot record
			node_communication(neighbor, node, stime, rtime, G);
			_dprintf("communicate @ %d\n\n", time);

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
		if(time < stime)
			continue;

//		i--; j--;
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
int get_max_obRev(int node, int stime, int wtime, int k, int *best_candidate,  int *ob_time, int *ob_candidate, PINFO *n, const MATRIX *G)
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
	long wtime_s = wtime;

	while((read = getline(&line, &len, f)) != -1) {
		int i, j, time;
		sscanf(line, "%d,%d,%d", &i, &j, &time);
		if(time < stime)	//not start yet...
			continue;

//		i--; j--;
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
			double rev = cal_mrev(G, n, -1, x, wtime * OB_WINDOW/*actual file transfer time window*/);	//we do not need to consider the path from souce node to candidate
			if(cnt <= k) {
			//within the observing time window
				ob_test[neighbor] = 1;
				if(ob_max < rev) {
					ob_max = rev;
					*ob_candidate = neighbor;
				}
			}
			else {	//start the real candidate selection
				if(*ob_time == -1)
					*ob_time = time - stime;
				
				can_test[neighbor] = 1;
				if(rev >= ob_max) {
					if(s_can == -1) {	//first node we meet with higher value than the ob_max after the first k meeting event
						s_can = neighbor;
						_dprintf("rev: %lf\n", rev);
					}
					else {
						*best_candidate = neighbor;
						res = max(res, rev);
					}
				}
			}
		}
	}

	_dprintf("max rev: %lf\n", res);

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

//		n1--; n2--;
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

//int *meeting_node = get_meetingNodes(source_node, stime, ob_events, &num);
int *get_meetingNodes(int source_node, int stime, int events, int *num, int *ob_time)
{
	FILE *f = fopen("./mobility.csv", "r");
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	int cnt = 2, cur = 0;
	int *res = (int *)calloc(cnt, sizeof(int));
	res[0] = -1;
	res[1] = -1;

	int i, j, time;
	while((read = getline(&line, &len, f)) != -1) {
		sscanf(line, "%d,%d,%d", &i, &j, &time);
		if(time < stime)
			continue;

		if(time - stime > events)
			break;

		if(i == source_node || j == source_node) {
			int key = (source_node - i == 0) ? j : i;
			int *find = (int *)bsearch(&key, res, cur, sizeof(int), int_cmp);
			//already meet this node
			if(find)
				continue;
			//new meeting node
			if(cur == cnt)
				res = (int *)realloc(res, ++cnt * sizeof(int));
			res[cur++] = key;
			qsort(res, cur, sizeof(int), int_cmp);
		}
	}
	*num = cur;
	*ob_time = time - stime;
	
	fclose(f);
	free(line);
	return res;
}

void merge_set(int *a1, int *a2, int n2)
{
	int i, j = 0;
	for(i=0; i<n2; i++) {
		if(a2[i] == -1)
			a2[i] = a1[j++];
	}
}

void map_set(int *s, int num, char *x)
{
	memset(x, 0, NODE_NUM * sizeof(char));
	int i;
	for(i=0; i<num; i++)
		x[s[i]] = 1;
}

//select_mcandidate(source_node, stime, wtime/OB_WINDOW, final.value, &ob_can, max_weight - 1, n, G);
int *select_mcandidate(int source_node, int stime, int events/*ob time*/, int wtime, double ob_max, int *ob_can, int num, int *cnt/*how many nodes we meet in the real candidates selection time window*/, int **meeting_node, PINFO *n, const MATRIX *G)
{
	int *s_can = (int *)calloc(num, sizeof(int));	//the selected candidates
	int cur = 0;
	*meeting_node = (int *)calloc(2, sizeof(int));
	int m_len = 2;
	*cnt = 0;
	FILE *f = fopen("./mobility.csv", "r");
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	char tested[NODE_NUM];
	char x[NODE_NUM];
	int tmp_s[num];
	memset(x, 0, NODE_NUM * sizeof(char));
	memset(tested, 0, sizeof(char) * NODE_NUM);
	int wtime_s = wtime;
	int m_events = 0;

	while((read = getline(&line, &len, f)) != -1) {
		int i, j, time;
		sscanf(line, "%d,%d,%d", &i, &j, &time);
		if(time < stime + events)
			continue;
		if(time > stime + wtime_s) 
			break;

		if(i == source_node || j == source_node) {
			//observing time has been passed
			int neighbor = (source_node - i == 0) ? j : i;
			if(tested[neighbor])
				continue;

			if(*cnt == m_len)
				*meeting_node = (int *)realloc(*meeting_node, ++m_len * sizeof(int));
			(*meeting_node)[(*cnt)++] = neighbor;
			tested[neighbor] = 1;

			if(cur == num /*we have already got enough candidates*/) 
				continue;
			
			int m, flag = 0;
			for(m=0; m<num; m++) {
				if(ob_can[m] == -1)
					continue;
				if(ob_can[m] == neighbor) {
					flag = 1;
					ob_can[m] = -1;
					s_can[cur++] = neighbor;
					break;
				}
			}
			if(flag)
				continue;

			for(m=0; m<num; m++) {
				if(ob_can[m] == -1)
					continue;
				
				memcpy(tmp_s, ob_can, sizeof(int) * num);
				tmp_s[m] = neighbor;
				merge_set(s_can, tmp_s, num);
				map_set(tmp_s, num, x);
				double tmp_res = cal_mrev(G, n, -1, x, wtime * OB_WINDOW - events);
				if(tmp_res >= ob_max) {
					s_can[cur++] = neighbor;
					ob_can[m] = -1;
					break;
				}
			}
		}
	}

	fclose(f);
	free(line);
	if(cur == num)
		return s_can;
	else {
		free(s_can);
		return NULL;
	}
}

int distributed_simulation(int source_node, int stime, int wtime, PINFO *n, const MATRIX * G, int *fail)
{
	int best_candidate = -1;
	char x[NODE_NUM];
	memset(x, 0, sizeof(char) * NODE_NUM);

	sim_delivery = 0;
	sim_delay = 0;
	sim_rev = 0;
	int total_events = get_meetingEvent(source_node, stime, wtime/OB_WINDOW/*the candidate selection time window*/);
	int ob_events = wtime / OB_WINDOW * DRATIO;	//time for observe candidates
	int ob_time = -1;
	int ob_candidate = -1;
	int type = 1;
#ifdef SINGLE_SELECT
	int candidate = get_max_obRev(source_node, stime, wtime/OB_WINDOW, total_events * DRATIO, &best_candidate, &ob_time, &ob_candidate, n, G);
	if(candidate >= 0) {
		x[candidate] = 1;
		if(candidate != best_candidate)
			*fail = 1;
	}
#if 0
	else {	//work around way if we do not get the candidate...
		candidate = ob_candidate;
		x[candidate] = 1;
		type = 0;
		ob_time = 0;
	}
#endif
	else {
		ob_time = -1;
		goto CLEANUP;
	}
	
	_dprintf("selected candidate: %d, best candidate: %d\n", candidate, best_candidate);
#else
	int num, i, j;
	int *meeting_node = get_meetingNodes(source_node, stime, ob_events, &num, &ob_time);
	int *meeting_node2 = NULL, *candidate = NULL;
	if(num == 0) {
		free(meeting_node);
		return ob_time;
	}
	
	int max_weight = CAN_NUM;	//the total num of nodes we could choose is (max_weight - 1)
	int ob_can[max_weight - 1];
	
	int *i_weight = (int *)calloc(num, sizeof(int));
	double *i_value = (double *)calloc(num, sizeof(double));
	item_init(&i_weight, &i_value, num, G, n, -1, wtime - ob_events);

	dp_item *dp = (dp_item *)calloc(max_weight * num, sizeof(dp_item));
	for(i=0; i<max_weight * num; i++) {
		dp[i].selection = (char *)calloc(NODE_NUM, sizeof(char));
		dp[i].id = meeting_node[i/max_weight];
	}

	knapsack(&dp, num, max_weight, i_weight, i_value, G, n, -1, wtime - ob_events);

	dp_item final = dp[i-1];
	j = 0;
	_dprintf("candidates selected from observing time: ");
	
	for(i=0; i<NODE_NUM; i++)
		if(final.selection[i]) {
			_dprintf("#%d\t", i);
			ob_can[j++] = i;
		}
	_dprintf("stime: %d\n", stime);
	_dprintf("max rev from observing time: %lf\n", final.value);
	if(j == 0) {
		ob_time = -1;
		goto CLEANUP;
	}

	int num2;
	candidate = select_mcandidate(source_node, stime, ob_events, wtime/OB_WINDOW, final.value, ob_can, j/*num of ob candidates*/, &num2, &meeting_node2, n, G);
	if(candidate == NULL) {
#ifdef _DEBUG
		_dprintf("no candidate found\n");
		for(i=0; i<num2; i++)
			_dprintf("+%d\t", meeting_node2[i]);
		_dprintf("\n");
#endif
		ob_time = -1;
		goto CLEANUP;
	}

	//we can compute what is the best choice based on meeting_node2 and num2...
	int *i_weight2 = (int *)calloc(num2, sizeof(int));
	double *i_value2 = (double *)calloc(num2, sizeof(double));
	item_init(&i_weight2, &i_value2, num2, G, n, -1, wtime - ob_events);

	dp_item *dp2 = (dp_item *)calloc(max_weight * num2, sizeof(dp_item));
	for(i=0; i<max_weight * num2; i++) {
		dp2[i].selection = (char *)calloc(NODE_NUM, sizeof(char));
		dp2[i].id = meeting_node2[i/max_weight];
	}

	knapsack(&dp2, num2, max_weight, i_weight2, i_value2, G, n, -1, wtime - ob_events);

	dp_item final2 = dp2[i-1];
	
#ifdef _DEBUG
	_dprintf("best candidates could be selected: ");
	for(i=0; i<NODE_NUM; i++)
		if(final2.selection[i]) 
			_dprintf("#%d\t", i);
	_dprintf("\n");
	_dprintf("max rev could obtain: %lf\n", final2.value);
#endif

	map_set(candidate, j, x);

#ifdef _DEBUG
	_dprintf("actual candidates: ");
	for(i=0; i<NODE_NUM; i++)
		if(x[i])
			_dprintf("#%d\t", i);
	_dprintf("\n");
#endif

	if(memcmp(x, final2.selection, NODE_NUM * sizeof(char)))
		*fail = 1;

	free(i_weight2);
	free(i_value2);
	for(i=0; i<max_weight * num2; i++) {
		if(dp2[i].selection)
			free(dp2[i].selection);
	}
	free(dp2);

#if 0
#ifndef SINGLE_SELECT
	if(candidate == NULL) {	//if we do not get the candidate set for multiple-selection distributed algorithm
		memcpy(x, final.selection, NODE_NUM);
		type = 0;
		ob_time = 0;
	}
#endif
#endif

	//start real file distributioin. Start time: stime+wtime; file validate time: wtime
	stime += ob_events;
	wtime -= ob_events;
#endif

	simulation_loop(source_node, stime, wtime, x, n, G, type);

CLEANUP:
#ifndef SINGLE_SELECT
	if(meeting_node)
		free(meeting_node);
	if(meeting_node2)
		free(meeting_node2);
	if(candidate)
		free(candidate);
	if(i_weight)
		free(i_weight);
	if(i_value)
		free(i_value);
	if(dp) {
		for(i=0; i<max_weight * num; i++) {
			if(dp[i].selection)
				free(dp[i].selection);
		}
		free(dp);
	}
#endif

	_dprintf("ob time: %d\n", ob_time);
	return ob_time;
//	simulation_start(source_node, stime, ob_time, wtime, n, G, ob_rev);
}

int get_start_time(int source_node, int stime)
{
	FILE *f = fopen("./mobility.csv", "r");
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	int res = 0;

	while((read = getline(&line, &len, f)) != -1) {
		int i, j, time;
		sscanf(line, "%d,%d,%d", &i, &j, &time);
//		i--; j--;
		if(i == source_node || j == source_node) {
			if(time >= stime) {
				res = time;
				break;
			}
		}
	}

	free(line);
	fclose(f);
	return (res - 30);
}

///////////////////////////////////////////////////////////////////////////////////////
MATRIX *build_graph(void)
{
	int i, j;
	bool flag;
	MATRIX *res = distribution_to_matrix(&flag);
	
	if(flag)
		folyd_warshall(res, NODE_NUM);
	else {
		double dist[NODE_NUM];
		for(i=0; i<NODE_NUM; i++)
			dijkstra(res, dist, i, NODE_NUM);
	}

	for(i=0; i<NODE_NUM; i++) {
		for(j=0; j<NODE_NUM; j++) {
			m_path(res, i, j, NODE_NUM) = path(flag, res, i, j, NODE_NUM);
		}
	};
	
	return res;
}

void write_peers_delayD(MATRIX *G)
{
	int i, j;
	P_DELAY *peer_D = info_collect(G);
	for(i=0; i<NODE_NUM; i++) {	//node id
		NODE *n = &node[i];
		for(j=0; j<NODE_NUM; j++) {	//neighbor id
			P_DELAY *tmp = matrix(peer_D, i, j, NODE_NUM);
			if(tmp->len == 0)
				continue;
			
			qsort(tmp->delay, tmp->cur, sizeof(int), unsigned_cmp);
			
			NEIGHBOR key, *res;
			key.id = j;
			res = bsearch(&key, n->nei, n->cur, sizeof(NEIGHBOR), nei_cmp);
			bool flag = false;
			if(!res) {
				if(expect_false(n->cur == n->num))
					array_needsize(NEIGHBOR, n->nei, n->num, n->num + 1, array_zero_init);

				res = &(n->nei[n->cur++]);
				res->id = j;
				flag = true;
			}
			
			{
				res->delay_average = cal_pdf(tmp->delay, tmp->cur);
				res->num = pdf_cur;
				if(res->delay_pdf)
					free(res->delay_pdf);
				
				res->delay_pdf = (double *)calloc(pdf_cur, sizeof(double));
				memcpy(res->delay_pdf, pdf, pdf_cur * sizeof(double));
			}
			if(flag)
				qsort(n->nei, n->cur, sizeof(NEIGHBOR), nei_cmp);

			free(tmp->delay);
		}
	}

	free(peer_D);
}

double double_getn_str(char *str, int n)
{
	char *token = NULL;
	int i = 0;
	double res;
	token = strtok(str, ",");

	while(token != NULL) {
		res = atof(token);
		if(i == n)
			break;

		token = strtok(NULL, ",");
		i++;
	}

	return res;
}

double read_intrfile(int id)
{
	FILE *f = fopen("../node_intrfile.csv", "r");
	while(lockf(fileno(f), F_TEST, 0L) == -1)	//test if the file is ready for use
		sleep(1);

	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	int i = 0;
	double res;

	while((read = getline(&line, &len, f)) != -1) {
		if(i == id) {
			res = double_getn_str(line, sim_no);
			break;
		}
		i++;
	}

	free(line);
	fclose(f);
	return res;
}

double *interest_gen(void)
{
	double *res = (double *)calloc(NODE_NUM, sizeof(double));
	FILE *f = fopen("node_interest.log", "w");
	
	srand(_seed);
	int i;
	for(i=0; i<NODE_NUM; i++) {
#ifndef INTR_TEST
		res[i] = r1();
#else
		res[i] = read_intrfile(i);
#endif
		fprintf(f, "%lf\n", res[i]);
	}

	fclose(f);
	return res;
}

bool init_var(void)
{
	FILE *f = fopen("sim.conf", "r");
	if(f == NULL)
		return false;
	
	char *line = NULL;
	size_t len = 0;
	ssize_t read;

	while((read = getline(&line, &len, f)) != -1) {
		char name[32] = {0};
		int value;
		sscanf(line, "%s = %d", name, &value);
		if(strcmp(name, "TIME_WINDOW") == 0)
			TIME_WINDOW = value;
		else if(strcmp(name, "CAN_NUM") == 0)
			CAN_NUM = value;
		else if(strcmp(name, "PRICE") == 0)
			PRICE = value;
		else if(strcmp(name, "COST") == 0)
			COST = value;
		else if(strcmp(name, "OB_WINDOW") == 0)
			OB_WINDOW = value;
		else if(strcmp(name, "DRATIO") == 0)
			DRATIO = (double)value / 10.;
		else
			printf("unknow parameter\n");
	}

	free(line);
	fclose(f);
	return true;
}

decl_shuffle(double);

#ifdef INTR_TEST
void make_interfile(void)
{

	FILE *f;
	if(access("../node_intrfile.csv", F_OK) != -1) 
		return;
	else {
		f = fopen("../node_intrfile.csv", "w");
		if(lockf(fileno(f), F_LOCK, 0L) == -1)
			perror("lockf LOCK");
	}
	
	const gsl_rng_type *T;
	gsl_rng *r;
	int i;
	double *rate = (double *)calloc(SIM_ROUND, sizeof(double));

	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	for(i=0; i<NODE_NUM; i++) {
		int j;
		for(j=0; j<10; j++) 
			rate[j] = r2(j, r);
		
		shuffle_double(rate, SIM_ROUND);
		
		for(j=0; j<SIM_ROUND-1; j++)
			fprintf(f, "%.5lf,", rate[j]);
		fprintf(f, "%.5lf\n", rate[j]);
	}

	gsl_rng_free(r);
	free(rate);
	
	rewind(f);
	if(lockf(fileno(f), F_ULOCK, 0L) == -1)
		perror("lockf ULOCK");
	fclose(f);
}
#endif

void sim_setup(const char *filename, MATRIX **G)
{
	init_var();

//do not change the functions seq ;-)
	cal_distribution(filename);

#ifdef INTR_TEST
	make_interfile();
#endif
	
	p_ccdf = (peerlist *)calloc(NODE_NUM, sizeof(peerlist));
	g_interest = interest_gen();
	
	*G = build_graph();

#ifndef FIXED_ROUTE
	write_peers_delayD(*G);
#endif

	write_distribution("./pdf.csv");

	write_record(*G);
}

void sim_unit_run(int src_node, const MATRIX *G)
{
	int i, wtime = TIME_WINDOW * 60;	//time window is xx min
	int source_node = src_node;
	char x[NODE_NUM];
	memset(x, 0, NODE_NUM * sizeof(char));
	
	int max_weight = CAN_NUM;	//the total num of nodes we could choose is (max_weight - 1)
	int *i_weight = (int *)calloc(NODE_NUM, sizeof(int));
	double *i_value = (double *)calloc(NODE_NUM, sizeof(double));
	dp_item *dp = (dp_item *)calloc(max_weight * NODE_NUM, sizeof(dp_item));
	dp_item final;
	int t_time = 0, tcnt = 0, mcnt = 0, cnt = 0, stime;
	double average_rev = 0, average_delivery = 0, average_delay = 0;

	PINFO *ni = build_node_info(p_ccdf, source_node, wtime);
	
#ifdef USE_SOLVER
	write_node_interest(ni);
	write_cdf(G, source_node, wtime);
#endif

#ifdef RANDOM_TEST
	//cal_mrev() func test...
	x[2] = 1;
	x[8] = 1;
	x[12] = 1;
	x[25] = 1;
	double rev = cal_mrev(G, ni, source_node, x, wtime);
	_dprintf("\n=============OPT RESULTS===================\n\n");
	_dprintf("random rev: %lf\n", rev);
	_dprintf("random candidates:\t");
	for(i=0; i<NODE_NUM; i++) {
		if(x[i] != 0)
			_dprintf("#%d\t", i);
	}
	_dprintf("\n\n");
#endif

	item_init(&i_weight, &i_value, NODE_NUM, G, ni, source_node, wtime);
	for(i=0; i<max_weight * NODE_NUM; i++) {
		dp[i].selection = (char *)calloc(NODE_NUM, sizeof(char));
		dp[i].id = i/max_weight;
	}

	knapsack(&dp, NODE_NUM, max_weight, i_weight, i_value, G, ni, source_node, wtime);
	final = dp[i-1];
	
	fprintf(f_src, "src %d\n", source_node);
	fprintf(f_src, "m_rev %lf\n", final.value);
	_dprintf("max rev: %lf\n", final.value);
	_dprintf("candidates:\t");
	fprintf(f_scan, "candidates for src %d:\t", source_node);
	for(i=0; i<NODE_NUM; i++) {
		if(final.selection[i] != 0)
			fprintf(f_scan, "#%d\t", i);
	}
	fprintf(f_scan, "\n");

	if(final.value == 0)
		goto END_RUN;
	
/////////////////////////SIMULATION///////////////////////////////////////////
#ifdef CENTRA_SIM	
	_dprintf("\n============Centralized SIMULATION RESULTS===================\n\n");
	sim_type = 0;		//centralized sim
	for(;;) {
		stime = get_start_time(source_node, t_time);
		if(stime < 0)
			break;
		
		simulation_loop(source_node, stime, wtime, final.selection, ni, G, 0);
		t_time = stime + 6 * TSLOT;	//roundup to the next time slot
		cnt++;

		if(sim_rev == 0) {
			_dprintf("!!!!!sim rev == 0... @ %d!!!!!\n", stime);
			mcnt++;
			continue;
		}
		
		average_rev += sim_rev;
		average_delivery += sim_delivery;
		average_delay += (double)sim_delay/(double)sim_delivery;
		tcnt++;

#ifdef _DEBUG
		if(sim_rev > final.value)
			fprintf(f_log, "\n$$$$sharings: %d @ %d (%d)$$$$\n\n", sim_delivery, stime, tcnt);
#endif

		_dprintf("@@@@@stime: %d@@@@@\n\n", stime);
	}
	c_runtime += tcnt;
	fprintf(f_src, "total_cnt %d\n", cnt);
	fprintf(f_src, "c_cnt %d\n", tcnt);
	fprintf(f_src, "c_rev %lf\n", average_rev/(double)tcnt - CAN_NUM*COST);
	fprintf(f_src, "c_sharings %lf\n", average_delivery/(double)tcnt);
	fprintf(f_src, "c_delay %lf\n", average_delay/(double)tcnt);
	
	_dprintf("running times: %d (missed: %d)\n", tcnt, mcnt);
	_dprintf("sim revenue: %lf\n", average_rev/(double)tcnt);
	_dprintf("total sharing: %lf\n", average_delivery/(double)tcnt);
	_dprintf("average delay: %lf\n", average_delay/(double)tcnt);
#endif

#ifdef DISTRI_SIM	
	_dprintf("\n==============Distributed SIMULATION RESULTS=====================\n\n");
	t_time = 0; tcnt = 0; mcnt = 0;
	average_rev = 0; average_delivery = 0; average_delay = 0;
	sim_type = 1;		//distributed sim
	for(;;) {
		int fail = 0;
		stime = get_start_time(source_node, t_time);
		if(stime < 0)
			break;
			
		int ob_delay = distributed_simulation(source_node, stime, wtime, ni, G, &fail);
		t_time = stime + 6 * TSLOT;

	
		if(ob_delay == -1 || sim_rev== 0) {
			_dprintf("sim rev == 0... @ %d\n", stime);
			mcnt++;
			continue;
		}

		average_rev += sim_rev;
		average_delivery += sim_delivery;
		average_delay += (double)sim_delay/(double)sim_delivery + (double)ob_delay;
		tcnt++;
		mcnt += fail;

		//if(sim_rev > final.value)
		//	printf("sharings: %d @ %d\n", sim_delivery, stime);
	}
	d_runtime += tcnt;
	fprintf(f_src, "d_cnt %d\n", cnt - mcnt);
	fprintf(f_src, "d_rev %lf\n", average_rev/(double)tcnt - CAN_NUM*COST);
	fprintf(f_src, "d_sharings %lf\n", average_delivery/(double)tcnt);
	fprintf(f_src, "d_delay %lf\n", average_delay/(double)tcnt);
	
	_dprintf("running times: %d (missed: %d)\n", tcnt, mcnt);
	_dprintf("sim revenue: %lf\n", average_rev/(double)tcnt);
	_dprintf("total sharing: %lf\n", average_delivery/(double)tcnt);
	_dprintf("average delay: %lf\n", average_delay/(double)tcnt);
#endif
	_dprintf("\n===============DONE===========================\n\n");
	
////////////////////////////////////////////CLEANUP//////////////////////////////
END_RUN:
	fprintf(f_src, "\n");
	free(i_weight);
	free(i_value);
	for(i=0; i<max_weight * NODE_NUM; i++) {
		if(dp[i].selection)
			free(dp[i].selection);
	}
	free(dp);
	free(ni);
}

void sim_clean(MATRIX *G)
{
	node_free();
	convolution_free();
	
	int i;
	for(i=0; i<NODE_NUM * NODE_NUM; i++) {
		PATH *tmp = G[i].path;
		if(tmp && tmp->path) {
			free(tmp->path);
			free(tmp);
		}
	}
	free(G);
	
	free_peerlist(p_ccdf, NODE_NUM);
	free(g_interest);
}

void sim_log_init(void)
{
	f_src = fopen("src.log", "w");
	f_can = fopen("can.log", "w");
	f_dst = fopen("dst.log", "w");
	f_scan = fopen("s_can.log", "w");

	f_log = fopen("debug.log", "w");

	c_can_log.comm_load= (double *)calloc(NODE_NUM, sizeof(double));
	c_can_log.storage_load= (double *)calloc(NODE_NUM, sizeof(double));
	
	c_dst_log.delay= (double *)calloc(NODE_NUM, sizeof(double));
	c_dst_log.receivings= (double *)calloc(NODE_NUM, sizeof(double));

	d_can_log.comm_load= (double *)calloc(NODE_NUM, sizeof(double));
	d_can_log.storage_load= (double *)calloc(NODE_NUM, sizeof(double));
	
	d_dst_log.delay= (double *)calloc(NODE_NUM, sizeof(double));
	d_dst_log.receivings= (double *)calloc(NODE_NUM, sizeof(double));
}

void sim_log_end(void)
{
	fclose(f_src);
	fclose(f_can);
	fclose(f_dst);
	fclose(f_log);
	fclose(f_scan);

	free(c_can_log.comm_load);
	free(c_can_log.storage_load);
	
	free(d_dst_log.delay);
	free(d_dst_log.receivings);

	free(d_can_log.comm_load);
	free(d_can_log.storage_load);
	
	free(d_dst_log.delay);
	free(d_dst_log.receivings);
}

void save_res(const char *name)
{
	char cmd[64] = {0};
	sprintf(cmd, "./mkres.sh %s", name);
	cmd_system(cmd);
}

int main(int argc, char *argv[])
{
	if(argc < 4) {
		printf("lack of parameters\n");
		return 1;
	}
/*
 * arg1 graph.csv
 * arg2 sim No.
 * arg3 sim res name
 * */
	
	int sn;
	MATRIX *G;

	sim_no = atoi(argv[2]);
	sim_setup(argv[1], &G);
	sim_log_init();

	printf("sim%s start...\n", argv[2]);
	for(sn=0; sn<NODE_NUM; sn++) {
		_dprintf("source node: %d\n", sn);
#ifdef _DEBUG
		fprintf(f_log, "\n\n$$$$NODE: %d$$$$$$\n", sn);
#endif
		sim_unit_run(sn, G);
	}

	write_can_log();
	write_dst_log();

	sim_clean(G);
	sim_log_end();

	save_res(argv[3]);
	
	printf("sim%s done...\n", argv[2]);
	return 0;
}

