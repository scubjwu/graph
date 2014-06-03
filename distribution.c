#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>

#include "common.h"
#include "shortest_path.h"
#include "convolution.h"

#include "distribution.h"

static const double INF = DBL_MAX/2 - 1;

static char *line = NULL;
static NODE *node;
static unsigned int NODE_NUM = 0;
static FILE *fp;
static unsigned int *delay_t = NULL;
static int dnum_t = 0;
static double *pdf = NULL;
static int pdf_len = 0;
static int pdf_cur = 0;

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
#define TSLOT	60
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
#undef TSLOT
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

//double *convolution(double *a1, int n1, double *a2, int n2)
void do_convolution(FILE *f, char *buff, PATH *p)
{
	static char conv_buff[10240] = {0};
	char *str;

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
	//write back a1 with len n1
	str = conv_buff;
	str[0] = 0;
	str += sprintf(str, "%d,%d,%d,", p->path[0], p->path[p->cur - 1], p->cur - 1);
	//make cdf
	for(i=1; i<n1; i++)
		a1[i] += a1[i-1];

	double_to_string(str, a1, n1);
	fwrite(conv_buff, sizeof(char), strlen(conv_buff), f);
	//free memory
	for(i=0; i<c; i++)
		free(D[i]);
}

int main(int argc, char *argv[])
{
	cal_distribution(argv[1], argv[2]);

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

	FILE *fc = fopen("./ccdf.csv", "w");
	FILE *fpa = fopen("./path.csv", "w");
	int i, j;
	
	char *buff;
	size_t buff_len = 1024;
	buff = (char *)calloc(buff_len, sizeof(char));
	char pbuff[1024] = {0};

	for(i=0; i<NODE_NUM; i++)
		for(j=0; j<NODE_NUM; j++) {
			m_path(G, i, j, NODE_NUM) = path(flag, G, i, j, NODE_NUM);
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
				double tmp_res[res->num];
				tmp_res[0] = res->delay_pdf[0];
				for(k=1; k<res->num; k++)
					tmp_res[k] = res->delay_pdf[k] + tmp_res[k - 1];

				double_to_string(str, tmp_res, res->num);
				fwrite(buff, sizeof(char), strlen(buff), fc);
				continue;
			}
			//need to recal the delay_cdf
			do_convolution(fc, buff, tmp);
		}

	free(buff);
	fclose(fc);
	fclose(fpa);
	node_free();
	return 0;
}

