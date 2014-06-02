#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>

#include "common.h"
#include "shortest_path.h"

#include "distribution.h"

static const double INF = DBL_MAX/2 - 1;

static char *line = NULL;
static NODE *node;
static unsigned int NODE_NUM = 0;
static FILE *fp;
static unsigned int *delay_t = NULL;
static int dnum_t = 0;
static double *cdf = NULL;
static int cdf_len = 0;
static int cdf_cur = 0;

static int int_cmp(const void *n1, const void *n2)
{
	return (*(unsigned int *)n1 - *(unsigned int *)n2);
}

static double cal_cdf(unsigned int *array, int num)
{
#define TSLOT	60
	if(expect_false(cdf_len < num))
		array_needsize(double, cdf, cdf_len, num, array_zero_init);
	
	int i, time = TSLOT, sum = 0;
	double total_delay = 0;
	cdf_cur = 0;
	for(i=0; i<num; i++) {
		if(array[i] > time) {
			cdf[cdf_cur++] = (double)sum/(double)num;
			time += TSLOT;
		}
		sum++;
		total_delay += array[i];
	}
	cdf[cdf_cur++] = (double)sum/(double)num;	//should grows to 1 at last

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
	p->delay_average = cal_cdf(delay, num);
	p->num = cdf_cur;
	p->delay_cdf = (double *)calloc(cdf_cur, sizeof(double));
	memcpy(p->delay_cdf, cdf, cdf_cur * sizeof(double));
}

static void get_node_info(unsigned int id)
{
#define NEIGHBOR_THRESHOLD	50
	size_t len = 0;
	ssize_t read;
	unsigned int p_neighbor = 0;
	int cur_t = 0;

	NODE *n = &(node[id-1]);
	array_needsize(NEIGHBOR, n->nei, n->num, 2, array_zero_init);

	while((read = getline(&line, &len, fp)) != -1) {
		unsigned int time, node, neighbor;
		sscanf(line, "%d,%d,%d", &time, &node, &neighbor);

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

	if(cdf) {
		free(cdf);
		cdf = NULL;
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
					if(node[i].nei[j].delay_cdf) {
						free(node[i].nei[j].delay_cdf);
						node[i].nei[j].delay_cdf = NULL;
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
				double_to_string(tmp, nei->delay_cdf, nei->num);
				fwrite(buff, sizeof(char), strlen(buff), f);
			}
		}
	}

	free(buff);
	buff = NULL;
	fclose(f);

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


	node_free();
	return 0;
}

