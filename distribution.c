#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include "common.h"

typedef struct neighbor_t {
	double *delay_cdf;
	double delay_average;
	unsigned int id;
	int num;
} NEIGHBOR;

typedef struct node_t {
	NEIGHBOR *nei;
	int num;
	int cur;
} NODE;

static char *line = NULL;
static NODE *node;
static unsigned int NODE_NUM = 0;
static FILE *fp;
static unsigned int *delay_t = NULL;
static int dnum_t = 0;
static double *cdf = NULL;
static int cdf_len = 0;
static int cdf_cur = 0;

static char *cmd_system(const char *cmd)
{
#define BUFLEN	512
	char *res = "";
	char buf[BUFLEN] = {0};
	FILE *f;
	
	f = popen(cmd, "r");
	while(fgets(buf, BUFLEN-1,f) != NULL)
		res = buf;

	if(f != NULL)
		pclose(f);

	return res;
#undef BUFLEN
}

static int int_cmp(const void *n1, const void *n2)
{
	return (*(unsigned int *)n1 - *(unsigned int *)n2);
}

double cal_cdf(unsigned int *array, int num)
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

void neighbor_wb(unsigned int *delay/*neighbor delay distribution*/, int num/*num of inter contact time record*/, int nei_id/*neighbor id*/, NODE *n/*node*/)
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

#define NEIGHBOR_THRESHOLD	50
void get_node_info(unsigned int id)
{
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
}

void exit_clean(void)
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

	fclose(fp);
}

void double_to_string(char *str, const double *array, int len)
{
	int i;
	for(i=0; i<len-1; i++)
		str += sprintf(str, "%.5lf,", array[i]);

	sprintf(str, "%.5lf\r\n", array[i]);
}

void write_distribution(void)
{
	char *buff;
	size_t buff_len = 1024;
	buff = (char *)calloc(buff_len, sizeof(char));

	FILE *f = fopen("./res.D", "w");
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
}

int main(int argc, char *argv[])
{
	char cmd[512] = {0};

	fp = fopen(argv[1], "r");
	if(fp == NULL) {
		perror("fopen");
		exit(1);
	}

	sprintf(cmd, "cut -d , -f 2 %s | sort | uniq | wc -l", argv[1]);
	NODE_NUM = atoi(cmd_system(cmd));	//num of nodes in the network
	node = (NODE *)calloc(NODE_NUM, sizeof(NODE));

	int i;
	for(i=1; i<=NODE_NUM; i++)
		get_node_info(i);

	write_distribution();

	exit_clean();

#if 0
	for(i=0; i<num; i++) {
		if(node[i].cur)
			printf("node %d has %d neighbors\n", i, node[i].cur);
	}
#endif
}

