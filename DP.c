#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>

#include "common.h"
#include "convolution.h"

static const unsigned int _seed = INT_MAX - 1;

typedef struct peer_t {
	int id;
	int stime;
	double *cdf;
} PEER;

typedef struct pinfo_t {
	double probability;
	double interest;
} PINFO;

typedef PEER * peerlist;

static inline double r1(void)
{
	return (double)rand() / (double)RAND_MAX;
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

peerlist *get_cdf(const char *filename, int num)
{
	FILE *f = fopen(filename, "r");
	if(f == NULL) {
		perror("fopen");
		return NULL;
	}

	peerlist *p = (peerlist *)calloc(num, sizeof(peerlist));
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	int id = 0, cur = 0;
	while((read = getline(&line, &len, f)) != -1) {
		char *token;
		int i, j;
		token = strtok(line, ",");	
		i = atoi(token);
		if(i > id) {
			p[id] = (PEER *)realloc(p[id], (cur + 1) * sizeof(PEER));
			//mark as the end of the list
			memset(p[id] + cur, 0, sizeof(PEER));

			id = i;
			cur = 0;
		}
		token = strtok(NULL, ",");

		cur++;
		p[i] = (PEER *)realloc(p[i], cur * sizeof(PEER));
		PEER *tmp = (p[i] + cur - 1);

		int p = 0, cn = 4, cp = 0;
		double *cdf = (double *)calloc(cn, sizeof(double));

		while(token != NULL) {
			if(p == 0)
				tmp->id = atoi(token);
			else if(p == 1)
				tmp->stime = atoi(token);
			else {
				if(cp == cn)
					cdf = (double *)realloc(cdf, ++cn * sizeof(double));
				cdf[cp++] = atof(token);
			}
			p++;
			token = strtok(NULL, ",");
		}
		tmp->cdf = cdf;
	}

	return p;
}

double get_probability(PEER *p, int time)
{
#define TUNIT	60
	int rtime = p->stime * TUNIT;
	double *res = p->cdf;
	for(;;) {
		if(*res == 1)
			return *res;

		if(rtime >= time)
			return *res;

		rtime += TUNIT;
		res++;
	} 
#undef TUNIT
}

int get_cdf_len(double *a)
{
	double *tmp = a;
	int res = 1;
	for(;;) {
		if(*tmp == 1)
			break;
		tmp++;
		res++;
	}

	return res;
}

double cal_mrev(PINFO *n, int s, char *X, int num, int time)
{
#define PRICE	50
	int i, j;
	double m, sum = 0;
	for(j=0; j<num; j++) {
		if(n[j].interest == 0)
			continue;
		//do the mulplication
		m = 1;
		for(i=0; i<num; i++) {
			if(X[i] == 0)
				continue;
			
#if 0
			double *si, *ij;
			PEER *p_si, *p_ij;
			p_si = peer_search(p[s], i);
			p_ij = peer_search(p[i], j);
			if(p_si == NULL || p_ij == NULL)
				continue;

			PEER p_sj = {0};
			p_sj.stime = p_si->stime + p_ij->stime;
			int si_len, ij_len;
			si_len = get_cdf_len(p_si->cdf);
			ij_len = get_cdf_len(p_ij->cdf);
			p_sj.cdf = convolution(p_si->cdf, si_len, p_ij->cdf, ij_len);
			double _m = 1 - get_probability(&p_sj, time);
#endif
			P_sij = get_probability(s, i, j, time);
			m = m * (1 - P_sij);
		}
		//sum together
		sum += (1 - m) * n[j].interest * PRICE;
	}

	return sum;
#undef PRICE
}

PINFO *build_node_info(peerlist *p, int s, int time, int num)
{
	PINFO *res = (PINFO *)calloc(num, sizeof(PINFO));

	int i;
	for(i=0; i<num; i++) {
		PEER *tmp = peer_search(p[s], i);
		if(tmp == NULL)
			continue;

		res[i].probability = get_probability(tmp, time);
		res[i].interest = r1();
	}

	return res;
}

//test...
int main(int argc, char *argv[])
{
	if(argc < 4) {
		printf("need para\n");
		exit(1);
	}

	srand(_seed);
	int num = atoi(argv[2]);
	int source_node = atoi(argv[1]);
	int live_time = atoi(argv[3]);

	peerlist *p = get_cdf("ccdf.csv", num);
	PINFO *node = build_node_info(p, source_node, live_time, num);

	char x[num];
	memset(x, 0, num * sizeof(char));
	x[6] = 1;
	cal_mrev(p, node, source_node, x, num, live_time);

	free_peerlist(p, num);
	free(node);
	return 0;
}

