#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "common.h"

typedef struct peer_t {
	int id;
	int stime;
	double *cdf;
} PEER;

typedef PEER * peerlist;

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

//test...
int main(void)
{
	peerlist *p = get_cdf("ccdf.csv", 76);
	PEER *res = peer_search(p[0], 68);

	printf("%lf\n", get_probability(res, 500));

	free_peerlist(p, 76);
	
	return 0;
}

