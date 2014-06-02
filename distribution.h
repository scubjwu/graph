#ifndef _DISTRIBUTION_H
#define _DISTRIBUTION_H

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

bool cal_distribution(const char *inputF, const char *outputF);

#endif

