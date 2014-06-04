#ifndef _DISTRIBUTION_H
#define _DISTRIBUTION_H

typedef struct neighbor_t {
	double *delay_pdf;
	double delay_average;
	unsigned int id;
	int num;
} NEIGHBOR;

typedef struct node_t {
	NEIGHBOR *nei;
	int num;
	int cur;
} NODE;

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

bool cal_distribution(const char *inputF, const char *outputF);

#endif

