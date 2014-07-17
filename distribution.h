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

typedef struct func_data_t {
	int snode, wtime;
	MATRIX *graph;
	PINFO *node;
} FUNC_DATA;

typedef struct DP_t {
	double value;
	char *selection;	//NODE_NUM
	int id;
} dp_item;

typedef PEER * peerlist;

bool cal_distribution(const char *inputF);
int find_next_hop(MATRIX *G, M_NODE *n, int dest);

#endif

