#ifndef _GRAPH_FW_H
#define _GRAPH_FW_H

typedef struct path_t {
	int *path;
	size_t num;
	size_t cur;
} PATH;

typedef struct matrix_t {
	double weight;
	int parent;
	PATH *path;
} MATRIX;

void folyd_warshall(MATRIX *graph, int num);
PATH *path(int type, MATRIX *graph, int src, int dst, int num);
void dijkstra(MATRIX *graph, double *dist, int s, int num);

#endif

