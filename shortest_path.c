#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include <float.h>

#include "common.h"
#include "shortest_path.h"

const double INF = DBL_MAX/2 - 1;

//parent(i, j) == -1 means no path between node i & j
//graph(i, j) == INF means no path betwen node i & j
void folyd_warshall(MATRIX *graph, int num)
{
	int i, j, k;
	for(i=0; i<num; i++)
		for(j=0; j<num; j++) {
			if(m_weight(graph, i, j, num) == 0 ||
				m_weight(graph, i, j, num) == INF)
				m_parent(graph, i, j, num) = -1;
			else
				m_parent(graph, i, j, num) = i;
		}

	for(k=0; k<num; k++)
		for(i=0; i<num; i++)
			for(j=0; j<num; j++) {
				double dist = m_weight(graph, i, k, num) + m_weight(graph, k, j, num);
				if(dist < m_weight(graph, i, j, num)) {
					m_weight(graph, i, j, num) = dist;
					m_parent(graph, i, j, num) = m_parent(graph, k, j, num);
				}
			}
}

static bool path_dijkstra(PATH *path, MATRIX *graph, int s, int src, int dst, int num)
{
	int p = m_parent(graph, s, dst, num);
	if(p == -1) {
		printf("no path\n");
		return false;
	}

	if(p == src) {
		if(expect_false(path->cur == path->num))
			array_needsize(int ,path->path, path->num, path->num + 1, array_zero_init);
		path->path[path->cur++] = p;
		return true;
	}

	path_dijkstra(path, graph, s, src, p, num);
	path_dijkstra(path, graph, s, p, dst, num);
}

static bool path_fw(PATH *path, MATRIX *graph, int src, int dst, int num)
{
	if(m_parent(graph, src, dst, num) == -1) {
		printf("no path\n");
		return false;
	}

	int i = m_parent(graph, src, dst, num);
	if(i == src) {
		if(expect_false(path->cur == path->num))
			array_needsize(int, path->path, path->num, path->num + 1, array_zero_init);
		path->path[path->cur++] = i;
		return true;
	}

	path_fw(path, graph, src, i, num);
	path_fw(path, graph, i, dst, num);
}

PATH *path(int type, MATRIX *graph, int src, int dst, int num)
{
	PATH *p = (PATH *)calloc(1, sizeof(PATH));
	array_needsize(int, p->path, p->num, 2, array_zero_init);

	switch(type) {
	case 0:
		if(path_fw(p, graph, src, dst, num)) {
			if(expect_false(p->cur == p->num))
				array_needsize(int, p->path, p->num, p->num + 1, array_zero_init);
			p->path[p->cur++] = dst;
			return p;
		}
		break;
	case 1:
		if(path_dijkstra(p, graph, src, src, dst, num)) {
			if(expect_false(p->cur == p->num))
				array_needsize(int, p->path, p->num, p->num + 1, array_zero_init);
			p->path[p->cur++] = dst;
			return p;
		}
		break;
	default:
		printf("unknown type\n");
		break;
	}

	free(p);
	return NULL;
}

void dijkstra(MATRIX *graph, double *dist, int s, int num)
{
	int i;
	bool done[num];

	for(i=0; i<num; i++) {
		dist[i] = INF;
		done[i] = false;
	}
	dist[s] = 0;

	for(;;) {
		int u = -1, v;
		double bestDist = DBL_MAX;
		for(i=0; i<num; i++) 
			if(!done[i] && dist[i] < bestDist) {
				u = i;
				bestDist = dist[i];
			}
		
		if(bestDist == DBL_MAX)
			break;

		for(v=0; v<num; v++) {
			if(!done[v] && m_weight(graph, u, v, num) != INF) {
				if(dist[v] > dist[u] + m_weight(graph, u, v, num)) {
					dist[v] = dist[u] + m_weight(graph, u, v, num);
					m_parent(graph, s, v, num) = u;
				}
			}
		}
		done[u] = true;
	}
}

//test
int main(void)
{
	double graph[5][5] = {{0,10,INF,5,INF},{INF,0,1,2,INF},{INF,INF,0,INF,4},{INF,3,9,0,2},{7,INF,6,INF,0}};

	int n = 5, i, j;
	MATRIX *G = (MATRIX *)calloc(n*n, sizeof(MATRIX));
	for(i=0; i<n; i++)
		for(j=0; j<n; j++) {
			m_weight(G, i, j, n) = graph[i][j];
			m_parent(G, i, j, n) = -1;
		}

	double dist[n];
	for(i=0; i<n; i++)
		dijkstra(G, dist, i, n);
	m_path(G, 0, 2, n) = path(1, G, 0, 2, n);

//	folyd_warshall(G, n);
//	m_path(G, 4, 0, n) = path(0, G, 4, 0, n);

	PATH *p = m_path(G, 0, 2, n);
	if(p == NULL) 
		return 0;

	for(i=0; i<p->cur; i++)
		printf("%d-", p->path[i]);
	printf("NIL\n");

	return 0;
}

