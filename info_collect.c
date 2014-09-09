#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>

#include "common.h"
#include "shortest_path.h"
#include "simulation.h"
#include "distribution.h"

#include "info_collect.h"

//extern NODE *node;
extern unsigned int NODE_NUM;

void generate_bmsg(M_NODE *node, int time)
{
	node->buff_len = NODE_NUM;
	node->buff_cur = NODE_NUM;
	node->buffer = (FDATA *)calloc(NODE_NUM, sizeof(FDATA));

	int i;
	for(i=0; i<NODE_NUM; i++) {
		if(i == node->id)
			continue;

		node->buffer[i].src = node->id;
		node->buffer[i].dest = i;
		node->buffer[i].stime = time;
		node->buffer[i].type = BMSG;
	}
}

void fill_bmsg(M_NODE *node, int time)
{
	int i;
	for(i=0; i<NODE_NUM; i++) {
		if(i == node->id)
			continue; 
		
		if(node->buffer[i].stime == -1)
			node->buffer[i].stime = time;
	}
}


void broadcast_data(M_NODE *node, FDATA *data, int time, P_DELAY *r)
{
	if(data->dest == node->id) {
		//reach the dest, record the delay
		P_DELAY *p = matrix(r, data->src, data->dest, NODE_NUM);

		if(expect_false(p->cur == p->len))
			array_needsize(unsigned int, p->delay, p->len, p->len + 1, array_zero_init);

		p->delay[p->cur++] = time - data->stime;
	}
	else {
		if(expect_false(node->buff_cur == node->buff_len))
			array_needsize(FDATA, node->buffer, node->buff_len, node->buff_len + 1, array_zero_init);

		memcpy(&(node->buffer[node->buff_cur++]), data, sizeof(FDATA));
	}
}

void start_boradcast(M_NODE *list[], int num, M_NODE *node, int time, MATRIX *G, P_DELAY *r)
{
	int i, j, next_hop;
	M_NODE *n;
	FDATA *b;
	
	for(i=0; i<num; i++) {
		n = list[i];	
		if(n->buff_len == 0)
			generate_bmsg(n, time);
		else
			fill_bmsg(n, time);

		for(j=0; j<n->buff_cur; j++) {
			b = &(n->buffer[j]);
			if(b->type == BMSG) {	
				//find possible next hop
				next_hop = find_next_hop(G, n, b->dest);
				if(next_hop != -1) {
					broadcast_data(&node[next_hop], b, time, r);

					//remove the data
					if(b->src == n->id) 
						b->stime = -1;
					else 
						remove_data(b, n, j);
				}
			}
		}
	}
}

void broadcast_com(char *neighbor, M_NODE *node, int time, MATRIX *G, P_DELAY *r)
{
	M_NODE *list[NODE_NUM];
	int i, cur = 0;
	for(i=0; i<NODE_NUM; i++) {
		if(neighbor[i]) {
			memcpy(node[i].neighbor, neighbor, sizeof(char) * NODE_NUM);
			node[i].neighbor[i] = 0;	//clear self bit
			list[cur++] = &(node[i]);
		}
	}

	start_boradcast(list, cur, node, time, G, r);
}

P_DELAY* info_collect(MATRIX *G)
{
	P_DELAY *res = (P_DELAY *)calloc(NODE_NUM * NODE_NUM, sizeof(P_DELAY));
	
	M_NODE *node = (M_NODE *)calloc(NODE_NUM, sizeof(M_NODE));
	int i;
	for(i=0; i<NODE_NUM; i++) {
		node[i].id = i;
		node[i].neighbor = (char *)calloc(NODE_NUM, sizeof(char));
	}
	
	char neighbor[NODE_NUM];
	memset(neighbor, 0, sizeof(char) * NODE_NUM);
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	FILE *f = fopen("./mobility.csv", "r");
	int stime = atoi(cmd_system("head -1 ./mobility.csv | cut -d , -f 3"));
	int rtime = stime;

	while((read = getline(&line, &len, f)) != -1) {
		int n1, n2, time;
		sscanf(line, "%d,%d,%d", &n1, &n2, &time);
		if(time < stime)	//not start yet
			continue;

//		n1--; n2--;
		if(time == rtime) {
			neighbor[n1] = 1;
			neighbor[n2] = 1;
			continue;
		}

		if(time > rtime) {
			//handle previous time slot record
			broadcast_com(neighbor, node, rtime, G, res);
			
			//reset record
			rtime = time;
			memset(neighbor, 0, sizeof(char) * NODE_NUM);
			neighbor[n1] = 1;
			neighbor[n2] = 1;
		}	
	}

	//free M_NODE
	for(i=0; i<NODE_NUM; i++) {
		if(node[i].neighbor)
			free(node[i].neighbor);
		if(node[i].buffer)
			free(node[i].buffer);
	}
	free(node);
	free(line);
	fclose(f);

	return res;
}

