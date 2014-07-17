#ifndef SIMULATIOIN_H
#define SIMULATION_H

typedef enum data_type_t {
	FILE_DISTRIBUTION = 0,
	FILE_ADV,
	FILE_REQ,
	FILE_TRANS,
	BMSG
} DATA_TYPE;

typedef struct node_data_t {
	int src;
	int dest;
	int stime;
	DATA_TYPE type;
} FDATA;

typedef struct mobile_node_t {
	int id;

	char *neighbor;	//store the neighbor at each time slot

	char *candidate_list;	//store the candidate info for requestor

	FDATA *buffer;	//data carried by the node
	int buff_cur;
	int buff_len;

	bool source;
	bool candidate;	//if the node is candidate
	bool have_file;	//if the node has received the file

	double interest;	//the probability the node is intereting in shared file
} M_NODE;

#define remove_data(data, node)	\
{	\
	FDATA *__last = &(node->buffer[node->buff_cur - 1]);	\
	FDATA *__tmp = data;	\
	data = __last;	\
	__last = __tmp;	\
	node->buff_cur--;	\
	j--;	\
}

#endif
