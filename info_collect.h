#ifndef _INFO_COL_H
#define _INFO_COL_H

typedef struct peer_delay {
	unsigned int *delay;
	int len;
	int cur;
} P_DELAY;

P_DELAY* info_collect(MATRIX *G);

#endif
