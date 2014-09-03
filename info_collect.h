#ifndef _INFO_COL_H
#define _INFO_COL_H

typedef struct peer_delay {
	unsigned int *delay;
	size_t len;
	size_t cur;
} P_DELAY;

P_DELAY* info_collect(MATRIX *G);

#endif
