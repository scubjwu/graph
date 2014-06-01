#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>

#include "common.h"

#define MALLOC_ROUND 1024

static inline size_t array_nextsize(size_t elem, size_t cur, size_t cnt)
{
	size_t ncur = cur + 1;

	do
		ncur <<= 1;
	while(cnt > ncur);

	if(elem * ncur > MALLOC_ROUND - sizeof(void *) * 4) {
		ncur *= elem;
		ncur = (ncur + elem + (MALLOC_ROUND - 1) + sizeof(void *) * 4) & ~(MALLOC_ROUND - 1);
		ncur = ncur - sizeof(void *) * 4;
		ncur /= elem;
	}

	return ncur;
}

//elem is the size of individual element; *cur is the current array size; cnt is the new size
void * declare_noinline array_realloc(size_t elem, void *base, size_t *cur, size_t cnt)
{
	*cur = array_nextsize(elem, *cur, cnt);
	base = realloc(base, elem * *cur);
	return base;
}

void array_zero_init(void *p, size_t op, size_t np, size_t elem)
{
	memset(p + (op * elem), 0, (np - op) * elem);
}

