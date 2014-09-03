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

char *cmd_system(const char *cmd)
{
#define BUFLEN	512
	char *res = "";
	static char buf[BUFLEN];
	FILE *f;
	
	f = popen(cmd, "r");
	memset(buf, 0, BUFLEN * sizeof(char));
	while(fgets(buf, BUFLEN-1, f) != NULL)
		res = buf;

	if(f != NULL)
		pclose(f);

	return res;
#undef BUFLEN
}

void double_to_string(char *str, const double *array, int len)
{
	int i;
	for(i=0; i<len-1; i++)
		str += sprintf(str, "%.5lf,", array[i]);

	sprintf(str, "%.5lf\r\n", array[i]);
}

void int_to_string(char *str, const int *array, int len)
{
	int i;
	for(i=0; i<len-1; i++)
		str += sprintf(str, "%d,", array[i]);

	sprintf(str, "%d\r\n", array[i]);
}

int irand(int n)
{
	int r, rand_max = RAND_MAX - (RAND_MAX % n);
	/* reroll until r falls in a range that can be evenly
	 * distributed in n bins.  Unless n is comparable to
	 * to RAND_MAX, it's not *that* important really. */
	while ((r = rand()) >= rand_max);
	return r / (rand_max / n);
}

