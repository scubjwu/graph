#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "common.h"

//complexity O(n1*n2)
static size_t _con_len = 0;
static double *_con_buff;

double __attribute__((optimize("O0"))) *convolution(double *a1, int n1, double *a2, int n2)
{
	int i, j, n = n1 - 1 + n2, p = 0;

	if(_con_len < n * n1) {
		_con_len = n * n1;
		_con_buff = (double *)realloc(_con_buff, _con_len * sizeof(double));
	}
	memset(_con_buff, 0, _con_len * sizeof(double));

	double *tmp = _con_buff;
//	double tmp[n][n1];	//this may overflow on stack...
//	memset(tmp, 0, n * n1 * sizeof(double));

	double *res = (double *)calloc(n, sizeof(double));

	for(i=0; i<n1; i++) {
		for(j=0; j<n2; j++) {
		//	tmp[j + i][i] = a1[i] * a2[j];
			*matrix(tmp, j + i, i, n1) = a1[i] * a2[j];
		}
	}

	for(i=0; i<n; i++) {
		double sum = 0;
		for(j=0; j<n1; j++) {
		//	sum += tmp[i][j];
			sum += *matrix(tmp, i, j, n1);
		}
		res[i] = sum;
	}

	//the len of res equals to (n1 + n2 - 1)
	return res;
}

void convolution_free(void)
{
	if(_con_buff) {
		free(_con_buff);
		_con_buff = NULL;
	}
}

#if 0
//test...
int main(void)
{
	char str[] = "0.2,0.4,1";
	char *token = NULL, *savestr;
	double a1[3] = {0};
	int i = 0;
	token = strtok_r(str, ",", &savestr);
	while(token != NULL) {
		a1[i++] = atof(token);
		token = strtok_r(NULL, ",", &savestr);
	}

	double a2[] = {0.3, 0.5, 0.6, 1};

	double *res;
	res = convolution(a2, 4, a1, 3);

	free(res);
	return 0;
}
#endif
