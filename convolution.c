#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//complexity O(n1*n2)
double __attribute__((optimize("O0"))) *convolution(double *a1, int n1, double *a2, int n2)
{
	int i, j, n = n1 - 1 + n2, p = 0;

	double tmp[n][n1];
//	memset(tmp, 0, n * n1 * sizeof(double));
	bzero(tmp, n * n1 * sizeof(double));

	double *res = (double *)calloc(n, sizeof(double));

	for(i=0; i<n1; i++) {
		for(j=0; j<n2; j++) {
			tmp[j + i][i] = a1[i] * a2[j];
		}
	}

	for(i=0; i<n; i++) {
		double sum = 0;
		for(j=0; j<n1; j++) {
			sum += tmp[i][j];
		}
		res[i] = sum;
	}

	//the len of res equals to (n1 + n2 - 1)
	return res;
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
