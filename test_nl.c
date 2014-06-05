#include <math.h>
#include <nlopt.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
	double a, b;
} my_constraint_data;

static double coe = 3.0;

double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
	if(grad) {
		grad[0] = 0.;
		grad[1] = 0.5 / sqrt(x[1]);
	}

	return sqrt(x[1]) * coe;
}

double myconstraint(unsigned n, const double *x, double *grad, void *data)
{
	my_constraint_data *d = (my_constraint_data *)data;
	double a = d->a, b = d->b;
	if(grad) {
		grad[0] = 3 * a * (a * x[0] + b) * (a * x[0] + b);
		grad[1] = -1.;
	}

	return ((a * x[0] + b) * (a * x[0] + b) * (a * x[0] + b) - x[1]);
}

int main(void)
{
	double lb[2] = {-HUGE_VAL, 0};
	double x[2] = {1.2, 5.4};
	double minf;
	my_constraint_data data[2] = {{2,0}, {-1,1}};
	nlopt_opt opt;

//	opt = nlopt_create(NLOPT_LD_MMA, 2);
	opt = nlopt_create(NLOPT_LN_COBYLA, 2);
	nlopt_set_lower_bounds(opt, lb);
	int test = 5;
	nlopt_set_min_objective(opt, myfunc, &test);

	nlopt_add_inequality_constraint(opt, myconstraint, &data[0], 1e-8);
	nlopt_add_inequality_constraint(opt, myconstraint, &data[1], 1e-8);

	nlopt_set_xtol_rel(opt, 1e-4);
	if(nlopt_optimize(opt, x, &minf) < 0)
		printf("nlopt failed\n");
	else
		printf("%lf, %lf-->min: %lf\n", x[0], x[1], minf);

	nlopt_destroy(opt);
	return 0;
}

