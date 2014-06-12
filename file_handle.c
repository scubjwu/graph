#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>

#include "common.h"

#define RANGE	90000

static const unsigned int seed = INT_MAX/2 - 1;

int main(void)
{
	cmd_system("cp ./graph.csv.bak ./graph.csv");

	FILE *fin = fopen("./graph.csv", "r");
	FILE *fout = fopen("./new_graph.csv", "w");

	int total_time = atoi(cmd_system("tail -1 ./mobility.csv | cut -d , -f 3"));
	total_time = total_time/4;

	char *line = NULL;
	size_t len = 0;
	ssize_t read;

	srand(seed);

	while((read = getline(&line, &len, fin)) != -1) {
		int n, i, j, time, rtime;
		sscanf(line, "%d,%d,%d", &i, &j, &time);
		for(n=0; n<3; n++) {
			rtime = (total_time + rand()%RANGE) * n;
			time += rtime;
			fprintf(fout, "%d,%d,%d\r\n", i, j, time);
		}
	}

	fclose(fin);
	fclose(fout);
	free(line);

	cmd_system("sort -t, -k1,1n -k2,2n -k3,3n ./new_graph.csv -o ./new_graph.csv.bak; mv ./new_graph.csv.bak ./new_graph.csv");
	cmd_system("sort -t, -k3,3n ./new_graph.csv -o ./new_mobility.csv");

	return 0;
}

