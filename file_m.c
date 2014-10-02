#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "common.h"

int main(void)
{
	FILE *fin = fopen("./mobility.csv", "r");
	FILE *fout = fopen("./new_mobility.csv", "w");

	int total_time = atoi(cmd_system("tail -1 ./mobility.csv | cut -d , -f 3"));
	char *line = NULL;
	size_t len = 0;
	ssize_t read;

	while((read = getline(&line, &len, fin)) != -1) {
		int n, i, j, time;
		sscanf(line, "%d,%d,%d", &i, &j, &time);
		for(n=0; n<4; n++)
			fprintf(fout, "%d,%d,%d\r\n", i-1, j-1, time + n * total_time);
	}

	fclose(fin);
	fclose(fout);
	free(line);

	cmd_system("sort -t, -k3,3n ./new_mobility.csv -o ./new_mobility_t.csv");
	cmd_system("sort -t, -k1,1n -k2,2n -k3,3n ./new_mobility.csv -o ./new_graph.csv");
	cmd_system("mv ./new_graph.csv ./graph.csv; mv ./new_mobility_t.csv ./mobility.csv");

	return 0;
}
