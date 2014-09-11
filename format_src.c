#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "uthash.h"
#include "dict.h"
#include "common.h"

int main(int argc, char *argv[])
{
	dict_construct(sim_conf);
	dict_put(&sim_conf, "TW", "TIME_WINDOW");
	dict_put(&sim_conf, "CN", "CAN_NUM");
	dict_put(&sim_conf, "COST", "COST");
	dict_put(&sim_conf, "OB", "OB_WINDOW");

	char res_name[32] = {0};
	sprintf(res_name, "%s_sim.csv", dict_get(sim_conf, argv[1]));
	FILE *fres = fopen(res_name, "w");
	fprintf(fres, "sim_t,id,m_rev,c_rev,c_sharings,c_delay,d_success,d_rev,d_sharings,d_delay\r\n");

	int i;
	for(i=0; i<10; i++) {
		char dst_dir[32] = {0};
		char file_name[64] = {0};
		char cmd[128] = {0};
		sprintf(dst_dir, "./sim%d/res_%s/", i, argv[1]);
		sprintf(cmd, "cat %s/sim.conf | grep %s | cut -d ' ' -f 3", dst_dir, dict_get(sim_conf, argv[1]));
		sprintf(file_name, "%s/src.log", dst_dir);

		FILE *f = fopen(file_name, "r");
		int sim_t = atoi(cmd_system(cmd));
		int src_id, cnt = 1;
		char name[32] = {0};
		double tmp, value;
		char *line = NULL;
		size_t len = 0;
		ssize_t read;
		while((read = getline(&line, &len, f)) != -1) {
			if(cnt == 1) {
				sscanf(line, "%s %d", name, &src_id);
				fprintf(fres, "%s=%d,%d,", dict_get(sim_conf, argv[1]), sim_t, src_id);
				cnt++;
			}
			else {
				if(cnt == 3 || cnt ==4) {
					sscanf(line, "%s %lf", name, &tmp);
					cnt++;
				}
				else if(cnt == 8) {
					sscanf(line, "%s %lf", name, &value);
					fprintf(fres, "%lf,", value/tmp);
					if(value == 0) {
						fprintf(fres, "0,0,0\r\n");
						while((read = getline(&line, &len, f)) != -1) {
							cnt++;
							if(cnt == 12) {
								cnt = 1;
								break;
							}
						}
						continue;
					}
					cnt++;
				}
				else if(cnt == 11) {
					sscanf(line, "%s %lf", name, &value);
					fprintf(fres, "%lf", value);
					cnt++;
				}
				else if(cnt == 12) {
					fprintf(fres, "\r\n");
					cnt = 1;
				}
				else {
					sscanf(line, "%s %lf", name, &value);
					fprintf(fres, "%lf,", value);
					if(cnt == 2 && value == 0) {
						fprintf(fres, "0,0,0,0,0,0,0\r\n");
						getline(&line, &len, f);
						cnt = 1;
						continue;
					}
					cnt++;
				}
			}
		}
		fclose(f);
		free(line);
	}
	fclose(fres);
}
