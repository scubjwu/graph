#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "uthash.h"
#include "dict.h"
#include "common.h"

int run_other(FILE *fres, int i, int cn, const char *s1, const char *s2, DICT *d)
{
	int m, n;
	double v[4][cn];
	char dst_dir[32] = {0};
	char file_name[64] = {0};
	char cmd[128] = {0};
	sprintf(dst_dir, "./sim%d/res_%s/", i, s1);
	sprintf(cmd, "cat %s/sim.conf | grep %s | cut -d ' ' -f 3", dst_dir, dict_get(d, s1));
	sprintf(file_name, "%s/%s.log", dst_dir, s2);
	
	FILE *f = fopen(file_name, "r");
	if(f == NULL) {
		perror("no such log file");
		return -1;
	}
	
	int sim_t = atoi(cmd_system(cmd));
	int step = -1, id = 0;
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	while((read = getline(&line, &len, f)) != -1) {
		double value;
		if(strstr(line, "###")) {
			step++;
			id = 0;
			continue;
		}
		if(strcmp(line, "\n") == 0)
			continue;
	
		sscanf(line, "%lf", &value);
		v[step][id] = value;
		id++;
	}
	fclose(f);
	free(line);
			
	//wirte back v[][]
	for(m=0; m<cn; m++) {
		fprintf(fres, "%d,%d,", sim_t, m);
		for(n=0; n<3; n++) 
			fprintf(fres, "%lf,", v[n][m]); 
				
		fprintf(fres, "%lf\r\n", v[n][m]);
	}

	return 0;
}

int get_node_num(const char *str)
{
	char tmp_name[64] = {0};
	sprintf(tmp_name, "./sim0/res_%s/dst.log", str);
	FILE *f = fopen(tmp_name, "r");
	if(f == NULL) 
		return -1;

	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	int res = 0;
	
	while((read = getline(&line, &len, f)) != -1) {
		if(strstr(line, "###") && res)
			break;
		res++;
	}
	fclose(f);
	free(line);

	return res - 2;
}

int run_src(FILE *fres, int i, const char *s1, const char *s2, DICT *d)
{
	char dst_dir[32] = {0};
	char file_name[64] = {0};
	char cmd[128] = {0};
	sprintf(dst_dir, "./sim%d/res_%s/", i, s1);
	sprintf(cmd, "cat %s/sim.conf | grep %s | cut -d ' ' -f 3", dst_dir, dict_get(d, s1));
	sprintf(file_name, "%s/%s.log", dst_dir, s2);

	FILE *f = fopen(file_name, "r");
	if(f == NULL) {
		perror("no such log file");
		return -1;
	}

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
			fprintf(fres, "%d,%d,", sim_t, src_id);
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

	return 0;
}

int main(int argc, char *argv[])
{
	if(argc < 3) {
		printf("lack of parameters.\nargv1: sim type\nargv2: sim obj\n");
		return 1;
	}

	dict_construct(sim_conf);
	dict_put(&sim_conf, "TW", "TIME_WINDOW");
	dict_put(&sim_conf, "CN", "CAN_NUM");
	dict_put(&sim_conf, "COST", "PRICE");
	dict_put(&sim_conf, "OB", "OB_WINDOW");
	dict_put(&sim_conf, "DRATIO", "DRATIO");

	char res_name[32] = {0};
	sprintf(res_name, "%s_%s.csv", dict_get(sim_conf, argv[1]), argv[2]);
	FILE *fres = fopen(res_name, "w");
	if(strcmp(argv[2], "can") == 0)
		fprintf(fres, "sim_t,id,c_comm_load,c_storage_load,d_comm_load,d_storage_load\r\n");
	else if(strcmp(argv[2], "dst") == 0)
		fprintf(fres, "sim_t,id,c_delay,c_receivings,d_delay,d_receivings\r\n");
	else if(strcmp(argv[2], "src") == 0)
		fprintf(fres, "sim_t,id,m_rev,c_rev,c_sharings,c_delay,d_success,d_rev,d_sharings,d_delay\r\n");
	else {
		printf("sim obj: can | dst | src\n");
		goto END;
	}

	int cn = get_node_num(argv[1]);
	if(cn == -1) {
		printf("sim type: TW | CN | COST | OB | DRATIO\n");
		goto END;
	}

	int i;
	for(i=0; i<10; i++) {
		if(strcmp(argv[2], "src") == 0) {
			if(run_src(fres, i, argv[1], argv[2], sim_conf) == -1)
				break;
		}
		else {
			if(run_other(fres, i, cn, argv[1], argv[2], sim_conf) == -1)
				break;
		}
	}

END:
	fclose(fres);
	dict_destroy(sim_conf);

	return 0;
}

