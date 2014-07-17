#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(void)
{
	FILE *f_path = fopen("./path.csv", "r");
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	int id = -1;
	char filename[128] = {0};
	char *buff = (char *)calloc(1024, sizeof(char));
	char *_buff = (char *)calloc(1024, sizeof(char));
	FILE *f = NULL;
	int mlen = strlen("\r\n,");

	while((read = getline(&line, &len, f_path)) != -1) {
		int i;
		sscanf(line, "%d,", &i);
		if(id != i) {
			if(f)
				fclose(f);
			
			filename[0] = 0;
			sprintf(filename, "%d.path", i);
			f = fopen(filename, "w");

			id = i;
		}

		char *token = NULL, *ptr;
		char dest[12] = {0};
		token = strtok(line, ",");
		int hop = 0;
		ptr = _buff;
		ptr[0] = 0;
		while(token != NULL) {
			int tmp = atoi(token);
			if(tmp != id) {
				ptr += sprintf(ptr, "%s,", token);
				hop++;
			}

			int l = strlen(token);
			int len = (l - mlen + 1) > 0 ? (l - mlen + 1) : l;

			memset(dest, 0, 12);
			memcpy(dest, token, len);

			token = strtok(NULL, ",");
		}

		memset(buff, 0, 1024);
		if(hop == 1) {
			ptr -= mlen;
			sprintf(ptr, ",self\r\n");
			memcpy(buff, _buff, strlen(_buff));
		}
		else {
			ptr = ptr - strlen(dest) - 1 - mlen;
			sprintf(ptr, "\r\n");
			sprintf(buff, "%s,%s", dest, _buff);
		}

		fwrite(buff, sizeof(char), strlen(buff), f);
	}
}

