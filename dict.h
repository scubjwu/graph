#ifndef _DICT_H
#define _DICT_H

typedef struct dict_t {
	char *key;
	char *value;
	UT_hash_handle hh;
} DICT;

#define dict_construct(name)	\
	DICT *(name) = NULL;

void dict_destroy(DICT *dict);
void dict_put(DICT **dict, const char *key, const char *value);
char *dict_get(DICT *dict, const char *key);
void dict_del(DICT *dict, const char *key);

#endif
