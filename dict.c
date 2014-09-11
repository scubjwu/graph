#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "uthash.h"
#include "dict.h"

#define _free(x) {	\
	if(x)	free(x);	\
}

void dict_destroy(DICT *dict)
{
	DICT *cur, *tmp;

	HASH_ITER(hh, dict, cur, tmp) {
		HASH_DEL(dict, cur);
		_free(cur->key);
		_free(cur->value);
		_free(cur);
	}
}

void dict_put(DICT **dict, const char *key, const char *value)
{
	DICT *s;
	int value_len = strlen(value);
	HASH_FIND_STR(*dict, key, s);
	if(s == NULL) {
		int key_len;
		key_len = strlen(key);
		s = (DICT *)calloc(1, sizeof(DICT));
		s->key = (char *)calloc(key_len, sizeof(char));
		s->value = (char *)calloc(value_len, sizeof(char));
	
		strncpy(s->key, key, key_len);
		strncpy(s->value, value, value_len);
		HASH_ADD_STR(*dict, key, s);
		return;
	}

	if(value_len > strlen(s->value))
		s->value = (char *)realloc(s->value, value_len);
	memset(s->value, 0, strlen(s->value) * sizeof(char));
	strncpy(s->value, value, value_len);
}

char *dict_get(DICT *dict, const char *key)
{
	DICT *s;

	HASH_FIND_STR(dict, key, s);

	if(s)
		return s->value;
	else
		return NULL;
}

void dict_del(DICT *dict, const char *key)
{
	DICT *s;
	HASH_FIND_STR(dict, key, s);
	if(s) {
		HASH_DEL(dict, s);
		_free(s->value);
		_free(s->key);
		_free(s);
	}
}

#if 0
int main(void)
{
	dict_construct(test);

	dict_put(&test, "1", "test");
	dict_put(&test, "test", "pp");	
	dict_put(&test, "1", "l");
	dict_del(test, "test");

	printf("%s\t%s\n", dict_get(test, "test"), dict_get(test, "1"));

	dict_destroy(test);

	return 0;
}
#endif
