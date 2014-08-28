#ifndef _COMMON_H
#define _COMMON_H

#define _expect_false(expr)		__builtin_expect(!!(expr), 0)
#define _expect_true(expr)		__builtin_expect(!!(expr), 1)
#define expect_false(cond)		_expect_false(cond)
#define expect_true(cond)		_expect_true(cond)

#define set_attribute(attrlist)		__attribute__(attrlist)
#define declare_noinline		set_attribute((__noinline__))
#define declare_unused			set_attribute ((__unused__))

#define prefetch(x)	__builtin_prefetch(x)
//#define offsetof(a,b)	__builtin_offsetof(a,b)

#define matrix(a, i, j, n)	((a) + (i) * (n) + (j))
#define m_weight(a, i, j, n)	(((MATRIX *)((a) + (i) * (n) + (j)))->weight)
#define m_parent(a, i, j, n)	(((MATRIX *)((a) + (i) * (n) + (j)))->parent)
#define m_path(a, i, j, n)	(((MATRIX *)((a) + (i) * (n) + (j)))->path)


#define array_needsize(type, base, cur, cnt, init)	\
	if(expect_false((cnt) > (cur))) {	\
		size_t ocur_ = (cur);	\
		(base) = (type *)array_realloc(sizeof(type), (base), &(cur), (cnt));	\
		init((base), (ocur_), (cur), sizeof(type));	\
	}

#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))

#define min(x,y) ({ \
	typeof(x) _x = (x);	\
	typeof(y) _y = (y);	\
	(void) (&_x == &_y);	\
	_x < _y ? _x : _y; })

#define max(x,y) ({ \
	typeof(x) _x = (x);	\
	typeof(y) _y = (y);	\
	(void) (&_x == &_y);	\
	_x > _y ? _x : _y; })

/*Durstenfeld's method*/
#define decl_shuffle(type)					\
void shuffle_##type(type *list, size_t len) {		\
	int j;									\
	type tmp;							\
	while(len) {							\
		j = irand(len);						\
		if (j != len - 1) {					\
			tmp = list[j];					\
			list[j] = list[len - 1];			\
			list[len - 1] = tmp;				\
		}								\
		len--;							\
	}									\
}

void * declare_noinline array_realloc(size_t elem, void *base, size_t *cur, size_t cnt);
void array_zero_init(void *p, size_t op, size_t np, size_t elem);
char *cmd_system(const char *cmd);
void double_to_string(char *str, const double *array, int len);
void int_to_string(char *str, const int *array, int len);
int irand(int n);

#endif

