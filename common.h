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
#define offsetof(a,b)	__builtin_offsetof(a,b)

#define matrix(array, i, j, n)	(*(array + i * n + j))
#define m_weight(a, i, j, n)	(((MATRIX *)(a + i * n + j))->weight)
#define m_parent(a, i, j, n)	(((MATRIX *)(a + i * n + j))->parent)
#define m_path(a, i, j, n)	(((MATRIX *)(a + i * n + j))->path)


#define array_needsize(type, base, cur, cnt, init)	\
	if(expect_false((cnt) > (cur))) {	\
		size_t ocur_ = (cur);	\
		(base) = (type *)array_realloc(sizeof(type), (base), &(cur), (cnt));	\
		init((base), (ocur_), (cur), sizeof(type));	\
	}

void * declare_noinline array_realloc(size_t elem, void *base, size_t *cur, size_t cnt);
void array_zero_init(void *p, size_t op, size_t np, size_t elem);
char *cmd_system(const char *cmd);
void double_to_string(char *str, const double *array, int len);
void int_to_string(char *str, const int *array, int len);

#endif

