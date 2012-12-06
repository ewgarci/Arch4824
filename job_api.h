#ifndef _JOB_API_H_
#define _JOB_API_H_

#ifdef SESC
#include "sescapi.h"
#else
#include <pthread.h>
#endif

typedef struct matrix {
	int size;
	double** rows;
	double* data;
} matrix;

struct work_struct {
	matrix a_1,a_2,b_1,b_2,c;
	int tid;
	int id;
	/*
	 * add here whatever fields you need. See job_api_example.[ch]
	 * for an example.
	 */
#ifndef SESC
	pthread_t pthread;
#endif /* !SESC */
};

#ifdef SESC

static inline void job_create(void (*func), void *arg, long flags)
{
	sesc_spawn(func, arg, flags);
}

static inline int job_join(struct work_struct *work)
{
	sesc_wait();
	return 0;
}

static inline void job_init(void)
{
	sesc_init();
}

static inline void job_exit(void)
{
	sesc_exit(0);
}

#else

static inline void job_create(void (*func), void *arg, long flags)
{
	struct work_struct *work = arg;

	pthread_create(&work->pthread, NULL, func, arg);
}

static inline int job_join(struct work_struct *work)
{
	return pthread_join(work->pthread, NULL);
}

static inline void job_init(void)
{
}

static inline void job_exit(void)
{
	pthread_exit(NULL);
}

#endif /* SESC */

#endif /* _JOB_API_H_ */
