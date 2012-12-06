#ifndef _SOMMELIER_H_
#define _SOMMELIER_H_

#include "smat.h"

void smat_mult(const struct smat *a, const struct smat *b, struct smat *r);
void smat_add(const struct smat *a, const struct smat *b, struct smat *r);
void smat_vect(const struct smat *a, const struct smat *v, struct smat *r);
void smat_scale(const struct smat *a, const struct smat *alpha, struct smat *r);

#endif /* _SOMMELIER_H_ */
