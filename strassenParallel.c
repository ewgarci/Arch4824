/*****************************************************************
 * Adapted from http://web.cs.sunyit.edu/~levyt/strassen_pthread
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "job_api.h"
#define BREAK 32

matrix new_matrix(int size){
	matrix ret;
	ret.size = size;
	ret.rows = malloc(sizeof(double*)*size);
	ret.data = calloc(size*size, sizeof(double));
	int i;
	for(i=0;i<size;i++){
		ret.rows[i] = ret.data + i*size;
	}
	return ret;
}

matrix *adjust_matrix(int size, double **mem){
	int i;
	matrix *ret = malloc(sizeof(ret));
	ret->size = size;
	ret->data = mem[0];
	ret->rows = malloc(size * sizeof(double *));
	for(i=0;i<size;i++){
		ret->rows[i] = mem[0] + i*size;
	}
		
	return ret;
}

void delete_matrix(matrix m){
	free(m.rows);
	if(m.data != NULL){
		free(m.data);
	}
}

matrix get00(matrix m){
	matrix ret;
	ret.size = m.size/2;
	ret.rows = malloc(sizeof(double*)*ret.size);
	int i;
	for(i=0;i<ret.size;i++){
		ret.rows[i] = m.rows[i];
	}
	ret.data = NULL;
	return ret;
}

matrix get01(matrix m){
	matrix ret;
	ret.size = m.size/2;
	ret.rows = malloc(sizeof(double*)*ret.size);
	int i;
	for(i=0;i<ret.size;i++){
		ret.rows[i] = m.rows[i] + ret.size;
	}
	ret.data = NULL;
	return ret;
}

matrix get10(matrix m){
	matrix ret;
	ret.size = m.size/2;
	ret.rows = malloc(sizeof(double*)*ret.size);
	int i;
	for(i=0;i<ret.size;i++){
		ret.rows[i] = m.rows[i+ret.size];
	}
	ret.data = NULL;
	return ret;
}

matrix get11(matrix m){
	matrix ret;
	ret.size = m.size/2;
	ret.rows = malloc(sizeof(double*)*ret.size);
	int i;
	for(i=0;i<ret.size;i++){
		ret.rows[i] = m.rows[i+ret.size] + ret.size;
	}
	ret.data = NULL;
	return ret;
}

void plus(matrix a, matrix b, matrix c){
	int i,j;
	for(i=0;i<a.size;i++){
		for(j=0;j<a.size;j++){
			c.rows[i][j] = a.rows[i][j] + b.rows[i][j];
		}
	}
}

void minus(matrix a, matrix b, matrix c){
	int i,j;
	for(i=0;i<a.size;i++){
		for(j=0;j<a.size;j++){
			c.rows[i][j] = a.rows[i][j] - b.rows[i][j];
		}
	}
}

//expects c to be initialized
void strassen_mult(matrix a, matrix b, matrix c){
	matrix quar[12];
	matrix m[7];
	matrix temp1, temp2;
	int size = a.size/2;
	int i,j;

//	if(size==0){
//		c.rows[0][0] = a.rows[0][0] * b.rows[0][0];
//		return;
//	}
	
	if (a.size <= BREAK) {
		int i, j, k;
		for (i = 0; i < (a.size); i++) {
			for (k = 0; k < (a.size); k++) {
				for (j = 0; j < (a.size); j++) {
					c.rows[i][j] += a.rows[i][k] * b.rows[k][j];
				}
			}
		}
        
    }
	

	quar[0]=get00(a);  quar[1]=get01(a);
	quar[2]=get10(a);  quar[3]=get11(a);

	quar[4]=get00(b);  quar[5]=get01(b);
	quar[6]=get10(b);  quar[7]=get11(b);

	quar[8]=get00(c);  quar[9]=get01(c);
	quar[10]=get10(c); quar[11]=get11(c);
	
	for(i=0;i<7;i++){
		m[i] = new_matrix(size);
	}
	
	if (a.size <= BREAK)
		goto group;
		
		
	temp1 = new_matrix(size);
	temp2 = new_matrix(size);

	

	plus(quar[0], quar[3], temp1);
	plus(quar[4], quar[7], temp2);
	strassen_mult(temp1,temp2,m[0]);

	plus(quar[2], quar[3], temp1);
	strassen_mult(temp1, quar[4], m[1]);

	minus(quar[5], quar[7], temp2);
	strassen_mult(quar[0], temp2, m[2]);

	minus(quar[6], quar[4], temp2);
	strassen_mult(quar[3], temp2, m[3]);

	plus(quar[0], quar[1], temp1);
	strassen_mult(temp1, quar[7], m[4]);

	minus(quar[2], quar[0], temp1);
	plus(quar[4], quar[5], temp2);
	strassen_mult(temp1, temp2, m[5]);

	minus(quar[1], quar[3], temp1);
	plus(quar[6], quar[7], temp2);
	strassen_mult(temp1, temp2, m[6]);

	delete_matrix(temp1);
	delete_matrix(temp2);


group:
	
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			quar[8].rows[i][j] = m[0].rows[i][j] + m[3].rows[i][j] - m[4].rows[i][j] + m[6].rows[i][j];
			quar[9].rows[i][j] = m[2].rows[i][j] + m[4].rows[i][j];
			quar[10].rows[i][j] = m[1].rows[i][j] + m[3].rows[i][j];
			quar[11].rows[i][j] = m[0].rows[i][j] - m[1].rows[i][j] + m[2].rows[i][j] + m[5].rows[i][j];
		}
	}


	for(i=0;i<12;i++){
		delete_matrix(quar[i]);
	}
	for(i=0;i<7;i++){
		delete_matrix(m[i]);
	}
}

static void thread_main(void *args){
	struct work_struct *s = args;
	int tid = s->tid;
	int size = s->a_1.size;
	//printf("Thread %d Starting\n", tid);

    switch(tid){
		case 0:{
			matrix temp1 = new_matrix(size);
			matrix temp2 = new_matrix(size);
			plus(s->a_1, s->a_2, temp1);
			plus(s->b_1, s->b_2, temp2);
			strassen_mult(temp1, temp2, s->c);
			delete_matrix(temp2);
			delete_matrix(temp1);
		}
		break;
		case 1:{
			matrix temp1 = new_matrix(size);
			plus(s->a_1, s->a_2, temp1);
			strassen_mult(temp1, s->b_1, s->c);
			delete_matrix(temp1);
		}
		break;
		case 2:{
			matrix temp2 = new_matrix(size);
			minus(s->b_1, s->b_2, temp2);
			strassen_mult(s->a_1, temp2, s->c);
			delete_matrix(temp2);
		}
		break;
		case 3:{
			matrix temp2 = new_matrix(size);
			minus(s->b_1, s->b_2, temp2);
			strassen_mult(s->a_1, temp2, s->c);
			delete_matrix(temp2);
		}
		break;
		case 4:{
			matrix temp1 = new_matrix(size);
			plus(s->a_1, s->a_2, temp1);
			strassen_mult(temp1, s->b_1, s->c);
			delete_matrix(temp1);
		}
		break;
		case 5:{
			matrix temp1 = new_matrix(size);
			matrix temp2 = new_matrix(size);
			minus(s->a_1, s->a_2, temp1);
			plus(s->b_1, s->b_2, temp2);
			strassen_mult(temp1, temp2, s->c);
			delete_matrix(temp2);
			delete_matrix(temp1);
		}
		break;
		case 6:{
			matrix temp1 = new_matrix(size);
			matrix temp2 = new_matrix(size);
			minus(s->a_1, s->a_2, temp1);
			plus(s->b_1, s->b_2, temp2);
			strassen_mult(temp1, temp2, s->c);
			delete_matrix(temp2);
			delete_matrix(temp1);
		}
		break;
	}
}

void strassen_8thread_mult(matrix a, matrix b, matrix c){
	matrix quar[12];
	matrix m[7];
	struct work_struct *jobs[7];
	int size = a.size/2;
	int i, j;
	
	if(size==0){
		c.rows[0][0] = a.rows[0][0] * b.rows[0][0];
		return;
	}

	/*
	if (a.size <= BREAK) {
		int i, j, k;
		for (i = 0; i < (a.size); i++) {
			for (k = 0; k < (a.size); k++) {
				for (j = 0; j < (a.size); j++) {
					c.rows[i][j] += a.rows[i][k] * b.rows[k][j];
				}
			}
		}
        return;
    }
	*/
	
	job_init();
	
	quar[0]=get00(a);  quar[1]=get01(a);
	quar[2]=get10(a);  quar[3]=get11(a);

	quar[4]=get00(b);  quar[5]=get01(b);
	quar[6]=get10(b);  quar[7]=get11(b);

	quar[8]=get00(c);  quar[9]=get01(c);
	quar[10]=get10(c); quar[11]=get11(c);

	for(i=0;i<7;i++){
		m[i] = new_matrix(size);
	}

	for (i = 0; i < 7; i++) {
		jobs[i] = malloc(sizeof(struct work_struct));
		if (jobs[i] == NULL) {
			perror(NULL);
			exit(1);
		}
	}
	
	jobs[0]->a_1 = quar[0];
	jobs[0]->a_2 = quar[3];
	jobs[0]->b_1 = quar[4];
	jobs[0]->b_2 = quar[7];
	jobs[0]->c = m[0];
	jobs[0]->tid = 0;
	job_create(thread_main, jobs[0], 0);


	jobs[1]->a_1 = quar[2];
	jobs[1]->a_2 = quar[3];
	jobs[1]->b_1 = quar[4];
	jobs[1]->c = m[1];
	jobs[1]->tid = 1;
	job_create(thread_main, jobs[1], 0);


	jobs[2]->a_1 = quar[0];
	jobs[2]->b_1 = quar[5];
	jobs[2]->b_2 = quar[7];
	jobs[2]->c = m[2];
	jobs[2]->tid = 2;
	job_create(thread_main, jobs[2], 0);


	jobs[3]->a_1 = quar[3];
	jobs[3]->b_1 = quar[6];
	jobs[3]->b_2 = quar[4];
	jobs[3]->c = m[3];
	jobs[3]->tid = 3;
	job_create(thread_main, jobs[3], 0);


	jobs[4]->a_1 = quar[0];
	jobs[4]->a_2 = quar[1];
	jobs[4]->b_1 = quar[7];
	jobs[4]->c = m[4];
	jobs[4]->tid = 4;
	job_create(thread_main, jobs[4], 0);


	jobs[5]->a_1 = quar[2];
	jobs[5]->a_2 = quar[0];
	jobs[5]->b_1 = quar[4];
	jobs[5]->b_2 = quar[5];
	jobs[5]->c = m[5];
	jobs[5]->tid = 5;
	job_create(thread_main, jobs[5], 0);


	jobs[6]->a_1 = quar[1];
	jobs[6]->a_2 = quar[3];
	jobs[6]->b_1 = quar[6];
	jobs[6]->b_2 = quar[7];
	jobs[6]->c = m[6];
	jobs[6]->tid = 6;
	job_create(thread_main, jobs[6], 0);


	for(i=0;i<7;i++){
		job_join(jobs[i]);
	}

	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			quar[8].rows[i][j] = m[0].rows[i][j] + m[3].rows[i][j] - m[4].rows[i][j] + m[6].rows[i][j];
			quar[9].rows[i][j] = m[2].rows[i][j] + m[4].rows[i][j];
			quar[10].rows[i][j] = m[1].rows[i][j] + m[3].rows[i][j];
			quar[11].rows[i][j] = m[0].rows[i][j] - m[1].rows[i][j] + m[2].rows[i][j] + m[5].rows[i][j];
		}
	}

	for(i=0;i<12;i++){
		delete_matrix(quar[i]);
	}
	for(i=0;i<7;i++){
		delete_matrix(m[i]);
	}
}
/*
void print_matrix(matrix m){
	int i,j;
	printf("{");
	for(i=0;i<m.size;++i,i<m.size && printf(",")){
		printf("{");
		for(j=0;j<m.size;j++,j<m.size && printf(",")){
			printf("%d", (int)m.rows[i][j]);
		}
		printf("}");
	}
	printf("}");
}
*/
void strassenParallel(double **a, double **b, double **c, int tam) {
	
	matrix * as = adjust_matrix(tam, a);
	matrix * bs = adjust_matrix(tam, b);
	matrix * cs = adjust_matrix(tam, c);


	/*print_matrix(a);
	printf("\n");
	print_matrix(b);
	printf("\n");*/

	strassen_8thread_mult(*as,*bs,*cs);

	/*
	print_matrix(as);
	printf("\n");
	print_matrix(bs);
	printf("\n");
	print_matrix(cs);
	printf("\n");
	*/
}
