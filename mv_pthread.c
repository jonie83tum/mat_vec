#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

int random_int(int max) {
	int i;
	while ((i = (int) (((double) rand()) / ((double) RAND_MAX) * (max + 1))) 
			>= max);
	return i;
}

double** make_matrix(int n, int m) {
	double** A;
	int i, k;
	A = (double**) malloc(n * sizeof(double*));
	for (i=0;i<n;i++) 
		A[i] = (double*) malloc(m * sizeof(double));
	for (i=0;i<n;i++) 
		for (k=0;k<m;k++) 
			A[i][k] = 0.0;
	return A;
}

void free_matrix(double** A, int n, int m) {
	int i;
	for (i=0;i<n; i++) 
		free(A[i]);
	free(A);
}

double **fill_random_matrix(double **A, int n, int m) {
	int i,k;
	for (i=0;i<n;i++) 
		for (k=0;k<m;k++) 
			A[i][k] = ((double) rand()) / ((double) RAND_MAX);
	return A;
}

void print_matrix(double** A, int n, int m) {
	int i,k;
	for (i=0; i<n; i++) {
		for (k=0; k<m; k++) {
			printf("%e ",A[i][k]);
		}
		printf("\n");
	}
	printf("\n\n");
}

double *multiply_matrix_vector(double **A, double *v, double *w, int n, int m) {
	int i, k;
	for (i=0;i<n;i++) 
		w[i] = 0.0;
	for (i=0;i<n;i++) {
		for(k=0;k<m;k++) {
			w[i] += A[i][k] * v[k];
		}
	}
	return w;
}
double *make_vector(int n) {
	return (double*) malloc(n * sizeof(double));
}

void free_vector(double *v) {
	free(v);
}
double *fill_random_vector(double *v, int n) {
	int i;
	for (i=0;i<n;i++) 
		v[i] = ((double) rand()) / ((double) RAND_MAX);
	return v;
}

void print_vector(double *v, int n) {
	int i;
	for (i=0;i<n;i++) 
		printf("%e\n",v[i]);
	printf("\n");
}
double norm(double *v, double *w, int n){
	// this function calculates the euclidean norm of v-w
	double sum = 0;
	int i;
	for (i=0; i<n; i++){
		sum = sum + pow(v[i]-w[i],2);
	}
	sum = sqrt(sum);
	return sum;
}

// create a struct because only one parameter can be passed to start_routine
typedef struct {
	double **matrix;
	double *vec;
	double *result;
	int dimension;
	int row;
}thread_parm_t;

void *multiply_matrix_vector_para(void *parm){
	// this function should be called by each thread and multiplies the ith row of A with v 
	thread_parm_t *p = (thread_parm_t *)parm;
	int k;
	double **A = p->matrix;
	double *v = p->vec;
	double *wp = p->result;
	int m = p->dimension;
	int i = p->row;
	wp[i]=0;

	for(k=0;k<m;k++) {

		wp[i] += A[i][k] * v[k];
	}
	return NULL;
}

int main(int argc, char** argv) {
	int n=20, i;
	long thread_count;
	double **A, *v, *ws, *wp, no;
	thread_parm_t *parm[n];
	pthread_t thread[n];
	A=make_matrix(n, n);
	v=make_vector(n);
	wp=make_vector(n);
	ws=make_vector(n);
	A=fill_random_matrix(A, n, n);
	v=fill_random_vector(v, n);
	// calculate A * v sequentially as a reference
	ws=multiply_matrix_vector(A, v, ws, n, n);

	for (i=0; i<n; i++){
	parm[i] = malloc(sizeof(thread_parm_t));
	parm[i]->matrix=A;
	parm[i]->vec=v;
	parm[i]->result=wp;
	parm[i]->dimension=n;
	parm[i]->row=i;
	}

	// create all threads
	for (i=0; i<n; i++){
	pthread_create(&thread[i], NULL, multiply_matrix_vector_para, (void *)parm[i]);
	}
	
	//wait for all threads to complete	
	for (i=0; i<n; i++){
	pthread_join(thread[i],NULL);
	}
	// calcualte norm of wp - ws
	no = norm(wp,ws,n);
	printf("norm=%f\n",no);

	return 0;

}
