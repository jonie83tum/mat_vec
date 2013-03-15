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

void fill_random_matrix(int n, int m, double A[n][m]) {
	int i,k;
	for (i=0;i<n;i++) 
		for (k=0;k<m;k++) 
			A[i][k] = ((double) rand()) / ((double) RAND_MAX);
	return;
}
void print_matrix(int n, int m, double A[n][m]) {
	int i,k;
	for (i=0; i<n; i++) {
		for (k=0; k<m; k++) {
			printf("%e ",A[i][k]);
		}
		printf("\n");
	}
	printf("\n\n");
}

void fill_random_vector(double *v, int n) {
	int i;
	for (i=0;i<n;i++) 
		v[i] = ((double) rand()) / ((double) RAND_MAX);
	return;
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
		sum += pow(v[i]-w[i],2);
	}
	sum = sqrt(sum);
	return sum;
}
void multiply_matrix_vector(int n, double A[n][n], double *v, double *w) {
	int i, k;
	for (i=0;i<n;i++) 
		w[i] = 0.0;
	for (i=0;i<n;i++) {
		for(k=0;k<n;k++) {
			w[i] += A[i][k] * v[k];
		}
	}
	return;
}
// create a struct because only one parameter can be passed to start_routine
typedef struct {
	double (*matrix)[];
	double *vec;
	double *result;
	int dimension;
	int row;
}thread_parm_t;

void *multiply_matrix_vector_para(void *parm){
	// this function is called by each thread and multiplies the ith row of A with v 
	thread_parm_t *p = (thread_parm_t *)parm;
	int k;
	int m = p->dimension;
	double (*A)[m] = p->matrix;
	double *v = p->vec;
	double *wp = p->result;
	int i = p->row;
	wp[i]=0;
	for(k=0;k<m;k++) {
		wp[i] += A[i][k] * v[k];
	}
	return NULL;
}
int main(int argc, char** argv) {
	int n=5, i;
	double *v, *ws, *wp, no;
	 // allocate memory for A according to C99 standard (variabel size array)
        double (*A)[n] = malloc(n*n*sizeof(double));
        v = (double*) malloc(n*sizeof(double));
        wp = (double*) malloc(n*sizeof(double));
        ws = (double*) malloc(n*sizeof(double));

	thread_parm_t *parm[n];
	pthread_t thread[n];
	fill_random_matrix(n, n, A);
	fill_random_vector(v, n);
	// calculate A * v sequentially as a reference
	multiply_matrix_vector(n, A, v, ws);
	// for each thread one parameter list is creadted
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
	free(v);
	free(ws);
	free(wp);
	free(A);
	return 0;
}
