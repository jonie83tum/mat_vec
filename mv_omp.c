#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <omp.h>

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
int main(int argc, char** argv) {
	int n;
	double **A, *v, *ws, *wp, no;

	n = 5;

	A=make_matrix(n, m);
	v=make_vector(n);
	ws=make_vector(n);
	wp=make_vector(n);
	A=fill_random_matrix(A, n, n);
	v=fill_random_vector(v, n);
	ws=multiply_matrix_vector(A, v, w, n, n);
	no = norm(ws,wp,n);
	printf("norm=%f\n",no);
	print_vector(v,n);
	print_matrix(A, n, n);
	print_vector(ws,n);

	return 0;
}
