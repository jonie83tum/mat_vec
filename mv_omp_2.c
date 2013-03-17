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
void multiply_matrix_vector_para(int n, double A[n][n], double *v, double *w, int l){
	int k, i, th_n;
	th_n = omp_get_thread_num();
	int start = th_n*l;
	for (i = start; i < start + l; i++){
		w[i]=0.0;
		for(k=0;k<n;k++) {
			w[i] += A[i][k] * v[k];
		}
	}
	return;
}
int main(int argc, char** argv) {
	double *v, *ws, *wp, no;
	int l;
	int n = atoi(argv[1]);
        int num_th = atoi(argv[2]);
	// check whether command line arguments are correct
        if (n < num_th){
                printf("num_th has to be smaller or equal than n\n");
                return 1;
        }
        if (n % num_th != 0){
                printf("n has to be a multiple of num_th\n");
                return 1;
        }
	l = n/num_th;
 	// allocate memory for A according to C99 standard (variabel size array)
        double (*A)[n] = malloc(n*n*sizeof(double));
        v = (double*) malloc(n*sizeof(double));
        wp = (double*) malloc(n*sizeof(double));
        ws = (double*) malloc(n*sizeof(double));
        // fill A and v with random numbers
        fill_random_matrix(n, n, A);
        fill_random_vector(v, n);
        multiply_matrix_vector(n, A, v, ws);
#	pragma omp parallel num_threads(num_th)	
	multiply_matrix_vector_para(n, A, v, wp, l);
	no = norm(ws,wp,n);
	printf("norm=%f\n",no);
	free(v);
        free(ws);
        free(wp);
        free(A);
	return 0;
}
