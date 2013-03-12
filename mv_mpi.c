#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>
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
	int n=5, rank, size, k;
	double **A, *v, *ws, *wp, no, w_local=0;
	// Matrix A and vector v are to be muliplied
	// The solution is ws for sequential multiplication and
	// wp for parallel multiplication
	A=make_matrix(n, n);
	v=make_vector(n);
	ws=make_vector(n);
	wp=make_vector(n);
	A=fill_random_matrix(A, n, n);
	v=fill_random_vector(v, n);
	ws=multiply_matrix_vector(A, v, ws, n, n);
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (size != n){
		if (rank == 0){
		printf("Number of processors must be equal to length of vector\n");
		}
		MPI_Finalize();
		return 0;
	}
	for(k=0;k<n;k++) {
		w_local += A[rank][k] * v[k];
	}
	// rank 0 receives from all other processors and calculates the norm
	if (rank == 0){
		wp[0] = w_local;
		for (k=1; k<n; k++){
			MPI_Recv(&wp[k], 1, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		no = norm(ws,wp,n);
		printf("norm=%f\n",no);
	}
	// All other processors have to send w_local to p0
	if (rank != 0){
		MPI_Send(&w_local, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	// Alternative solution would be using  MPI_Gather
	free_vector(v);
	free_vector(ws);
	free_vector(wp);
	free_matrix(A, n, n);
	MPI_Finalize();
	return 0;
}
