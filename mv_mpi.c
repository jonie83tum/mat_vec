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
	return;
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
int main(int argc, char** argv) {
	MPI_Init(NULL, NULL);
	int n=5, rank, size, k;
	double *v, *ws, *wp, no, w_local=0;
/*      Matrix A and vector v are to be muliplied
        The solution is ws for sequential multiplication and
        wp for parallel multiplication 
*/
        // allocate memory for A according to C99 standard (variabel size array)
        double (*A)[n] = malloc(n*n*sizeof(double));
        v = (double*) malloc(n*sizeof(double));
        wp = (double*) malloc(n*sizeof(double));
        ws = (double*) malloc(n*sizeof(double));
        // fill A and v with random numbers
        fill_random_matrix(n, n, A);
        fill_random_vector(v, n);
        multiply_matrix_vector(n, A, v, ws);
	
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
	// Alternative solution would be using MPI_Gather or MPI_Reduce
	free(v);
	free(ws);
	free(wp);
	free(A);
	MPI_Finalize();
	return 0;
}
