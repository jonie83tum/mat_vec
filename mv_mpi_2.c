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
	int l, rank, size, k, i;
	double *v, *ws, *wp, no;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
/*      Matrix A and vector v are to be muliplied
        The solution is ws for sequential multiplication and
        wp for parallel multiplication 
*/
	int n = atoi(argv[1]);
        // check whether command line argument is correct
        if (n < size){
		if (rank == 0){
		printf("Number of processors must be smaller or qual than n\n");
		}
		MPI_Finalize();
		return 1;
        }
        if (n % size != 0){
		if (rank == 0){
		printf("n has to be a multiple of number of processors\n");
		}
		MPI_Finalize();
		return 1;
        }
        l = n/size;

        // allocate memory for A according to C99 standard (variabel size array)
	// allocate only the required memory
	double (*A)[n];
	if (rank == 0){
	        A = malloc(n*n*sizeof(double));
       		v = (double*) malloc(n*sizeof(double));
       		wp = (double*) malloc(n*sizeof(double));
        	ws = (double*) malloc(n*sizeof(double));
	}
	if (rank != 0){
	        A = malloc(l*n*sizeof(double));
       		v = (double*) malloc(n*sizeof(double));
       		wp = (double*) malloc(l*sizeof(double));
	}
        // fill and multiply A and v with random numbers only on rank 0
	if (rank == 0){
	        fill_random_matrix(n, n, A);
       		fill_random_vector(v, n);
        	multiply_matrix_vector(n, A, v, ws);
		// send A and v to the other nodes
		for (i=1; i<size; i++){
			for (k=0; k<l; k++){
				MPI_Send(&A[i*l+k][0], n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}	
			MPI_Send(&v[0], n, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
		}
	
	}
	// receive A and v from the other nodes
	if (rank != 0){
			for (k=0; k<l; k++){
				MPI_Recv(&A[k][0], n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}	
			MPI_Recv(&v[0], n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	// calculate the (local) values of w_p
	for (i=0; i<l; i++){
		wp[i] = 0;
		for(k=0; k<n; k++) {
			wp[i] += A[i][k] * v[k];
		}	
	}
	// rank 0 receives from all other processors and calculates the norm
	if (rank == 0){
		for (i=1; i<size; i++){
			MPI_Recv(&wp[i*l], l, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		no = norm(ws,wp,n);
		printf("norm=%f\n",no);
	}
	// All other processors have to send wp to rank 0
	if (rank != 0){
		MPI_Send(&wp[0], l, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	// Alternative solution would be using MPI_Gather or MPI_Reduce
	free(v);
	free(wp);
	free(A);
	if (rank == 0){
		free(ws);
	}
	MPI_Finalize();
	return 0;
}
