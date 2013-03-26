#include <stdio.h>
#include <stdlib.h>
#define f(A) (4.0/(1.0+A*A))
int n, thread_count;
double sum, w;

int main(int argc, char* argv[]) {
	int i;
	double pi, x;
	thread_count = strtol(argv[1], NULL, 10);  
	n = strtol(argv[2], NULL, 10);
	// check whether command line arguments are correct
	if (n % thread_count != 0){
		printf("n has to be a multiple of the number of threads\n");
		return 0;
	}

	sum = 0.0;
	
	w=(double)1.0/(double)n;
	// calculate the local sum
	for (i=0; i<n; i++){
		x = w*((double)i-0.5);
		sum += f(x);
	}


	pi = w*sum;
	printf("pi = %f \n", pi);
	
	return 0;
}
