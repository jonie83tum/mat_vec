#include <stdio.h>
#include <stdlib.h>
#define f(A) (4.0/(1.0+A*A))
int n, thread_count;
double sum, w;

int main(int argc, char* argv[]) {
	int i;
	double pi, x;
	n = strtol(argv[1], NULL, 10);

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
