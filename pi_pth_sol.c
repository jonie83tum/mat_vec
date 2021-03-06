#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#define f(A) (4.0/(1.0+A*A))
int n, thread_count, flag;
double sum, w;
void* thread_sum(void *rank) {
	// cast the void pointer back into integer	
	int *ranki = (int *) rank;	
	int start = *ranki*n/thread_count;
	int end = (*ranki+1)*n/thread_count;
	double local_sum = 0.0, x;
	int i;

	w=(double)1.0/(double)n;
	// calculate the local sum
	for (i=start; i<end; i++){
		x = w*((double)i-0.5);
		local_sum += f(x);
	}
	// update the global sum (critical section) using busy waiting 
	while (flag != *ranki);
	sum += local_sum;  
	flag++;

	return NULL;
}

int main(int argc, char* argv[]) {
	int i;
	double pi;
	thread_count = strtol(argv[1], NULL, 10);  
	n = strtol(argv[2], NULL, 10);
	// check whether command line arguments are correct
	if (n % thread_count != 0){
		printf("n has to be a multiple of the number of threads\n");
		return 0;
	}

	pthread_t* thread_handles;
	// create an array in order to pass a thread number to each thread
	int *threads;
	threads = (int *)malloc(thread_count*sizeof(int));
	thread_handles = (pthread_t*) malloc (thread_count*sizeof(pthread_t)); 
	
	sum = 0.0;
	flag = 0;

	for (i = 0; i < thread_count; i++){  
		threads[i] = i;
		pthread_create(&thread_handles[i], NULL, thread_sum, (void*)&threads[i]);  
		// Attention: If <(void*)&i> is passed to the start_routine instead of 
		// <(void*)&threads[]> this result is a race condition!
	}

	for (i = 0; i < thread_count; i++){ 
		pthread_join(thread_handles[i], NULL);
	} 
	pi = w*sum;
	printf("pi = %f \n", pi);
	
	free(threads);
	free(thread_handles);
	return 0;
}
