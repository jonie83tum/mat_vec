#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#define f(A) (4.0/(1.0+A*A))
int n, flag, thread_count;
double sum, w;

typedef struct {
	int rank;
}thread_parm_t;

void* thread_sum(void *parm) {
	thread_parm_t *p = (thread_parm_t *) parm;
	int rank = p->rank;
	int start = rank*n/thread_count;
	int end = (rank+1)*n/thread_count;
	double local_sum = 0.0, x;
	int i;

	w=(double)1.0/(double)n;
	for (i=start; i<end; i++){
		x = w*((double)i-0.5);
		local_sum += f(x);
	}

	while (flag != rank);
	sum += local_sum;  
	flag++;

	return NULL;
}

int main(int argc, char* argv[]) {
	int i;
	double pi;
	thread_count = strtol(argv[1], NULL, 10);  
	n = strtol(argv[2], NULL, 10);
	pthread_t* thread_handles;
	thread_parm_t *parm[thread_count];
	for (i=0; i<thread_count; i++){
		parm[i] = malloc(sizeof(thread_parm_t));
		parm[i]->rank = i;
	}

	thread_handles = (pthread_t*) malloc (thread_count*sizeof(pthread_t)); 

	sum = 0.0;
	flag = 0;

	for (i = 0; i < thread_count; i++){  
		pthread_create(&thread_handles[i], NULL, thread_sum, (void*)parm[i]);  
	}

	for (i = 0; i < thread_count; i++) 
		pthread_join(thread_handles[i], NULL); 
	pi = w*sum;
	printf("pi = %f \n", pi);	

	free(thread_handles);
	return 0;
}
double Serial_pi(long long n) {
	double sum = 0.0;
	long long i;
	double factor = 1.0;

	for (i = 0; i < n; i++, factor = -factor) {
		sum += factor/(2*i+1);
	}
	return 4.0*sum;

} 
