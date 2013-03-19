#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

typedef struct{
	int thread_num;
	int thread_count;
}thread_parm_t;

void* hello(void* parm){
	thread_parm_t *p = (thread_parm_t *)parm;
	int n = p->thread_num;
	int c = p->thread_count;
	printf("Hello World form thread %d out of %d.\n", n, c);
	return NULL;
}

int main(int argc, char* argv[]){
	// define number of threads
	int n = 5, i;
	// parameter list for each thread
	thread_parm_t *parm[n];
	// thread handle for each thread
	pthread_t thread_handles[n];
	// fill the parameter list
	for (i=0; i<n; i++){
		parm[i] = malloc (sizeof(thread_parm_t));
		parm[i]->thread_num = i;
		parm[i]->thread_count = n;
	}
	// create all threads
	for (i=0; i<n; i++){
		pthread_create(&thread_handles[i], NULL, hello, (void *) parm[i]);
	}
	// wait for all threads to complete
	for (i=0; i<n; i++){
		pthread_join(thread_handles[i], NULL);
	}	
	return 0;
}
