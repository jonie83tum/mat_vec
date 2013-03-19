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
	printf("Hello World from thread %d out of %d.\n", n, c);
	return NULL;
}

int main(int argc, char* argv[]){
	// parameter list 
	thread_parm_t *parm;
	// fill the parameter list
	parm = malloc (sizeof(thread_parm_t));
	parm->thread_num = 0;
	parm->thread_count = 1;
	hello(parm);	

	return 0;
}
