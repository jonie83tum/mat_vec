#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#define M 400000

int thread_count; 
double dx = 0.1e-2;
double dt = 0.1e-6;
int timcnt = 250;
double sum = 0;
pthread_mutex_t mutex;

typedef struct{
	int rank;
	double *uc;
	double *un;
	double *sum;
}thread_parm_t;

void* uc_init(void *parm){
	thread_parm_t *p = (thread_parm_t *)parm;
	int i, start, end;
	int rank = p->rank;
	double *uc = p->uc;
	start = rank*M/thread_count;
	end = (rank + 1) *M/thread_count;
	for (i=start; i<end; i++){
		uc[i] = 1.0;
	}

	return NULL;
}
void* uc_r(void *parm){
	thread_parm_t *p = (thread_parm_t *)parm;
	int i, start, end;
	int rank = p->rank;
	double *uc = p->uc;
	double *un = p->un;
	double r = dt/(dx*dx);
	double r1 = 1-2*r;
	start = rank*M/thread_count;
	end = (rank + 1) *M/thread_count;
	if (rank == 0){
		start = 1;
	}
	for( i=start; i<end; i++ ) {
		un[i] = r*(uc[i-1]+uc[i+1])+r1*uc[i];
	}
	return NULL;
}
void* uc_un(void *parm){
	thread_parm_t *p = (thread_parm_t *)parm;
	int i, start, end;
	int rank = p->rank;
	double *uc = p->uc;
	double *un = p->un;
	start = rank*M/thread_count;
	end = (rank + 1) *M/thread_count;
	if (rank == 0){
		start = 1;
	}
	for( i=start; i<end; i++ ) {
		uc[i] = un[i];
	}
	return NULL;
}
void* sum_1(void *parm){
	thread_parm_t *p = (thread_parm_t *)parm;
	int i, start, end;
	int rank = p->rank;
	double *uc = p->uc;
	double local_sum = 0;
	start = rank*M/thread_count;
	end = (rank + 1) *M/thread_count;
	if (rank == 0){
		start = 1;
	}
	if (rank == thread_count -1){
		end = M-1;
	}
	for( i=start; i<end; i++ ) {
		local_sum += uc[i];
	}
	pthread_mutex_lock(&mutex);
	sum += local_sum;
	pthread_mutex_unlock(&mutex); 

	return NULL;
}

int main(int argc, char **argv) {
	double uc[M+1], un[M+1];
	int i,j;

	thread_count = strtol(argv[1], NULL, 10);  
	if (M % thread_count != 0){
		printf("M (400000) has to be a multiple of the number of threads\n");
		return 0;
	}

	thread_parm_t *parm[thread_count];
	pthread_t thread[thread_count];

	for (i=0; i<thread_count; i++){
		parm[i] = (thread_parm_t*) malloc(sizeof(thread_parm_t));
		parm[i]->rank = i;
		parm[i]->uc = uc;
		parm[i]->un = un;
	}


	for (i=0; i<thread_count; i++){
		pthread_create(&thread[i], NULL, uc_init, (void *)parm[i]);
	}
	for (i=0; i<thread_count; i++){
		pthread_join(thread[i],NULL);
	}

	uc[0] = 0.0;
	uc[M] = 0.0;

	for( j=0; j<timcnt; j++ ) {
		for (i=0; i<thread_count; i++){
			pthread_create(&thread[i], NULL, uc_r, (void *)parm[i]);
		}

		for (i=0; i<thread_count; i++){
			pthread_join(thread[i],NULL);
		}
		for (i=0; i<thread_count; i++){
			pthread_create(&thread[i], NULL, uc_un, (void *)parm[i]);
		}

		for (i=0; i<thread_count; i++){
			pthread_join(thread[i],NULL);
		}

	}

	pthread_mutex_init(&mutex, NULL);
	for (i=0; i<thread_count; i++){
		pthread_create(&thread[i], NULL, sum_1, (void *)parm[i]);
	}

	for (i=0; i<thread_count; i++){
		pthread_join(thread[i],NULL);
	}
	printf( "sum  %9.8f\n", sum );

	return 0;
}



