#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>

typedef struct {
	int cur_x;
	int timer;
	int is_end;
}particle;

typedef unsigned int uint;
omp_lock_t lock;

int step(int a, int b, float p, particle *part, uint seed){
	//для имитации выбора с вероятностью p будем выбирать рандомно числа до RAND_MAX и если число лежит
	//в первых p*RAND_MAX, то выбираем ход направо

	while(!part->is_end){
		//выбираем шаг направо
		if(((double) rand_r(&seed)) <=  (double) RAND_MAX * p) {
			part->cur_x++;
		}
		//налево
		else{
			part->cur_x--;
		}
		part->timer++;

		//если дошли до завершающих состояний, возвращаем 1, если дошли до а, 2, если дошли до b
		if(part->cur_x == a){
			part->is_end = 1;
			return 1;
		}
		else if(part->cur_x == b){
			part->is_end = 1;
			return 2;
		}
	}
	return 0;
}




int main(int argc, char **argv){

	struct timeval begin, end;

	FILE *f = fopen("stats_python_16.txt", "w");

	int N, a, b, x, P;
	float p;
	int frequency_a = 0;
	int frequency_b = 0;

	a = atoi(argv[1]);
	b = atoi(argv[2]);
	x = atoi(argv[3]);
	N = atoi(argv[4]);
	p = atof(argv[5]);
	P = atoi(argv[6]);

	srand(time(NULL));

	omp_lock_t lock;
	omp_init_lock(&lock);

  for(int n=10; n <= 10000; n=n+50){

    particle* particles = (particle*)malloc(sizeof(particle)*n);
  	for(int i = 0; i < n; ++i){
  		particles[i].cur_x = x;
  		particles[i].timer = 0;
  		particles[i].is_end = 0;
  	}
		frequency_a = 0;
		frequency_b = 0;

		uint* seed = (uint*)malloc(sizeof(uint)*n);
		for(int i = 0; i < n; ++i) {
			seed[i] = rand();
		}

    assert(gettimeofday(&begin, NULL) == 0);
		for(int i =0; i < n; i++){
			int res = step(a, b, p, &particles[i], seed[i]);
			if(res == 1){
				frequency_a ++;
			}
			else if(res == 2){
				frequency_b ++;
			}
		}
    assert(gettimeofday(&end, NULL) == 0);
    double program_time_1thread  = ((end.tv_sec - begin.tv_sec) * 1000000u + end.tv_usec - begin.tv_usec) / 1.e6;

    for(int i = 0; i < n; ++i){
  		particles[i].cur_x = x;
  		particles[i].timer = 0;
  		particles[i].is_end = 0;
  	}
		frequency_a = 0;
		frequency_b = 0;

    omp_set_num_threads(P);
    assert(gettimeofday(&begin, NULL) == 0);
  	//захватывать lock будем только N раз, когда каждая частица завершила свое блуждание
		#pragma omp parallel for
		for(int i =0; i < n; i++){
			int res = step(a, b, p, &particles[i], seed[i]);
			if(res == 1){
				omp_set_lock(&lock);
				frequency_a ++;
				omp_unset_lock(&lock);
			}
			else if(res == 2){
				omp_set_lock(&lock);
				frequency_b ++;
				omp_unset_lock(&lock);
			}
		}
    assert(gettimeofday(&end, NULL) == 0);
    double program_time  = ((end.tv_sec - begin.tv_sec) * 1000000u + end.tv_usec - begin.tv_usec) / 1.e6;
    double S_p = program_time_1thread / program_time;
    double E_p = S_p / P;
    fprintf(f, "%d %f %f %f\n", n, program_time, S_p, E_p);
		free(particles);
		free(seed);

  }


	omp_destroy_lock(&lock);


	fclose(f);
	return 0;
}
