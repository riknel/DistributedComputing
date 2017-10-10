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

int step(int a, int b, float p, particle *part){
	//для имитации выбора с вероятностью p будем выбирать рандомно числа до RAND_MAX и если число лежит
	//в первых p*RAND_MAX, то выбираем ход направо

	int rand_number = rand_r();

	//выбираем шаг направо
	if(rand_number < p * RAND_MAX) {
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
	if(part->cur_x == b){
		part->is_end = 1;
		return 2;
	}
	return 0;

}


int main(int argc, char **argv){

	srand(time(NULL));
	struct timeval begin, end;

	FILE *f = fopen("stats.txt", "a");

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

	particle* particles = (particle*)malloc(sizeof(particle)*N);
	for(int i = 0; i < N; ++i){
		particles[i].cur_x = x;
		particles[i].timer = 0;
		particles[i].is_end = 0;
	}

	omp_set_num_threads(P);
	omp_lock_t lock;
	omp_init_lock(&lock);

	//начинаем подсчет времени
	assert(gettimeofday(&begin, NULL) == 0);


	//захватывать lock будем только N раз, когда каждая частица завершила свое блуждание
	#pragma omp parallel for
	for(int i =0; i < N; ++i){
		while(!particles[i].is_end){
			int res = step(a, b, p, &particles[i]);
			if (res == 1) {
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
	}

	//заканчиваем подсчет времени
	assert(gettimeofday(&end, NULL) == 0);

	omp_destroy_lock(&lock);

	float mean_time_walk = 0;
	for(int i = 0; i < N; ++i){
		mean_time_walk += particles[i].timer;

	}
	mean_time_walk = mean_time_walk / N;
	double program_time = ((end.tv_sec - begin.tv_sec) * 1000000u + end.tv_usec - begin.tv_usec) / 1.e6;
	double probability_b = (double)frequency_b / N;

	fprintf(f, "%f, %f, %f, %d, %d, %d, %d, %f, %d\n", probability_b, mean_time_walk, program_time, a, b, x, N, p, P);
	fclose(f);
	return 0;
}
