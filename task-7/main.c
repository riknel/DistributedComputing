#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <pthread.h>
#include <assert.h>

#define L 0
#define R 1
#define U 2
#define D 3
#define OK 4

typedef unsigned int uint;

typedef struct particle {
    //координаты внутри текущего квадрата
    int x;
    int y;
    int r;
} particle;

typedef struct arguments{
    int rang;
    int count_proc;
    int l;
    int a;
    int b;
    int N;
}arguments;

void walk(int rang, int count_proc, int l, int a, int b, int N);

void* walk_void_args(void* args_){
    arguments* args = (arguments*) args_;
    walk(args->rang, args->count_proc, args->l, args->a, args->b, args->N);
    return 0;
}


void walk(int rang, int count_proc, int l, int a, int b, int N) {
    double start = MPI_Wtime();

    //определяем квадрат процесса с номером rang
    int rang_x = rang % a;
    int rang_y = rang / a;

    int particles_count = N;
    particle* particles = (particle*) malloc(sizeof(particle) * particles_count);

    //для использовония rand_r генерируем seeds из главного процесса и отправляем по одному seed в остальные процессы
    uint seed;
    uint* seeds = (uint*)malloc(sizeof(uint)* count_proc);
    if (rang == 0) {
        srand(time(NULL));
        for (int i = 0; i < count_proc; ++i) {
            seeds[i] = rand();
        }
    }
    MPI_Scatter(seeds, 1, MPI_UNSIGNED, &seed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    free(seeds);

    //генерируем N частиц в квадрате
    for(int i = 0; i < N; ++i) {
        particles[i].x = rand_r(&seed) % l;
        particles[i].y = rand_r(&seed) % l;
        particles[i].r = rand_r(&seed) % (a * b);
    }

//    если завершились, то каждый процесс для каждой точки из своего квадрата подсчитывает нужный вектор
    int** mask = (int**)malloc(l*sizeof(int*));
    for(int i = 0; i < l; ++i) {
        mask[i] = (int*)malloc(l*N*sizeof(int));
        for(int j = 0; j < l*N; ++j) {
            mask[i][j] = 0;
        }
    }

    MPI_File f_bin;
    MPI_File_delete("data.bin", MPI_INFO_NULL);
    MPI_File_open(MPI_COMM_WORLD, "data.bin", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &f_bin);

    //проходимся по массиву частиц
    for(int i = 0; i < N; ++i) {
        int part_x = particles[i].x;
        int part_y = particles[i].y;
        int part_r = particles[i].r;
        mask[part_y][part_x * count_proc + part_r] += 1;
    }

    //теперь каждый процесс выводит в файл слепок своей области
    for (int i = 0; i < l; ++i) {
        for(int j = 0; j < l; ++j) {
            //выразим координаты точки на торе через координаты его узла и координаты точки в узле
            int general_x = rang_x * l + j;
            int general_y = rang_y * l + i;
            int line_length = l * N * a;

            //находим насколько нужно сдвитнуь каретку
            int place = general_y * line_length + general_x * N;

            MPI_File_set_view(f_bin, place*sizeof(int), MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
            MPI_File_write(f_bin , &mask[i][j * N], N,  MPI_INT, MPI_STATUS_IGNORE);
        }

    }

    MPI_File_close(&f_bin);

    for (int i=0;i < l; i++) {
        free(mask[i]);
    }
    free(mask);

    //заканчиваем подсчет времени
    double end = MPI_Wtime();

    // главный процесс выводит в файл входные данные и время работы программы
    if (rang==0) {
        double program_time = end - start;
        FILE *f_stats;
        f_stats = fopen("stats.txt", "w");
        fprintf(f_stats, "%d %d %d %d %fs", l, a, b, N, program_time);
        fclose(f_stats);
    }

    free(particles);
}


int main(int argc, char * argv[]) {
    MPI_Init(&argc, &argv);

    int rang, count_proc;

    MPI_Comm_rank(MPI_COMM_WORLD, &rang);
    MPI_Comm_size(MPI_COMM_WORLD, &count_proc);

    int l, a, b, N;

    l = atoi(argv[1]);
    a = atoi(argv[2]);
    b = atoi(argv[3]);
    N = atoi(argv[4]);

    //кол-во процессов должно быть равно кол-ву узлов
    assert (count_proc == a * b);

    arguments args;
    args.rang = rang;
    args.count_proc = count_proc;
    args.l = l;
    args.a = a;
    args.b = b;
    args.N = N;

    walk_void_args(&args);

    MPI_Finalize();
    return 0;
}
