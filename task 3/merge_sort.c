#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>
#include <pthread.h>

//кол-во доступных потоков
int count_threads = 0;

pthread_mutex_t mutex;

void merge_parallel(int* array, int* copy_array, int left_1, int right_1, int left_2, int right_2, int left_3, int m);
void sort_parallel(int* array, int* copy_array, int left, int right, int m);

//структура для хранения аргументов функций merge_parellel и sort_parallel
// в функцию будем передовать void* arguments
typedef struct {
    int *array;
    int *copy_array;
    int left_1;
    int right_1;
    int left_2;
    int right_2;
    int left_3;
    int m;
} arguments;


void swap(int* a, int* b) {
  int temp = *a;
  *a = *b;
  *b = temp;
}

// функция, находящая индекс первого элемента в отсортированном array[left:right - 1], который >= x
// если такого нет, то возращаем следующий элемент за подмассивом, то есть right
int binSearch(int* array, int x, int left, int right) {
    if(right <= left) {
        return right;
    }

    int middle = (left + right) / 2;
    if (array[middle] < x) {
        return binSearch(array, x, middle + 1, right);
    }
    return binSearch(array, x, left, middle);
}

//обычный однопоточный merge
void merge(int* array, int* copy_array, int left_1, int right_1, int left_2, int right_2, int left_3) {
  int index_1 = 0;
  int index_2 = 0;
  int index_3 = 0;
  while (left_1 + index_1 <= right_1 && left_2 + index_2 <= right_2) {
    if(array[left_1 + index_1] <= array[left_2 + index_2]) {
      copy_array[left_3 + index_3] = array[left_1 + index_1];
      index_1++;
    }
    else{
      copy_array[left_3 + index_3] = array[left_2 + index_2];
      index_2++;
    }
    index_3++;
  }
  while (left_1 + index_1 <= right_1) {
    copy_array[left_3 + index_3] = array[left_1 + index_1];
    index_1++;
    index_3++;
  }
  while (left_2 + index_2 <= right_2) {
    copy_array[left_3 + index_3] = array[left_2 + index_2];
    index_2++;
    index_3++;
  }
}

//вспомогательная функция для merge_parallel, принимающая в аргументах void*
void* merge_args_void(void* args_){
  arguments* args = (arguments*) args_;
  merge_parallel(args->array, args->copy_array, args->left_1, args->right_1, args->left_2, args->right_2, args->left_3, args->m);
  return NULL;
}


//функция, сливающая массивы array[left1:right1] и array[left2:right2] в copy_array[left_3:...]
//left_3 - текущий индекс в массиве copy_array, начиная с которого будем записывать полученные слитием данные в copy_array
void merge_parallel(int* array, int* copy_array, int left_1, int right_1, int left_2, int right_2, int left_3, int m){
    int n_1 = right_1 - left_1 + 1;
    int n_2 = right_2 - left_2 + 1;
    //делаем так, чтоб первый подмассив был больше по размеру
    if (n_1 < n_2) {
        swap(&left_1, &left_2);
        swap(&right_1, &right_2);
        swap(&n_1, &n_2);
    }

    if (n_1 == 0) {//Следовательно n_2 = 0
        return;
    }

    //когда подмассивы стали слишком малы быстрей слить их однопоточно, чем плодить много подзадач слития
    if (n_1 <= m) {
      merge(array, copy_array, left_1, right_1, left_2, right_2, left_3);
      return;
    }

    int middle_1 = (right_1 + left_1) / 2;
    int x = array[middle_1];
    //бинпоиском ищем во втором подмассиве первый элемент, который >= x
    int middle_2 = binSearch(array, x, left_2, right_2 + 1);

    //middle_1 - left_1 первых элементов в первом подмассиве <= x
    //middle_2 - left_2 первых элементов в втором подмассиве < x
    int middle_3 = left_3 + (middle_1 - left_1 + middle_2 - left_2);
    copy_array[middle_3] = x;

    pthread_t thread;
    //отвечает за то, был ли создан поток
    int flag = 0;
    arguments* args = (arguments *) malloc(sizeof(arguments));

    pthread_mutex_lock(&mutex);
    //если можем создать поток, то создаем
    if (count_threads > 0){
      count_threads--;
      pthread_mutex_unlock(&mutex);
      flag++;

      //переводим аргументы к типу void*
      args->array = array;
      args->copy_array = copy_array;
      args->left_1 = left_1;
      args->right_1 = middle_1 - 1;
      args->left_2 = left_2;
      args->right_2 = middle_2 - 1;
      args->left_3 = left_3;
      args->m = m;

      //создаем поток и вызываем вспомогательную функцию
      pthread_create(&thread, NULL, merge_args_void, (void*) args);
    }

    //иначе выпоняет текущий поток
    else{
      pthread_mutex_unlock(&mutex);
      merge_parallel(array, copy_array, left_1, middle_1 - 1, left_2, middle_2 - 1, left_3, m);
    }

    //выполняет текущий поток
    merge_parallel(array, copy_array, middle_1 + 1, right_1, middle_2, right_2, middle_3 + 1, m);

    if(flag){
      pthread_join(thread, NULL);
      pthread_mutex_lock(&mutex);
      count_threads++;
      pthread_mutex_unlock(&mutex);
    }
    free(args);
}


int comp (const void * a, const void * b){
    return (*(int *)a - *(int *)b);
}

//вспомогательная функция для sort_parallel, принимающая в аргументах void*
void* sort_args_void(void* args_){
  arguments* args = (arguments*) args_;
  sort_parallel(args->array, args->copy_array, args->left_1, args->right_2, args->m);
  return NULL;
}

void sort_parallel(int* array, int* copy_array, int left, int right, int m) {
    int n  = right - left + 1;

    if(n < m) {
        qsort(array + left, n, sizeof(int), comp);
    }
    else {

        int middle = (right + left) / 2;

        //абсолютно так же как и в merge_parallel
        pthread_t thread;
        int flag = 0;
        arguments* args = (arguments *) malloc(sizeof(arguments));
        pthread_mutex_lock(&mutex);
        if (count_threads > 0){
          count_threads--;
          pthread_mutex_unlock(&mutex);
          flag++;

          args->array = array;
          args->copy_array = copy_array;
          args->left_1 = left;
          args->right_1 = 0;
          args->left_2 = 0;
          args->right_2 = middle;
          args->left_3 = 0;
          args->m = m;

          pthread_create(&thread, NULL, sort_args_void, (void*) args);
        }
        else{
          pthread_mutex_unlock(&mutex);
          sort_parallel(array, copy_array, left, middle, m);
        }

        sort_parallel(array, copy_array, middle + 1, right, m);

        if(flag){
          pthread_join(thread, NULL);
          pthread_mutex_lock(&mutex);
          count_threads++;
          pthread_mutex_unlock(&mutex);
        }
        free(args);
        //сливаем эти уже отсортированные 2 подмассива, не зыбывая о том, что после каждого слития
        // нужно обновить массив array, на который будем опираться при последующих слитиях
        merge_parallel(array, copy_array, left, middle, middle + 1, right, left, m);
        memcpy(array + left, copy_array + left, sizeof(int) * n);
    }
}


int main(int argc, char **argv){

  srand(time(NULL));

  int n, m, P;
  struct timeval begin, end;

	FILE *f_stats = fopen("stats.txt", "a");
  FILE *f_data = fopen("data.txt", "a");

  n = atoi(argv[1]);
	m = atoi(argv[2]);
	P = atoi(argv[3]);
  count_threads = P - 1;

  int* array = (int*)malloc(sizeof(int) * n);

  for (int i = 0; i < n; ++i) {
    array[i] = rand() % 100;
  }
  for (int i = 0; i < n; ++i) {
     fprintf(f_data, "%d ", array[i]);
  }


  int* copy_array = (int*)malloc(sizeof(int) * n);
  memcpy(copy_array, array, sizeof(int) * n);

  assert(gettimeofday(&begin, NULL) == 0);

  sort_parallel(array, copy_array, 0, n - 1, m);

  assert(gettimeofday(&end, NULL) == 0);

  //проверка на корректность
  for (int i = 0; i < n-1;i++){
    assert(copy_array[i] <= copy_array[i+1]);
  }
  double program_time = ((end.tv_sec - begin.tv_sec) * 1000000u + end.tv_usec - begin.tv_usec) / 1.e6;

  fprintf(f_data, "\n");
  for (int i = 0; i < n; ++i) {
    fprintf(f_data, "%d ", copy_array[i]);
  }
  fprintf(f_data, "\n\n");

  fprintf(f_stats, "%f %d %d %d", program_time, n, m, P);
  fprintf(f_stats, "\n");
  fclose(f_stats);
  fclose(f_data);
  free(array);
  free(copy_array);
	return 0;
}
