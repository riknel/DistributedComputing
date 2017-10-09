#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>


void swap(int* a, int* b) {
  int temp = *a;
  *a = *b;
  *b = temp;
}

//функция, находящая индекс первого элемента в отсортированном array[left:right - 1], который >= x
//если такого нет, то возращаем следующий элемент за подмассивом, то есть right
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


//функция, сливающая массивы array[left1:right1] и array[left2:right2] в copy_array
//left_3 - текущий индекс в массиве copy_array, начиная с которого будем записывать полученные данные в copy_array
void merge_parallel(int* array, int* copy_array, int left_1, int right_1, int left_2, int right_2, int left_3){

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

    int middle_1 = (right_1 + left_1) / 2;
    int x = array[middle_1];
    //бинпоиском ищем во втором подмассиве первый элемент, который >= x
    int middle_2 = binSearch(array, x, left_2, right_2 + 1);

    //middle_1 - left_1 первых элементов в первом подмассиве <= x
    //middle_2 - left_2 первых элементов в втором подмассиве < x
    int middle_3 = left_3 + (middle_1 - left_1 + middle_2 - left_2);
    copy_array[middle_3] = x;

#pragma omp parallel
    {
        //single так как каждую задачу(слитие) нужно положить в пул лишь один раз
#pragma omp single nowait
        {
            //кладем в пул 2 задачи слития новых разрезанных подмассивов
#pragma omp task
            {
                merge_parallel(array, copy_array, left_1, middle_1 - 1, left_2, middle_2 - 1, left_3);
            }
#pragma omp task
            {
                merge_parallel(array, copy_array, middle_1 + 1, right_1, middle_2, right_2, middle_3 + 1);
            }
        }
    }
}


int comp (const void * a, const void * b){
    return (*(int *)a - *(int *)b);
}


int main(int argc, char **argv){

	srand(time(NULL));
  int n, m, P;
  clock_t begin;
  clock_t end;

	FILE *f_stats = fopen("stats.txt", "a");
  FILE *f_data = fopen("data.txt", "a");

  n = atoi(argv[1]);
	m = atoi(argv[2]);
	P = atoi(argv[3]);

  int* array = (int*)malloc(sizeof(int) * n);

  for (int i = 0; i < n; ++i) {
    array[i] = rand() % 50;
  }
  for (int i = 0; i < n; ++i) {
     fprintf(f_data, "%d ", array[i]);
  }
  //при таком кол-ве потоков мы более-менее равномерно разделим массив на потоки и при этом учтем,
  //что каждый массив должен быть не длинней m
  int count_arr = n / m + 1;

  omp_set_num_threads(P);
  begin = clock();

  //каждый подмассив сортируем параллельно встроенным qsort
  #pragma omp parallel for
  for (int i = 0; i < count_arr; ++i) {
    //все подмассивы имеют размер m, кроме последнего
    if(i == count_arr - 1) {
      qsort(array + (i * m), n % m, sizeof(int), comp);
    }
    else {
      qsort(array + i * m, m, sizeof(int), comp);
    }
  }

  int* copy_array = (int*)malloc(sizeof(int) * n);
  memcpy(copy_array, array, sizeof(int) * n);

  //сливаем отсортированные массивы следующим образом:
  //сначала первый со вторым, потом полученный с третьим, полученный с четвертым и так далее
  for(int i = 1; i < count_arr; ++i) {
    if(i == count_arr - 1) {
        merge_parallel(array, copy_array, 0, i*m - 1, i*m, n - 1, 0);
    }
    else {
      merge_parallel(array, copy_array, 0, i*m - 1, i*m, (i + 1)* m - 1, 0);
    }
    memcpy(array, copy_array, sizeof(int) * n);
  }

  end = clock();
  double program_time = (double)(end - begin) / CLOCKS_PER_SEC;

  fprintf(f_data, "\n");
  for (int i = 0; i < n; ++i) {
    fprintf(f_data, "%d ", copy_array[i]);
  }
  fprintf(f_data, "\n\n");

  fprintf(f_stats, "%f %d %d %d", program_time, n, m, P);
  fprintf(f_stats, "\n");
  fclose(f_stats);
  fclose(f_data);
	return 0;
}
