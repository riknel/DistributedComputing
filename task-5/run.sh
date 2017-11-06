mpicc main.c -o a.out
mpirun -np 4 --oversubscribe ./a.out 10 2 2 100 10 0.25 0.25 0.25 0.25
