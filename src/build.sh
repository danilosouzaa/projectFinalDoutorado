nvcc --compiler-options "-Wall -O2 -pg -g -std=c99 -lpthread -fopenmp -lgomp -DDEBUG -DGRB -O3"  lp.cpp main.cpp preProcessing.c constraintsManipulation.c coverLifted.c -I/opt/gurobi811/linux64/include/ -L/opt/gurobi811/linux64/lib/ -lgurobi81  -o projectFinal2

