MPI_COMPILE = mpic++
COMPILE = g++
MAKE_FLAGS = -Wall -Wextra -pedantic -std=c++11 -O3 -Wshadow -Wfloat-equal -Wconversion -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -fsanitize=address -fsanitize=undefined -fno-sanitize-recover
.PHONY: all clean
all: seq_main par_main mpi_main
clean:
	rm -rf seq_main par_main mpi_main *.out *.data *.txt
seq_main: seq_main.cpp
	$(COMPILE) ${MAKE_FLAGS} -fopenmp -o seq_main seq_main.cpp 
par_main: par_main.cpp
	$(COMPILE) ${MAKE_FLAGS} -fopenmp -o par_main par_main.cpp
mpi_main: mpi_main.cpp
	$(MPI_COMPILE) ${MAKE_FLAGS} -fopenmp -o mpi_main mpi_main.cpp

