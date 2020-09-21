COMPILE = g++
MAKE_FLAGS = -Wall -Wextra -pedantic -std=c++1z -O3 -Wshadow -Wfloat-equal -Wconversion -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -fsanitize=address -fsanitize=undefined -fno-sanitize-recover
.PHONY: all clean
all: main
clean:
	rm -rf main *.out *.data
main: main.cpp
	$(COMPILE) ${MAKE_FLAGS} -fopenmp -o main main.cpp 
