main:main.cpp
	mpicxx -O3 -std=c++14 main.cpp -o main -fopenmp


main_nonparallel:main_nonparallel.cpp
	mpicxx -O3 -std=c++14  main_nonparallel.cpp -o main_nonparallel -fopenmp

clean:
	rm -f main main_nonparallel output.txt output_non_parallel.txt output*.txt