main:main.cpp
	g++ -O3 -std=c++14 -Wall -Wextra -pedantic -Wshadow -Wfloat-equal -Wconversion -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -fsanitize=address -fsanitize=undefined -fno-sanitize-recover main.cpp -o main -fopenmp


main_nonparallel:main_nonparallel.cpp
	g++ -O3 -std=c++14 -Wall -Wextra -pedantic -Wshadow -Wfloat-equal -Wconversion -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -fsanitize=address -fsanitize=undefined -fno-sanitize-recover main_nonparallel.cpp -o main_nonparallel -fopenmp

clean:
	rm -f main main_nonparallel output.txt output_non_parallel.txt output*.txt