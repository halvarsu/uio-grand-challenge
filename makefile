all:
	g++ -g -std=c++11 main.cpp Blocks.cpp Params.cpp -o main

main:
	g++ -g main.cpp -o main

blocks:
	g++ -g Blocks.cpp -o Blocks.o

fast:
	g++ -std=c++11 -Ofast -march=native -flto -fwhole-program main.cpp Blocks.cpp Params.cpp -o main
clean:
	rm *o main
