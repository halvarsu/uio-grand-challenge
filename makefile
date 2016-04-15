all:
	g++ -g main.cpp Blocks.cpp -o main

main:
	g++ -g main.cpp -o main

blocks:
	g++ -g Blocks.cpp -o Blocks.o
clean:
	rm *o main
