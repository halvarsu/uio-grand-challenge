all:
	g++ -g -std=c++11 main.cpp Block.cpp System.cpp Params.cpp Vector.cpp -o main
paranoid:
	g++ -g -std=c++11 -pedantic -Wall -Wextra -Wcast-align -Wcast-qual \
	-Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self \
	-Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept \
	-Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow \
	-Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 \
	-Wswitch-default -Wundef -Werror -Wno-unused main.cpp Block.cpp System.cpp \
	Params.cpp Vector.cpp -o main
fast:
	g++ -std=c++14 -Ofast -march=native -flto -fwhole-program \
	main.cpp Block.cpp Params.cpp System.cpp Vector.cpp -o f
clean:
	rm *o main
	-flto to fast
	-fwhole-program
