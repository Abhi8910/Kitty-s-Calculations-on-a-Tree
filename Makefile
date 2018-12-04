all: main.cpp
	g++ --std=c++1z main.cpp -o main

test0: all
	cat test0.txt | ./main
