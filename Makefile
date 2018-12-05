main: main.cpp
	g++ --std=c++1z main.cpp -o main

all: main



test0: main
	cat test0.txt | ./main

test1: main
	cat test1.txt | ./main
