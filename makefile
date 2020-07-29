main:
	g++ -std=c++2a -O2 -Wall -Wfatal-errors -Wall -Wno-sign-compare main.cpp -o main

results.txt: main
	./main | tee results.txt
