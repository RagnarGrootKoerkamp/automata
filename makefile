all: results.txt verify_results

main: main.cpp
	g++ -std=c++2a -O2 -Wall -Wfatal-errors -Wall -Wno-sign-compare main.cpp -o main

results.txt: main
	./main | tee results.txt

verify: verify_finite_order.sage
	sage verify_finite_order.sage < results.txt
