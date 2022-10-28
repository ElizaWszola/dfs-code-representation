graph: graph.cpp
	g++ -std=c++11 -O3 graph.cpp -o graph

main:   main.cpp graph.cpp graph.hpp common.hpp common.cpp
	g++ -std=c++11 -O3 graph.cpp main.cpp common.cpp -o main

fast:   main.cpp graph_fast.cpp graph_fast.hpp common.hpp common.cpp
	g++ -std=c++11 -O3 -DFAST graph_fast.cpp main.cpp common.cpp -o main
