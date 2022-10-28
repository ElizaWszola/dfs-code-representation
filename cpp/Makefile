graph: graph_fast.cpp
	g++ -std=c++17 -O3 graph_fast.cpp -o graph

main:   main.cpp graph.cpp graph.hpp common.hpp common.cpp
	g++ -std=c++17 -O3 graph.cpp main.cpp common.cpp -o main

fast:   main.cpp graph_fast.cpp graph_fast.hpp correlation.hpp correlation.cpp common.hpp common.cpp
	g++ -std=c++17 -O3 -DFAST graph_fast.cpp main.cpp correlation.cpp common.cpp -o main
