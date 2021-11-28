proj3: flow.o utils.o
	nvcc  ./src/proj3.cu csv.o flow.o utils.o -o proj3 -I ./include

flow.o: csv.o
	g++ -c ./src/flow.cpp -I ./include

utils.o:
	g++ -c ./src/utils.cpp -I ./include

csv.o:
	g++ -c ./src/csv.cpp -I ./include

clean:
	rm -f proj3 csv.o utils.o flow.o