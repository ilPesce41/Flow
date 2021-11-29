proj3: flow.o utils.o knn.o
	nvcc ./src/flow_k.cu ./src/spatial_distance.cu ./src/proj3.cu csv.o flow.o utils.o knn.o -lineinfo -o proj3 -I ./include

flow.o: csv.o
	g++ -c ./src/flow.cpp -I ./include

utils.o:
	g++ -c ./src/utils.cpp -I ./include

csv.o:
	g++ -c ./src/csv.cpp -I ./include

knn.o:
	g++ -c ./src/knn.cpp -I ./include


clean:
	rm ./*.o
	rm ./proj3