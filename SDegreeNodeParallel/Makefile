CFLAG = -O3 -std=c++11 -mavx2 -fopenmp -Wall 
OBJECTS = main.o SDegreeNodeParallel.o set_intersection.o

SDegreeNodeParallel : main.o SDegreeNodeParallel.o set_intersection.o
	g++ $(CFLAG) -o $@ $^
 
main.o : main.cpp SDegreeNodeParallel.h
	g++ $(CFLAG) -c -o $@ main.cpp 
	
SDegreeNodeParallel.o : SDegreeNodeParallel.cpp SDegreeNodeParallel.h
	g++ $(CFLAG) -c -o $@ SDegreeNodeParallel.cpp 

set_intersection.o : set_intersection.cpp set_intersection.h
	g++ $(CFLAG) -c -o $@ set_intersection.cpp
 
.PHONY : clean
clean:
	rm -rf *.o
	


