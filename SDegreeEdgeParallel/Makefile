CFLAG = -O3 -std=c++11 -mavx2 -fopenmp -Wall 
OBJECTS = main.o SDegreeEdgeParallel.o set_intersection.o

SDegreeEdgeParallel : main.o SDegreeEdgeParallel.o set_intersection.o
	g++ $(CFLAG) -o $@ $^
 
main.o : main.cpp SDegreeEdgeParallel.h
	g++ $(CFLAG) -c -o $@ main.cpp 
	
SDegreeEdgeParallel.o : SDegreeEdgeParallel.cpp SDegreeEdgeParallel.h
	g++ $(CFLAG) -c -o $@ SDegreeEdgeParallel.cpp 

set_intersection.o : set_intersection.cpp set_intersection.h
	g++ $(CFLAG) -c -o $@ set_intersection.cpp
 
.PHONY : clean
clean:
	rm -rf *.o
	


