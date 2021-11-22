CFLAG = -O3 -std=c++11 -mavx2 -fopenmp -Wall 
OBJECTS = main.o SDegreeParallel.o set_intersection.o

SDegreeParallel : main.o SDegreeParallel.o set_intersection.o
	g++ $(CFLAG) -o $@ $^
 
main.o : main.cpp SDegreeParallel.h
	g++ $(CFLAG) -c -o $@ main.cpp 
	
SDegreeParallel.o : SDegreeParallel.cpp SDegreeParallel.h
	g++ $(CFLAG) -c -o $@ SDegreeParallel.cpp 

set_intersection.o : set_intersection.cpp set_intersection.h
	g++ $(CFLAG) -c -o $@ set_intersection.cpp
 
.PHONY : clean
clean:
	rm -rf *.o
	


