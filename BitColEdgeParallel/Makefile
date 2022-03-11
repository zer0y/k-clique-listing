CFLAG = -O3 -std=c++11 -fopenmp -Wall 
OBJECTS = BitColEdgeParallel.o

BitColEdgeParallel : BitColEdgeParallel.o
	g++ $(CFLAG) -o $@ $^
 
BitColEdgeParallel.o : BitColEdgeParallel.cpp
	g++ $(CFLAG) -c -o $@ BitColEdgeParallel.cpp 
 
.PHONY : clean
clean:
	rm -rf *.o
	


