CFLAG = -O3 -std=c++11 -fopenmp -Wall 
OBJECTS = BitColParallel.o

BitColParallel : BitColParallel.o
	g++ $(CFLAG) -o $@ $^
 
BitColParallel.o : BitColParallel.cpp
	g++ $(CFLAG) -c -o $@ BitColParallel.cpp 
 
.PHONY : clean
clean:
	rm -rf *.o
	


