CFLAG = -O3 -std=c++11 -fopenmp -Wall 
OBJECTS = BitColNodeParallel.o

BitColNodeParallel : BitColNodeParallel.o
	g++ $(CFLAG) -o $@ $^
 
BitColNodeParallel.o : BitColNodeParallel.cpp
	g++ $(CFLAG) -c -o $@ BitColNodeParallel.cpp 
 
.PHONY : clean
clean:
	rm -rf *.o
	


