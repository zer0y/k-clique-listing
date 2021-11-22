CFLAG = -O3 -std=c++11 -fopenmp -Wall 
OBJECTS = BitCol.o

BitCol : BitCol.o
	g++ $(CFLAG) -o $@ $^
 
BitCol.o : BitCol.cpp
	g++ $(CFLAG) -c -o $@ BitCol.cpp 
 
.PHONY : clean
clean:
	rm -rf *.o
	


