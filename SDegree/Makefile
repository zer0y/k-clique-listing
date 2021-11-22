CFLAG = -O3 -std=c++11 -mavx2 -Wall 
OBJECTS = main.o SDegree.o set_intersection.o

SDegree : main.o SDegree.o set_intersection.o
	g++ $(CFLAG) -o $@ $^
 
main.o : main.cpp SDegree.h
	g++ $(CFLAG) -c -o $@ main.cpp 
	
SDegree.o : SDegree.cpp SDegree.h
	g++ $(CFLAG) -c -o $@ SDegree.cpp 

set_intersection.o : set_intersection.cpp set_intersection.h
	g++ $(CFLAG) -c -o $@ set_intersection.cpp
 
.PHONY : clean
clean:
	rm -rf *.o
	


