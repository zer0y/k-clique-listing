CFLAG = -O3 -std=c++11 -Wall 
Preprocess : Pre.o Graph.o
	g++ $(CFLAG) -o $@ $^
Pre.o : Pre.cpp Graph.h
	g++ $(CFLAG) -c -o $@ Pre.cpp 
	
Graph.o : Graph.cpp Graph.h
	g++ $(CFLAG) -c -o $@ Graph.cpp 
   
.PHONY : clean
clean:
	rm -rf *.o
	


