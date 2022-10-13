CXX=g++
CPPFLAGS=-Werror -std=c++11

all: main ringTest

main: main.o
	$(CXX) $(CPPFLAGS) -lgmp -lgmpxx -o main main.o

main.o: ring.cpp

ringTest: test/ringTest.o
	$(CXX) $(CPPFLAGS) -lgmp -lgmpxx -o ringTest test/ringTest.o

test/ringTest.o: ring.cpp

clean:
	rm -f *.o test/*.o main ringTest