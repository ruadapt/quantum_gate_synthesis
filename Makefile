CXX=g++
CPPFLAGS=-Werror -std=c++17

all: main ringTest

main: main.o
	$(CXX) $(CPPFLAGS) -lgmp -lgmpxx -o main main.o

main.o: ring.h ring.cpp

# Use c++17 for tests.
ringTest: test/ringTest.o
	$(CXX) $(CPPFLAGS) -lgmp -lgmpxx -o ringTest test/ringTest.o

# Use c++17 for tests.
test/ringTest.o: ring.h ring.cpp test/ringTest.cpp
	$(CXX) $(CPPFLAGS) -c -o test/ringTest.o test/ringTest.cpp

clean:
	rm -f *.o test/*.o main ringTest