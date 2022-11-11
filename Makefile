CXX=g++
CPPFLAGS=-Werror -Wall -Wextra -Wpedantic -std=c++17

all: main ringTest gridProblemsTest

main: main.o
	$(CXX) $(CPPFLAGS) -lgmp -lgmpxx -o main main.o

main.o: ring.h ring.cpp gridproblems.h gridproblems.cpp

gridProblemsTest: test/gridproblemsTest.o
	$(CXX) $(CPPFLAGS) -lgmp -lgmpxx -o gridProblemsTest test/gridProblemsTest.o

test/gridproblemsTest.o: ring.h ring.cpp gridproblems.h gridproblems.cpp

ringTest: test/ringTest.o
	$(CXX) $(CPPFLAGS) -lgmp -lgmpxx -o ringTest test/ringTest.o

test/ringTest.o: ring.h ring.cpp

clean:
	rm -f *.o test/*.o main ringTest