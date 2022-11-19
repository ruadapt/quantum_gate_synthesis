CC = g++
CXX = g++
CPPFLAGS= -isystem /usr/local/include -std=c++17 -Werror -Wall -Wextra -Wpedantic
LDLIBS = -lgmp -lgmpxx

headers = ring.h gridproblems.h matrix.h quadratic.h toReal.h
modules = ring.cpp gridproblems.cpp matrix.cpp quadratic.cpp toReal.cpp
tests = test/ringTest test/gridproblemsTest test/matrixTest test/quadraticTest test/toRealTest

all: main tests

tests: $(tests)

main: main.o
main.o: $(headers) $(modules)

test/ringTest: test/ringTest.o
test/gridproblemsTest: test/gridproblemsTest.o
test/matrixTest: test/matrixTest.o
test/quadraticTest: test/quadraticTest.o
test/toRealTest: test/toRealTest.o

test/ringTest.o: $(headers) $(modules)
test/gridproblemsTest.o: $(headers) $(modules)
test/matrixTest.o: $(headers) $(modules)
test/quadraticTest.o: $(headers) $(modules)
test/toRealTest.o: $(headers) $(modules)

clean:
	rm -f *.o test/*.o main $(tests)