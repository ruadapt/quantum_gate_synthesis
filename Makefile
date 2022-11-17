CC = g++
CXX = g++
CPPFLAGS= -isystem /usr/local/include -std=c++17 -Werror -Wall -Wextra -Wpedantic
LDLIBS = -lgmp -lgmpxx

headers = ring.h gridproblems.h matrix.h
modules = ring.cpp gridproblems.cpp matrix.cpp
tests = test/ringTest test/gridproblemsTest test/matrixTest

all: main tests

tests: $(tests)

main: main.o
main.o: $(headers) $(modules)

test/ringTest: test/ringTest.o
test/gridproblemsTest: test/gridproblemsTest.o
test/matrixTest: test/matrixTest.o

test/ringTest.o: $(headers) $(modules)
test/gridproblemsTest.o: $(headers) $(modules)
test/matrixTest.o: $(headers) $(modules)

clean:
	rm -f *.o test/*.o main $(tests)