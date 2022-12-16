CC = g++
CXX = g++
CPPFLAGS= -isystem /usr/local/include -std=c++17 -Werror -Wall -Wextra -Wpedantic -Wconversion
LDLIBS = -lgmp -lgmpxx

headers = types.h ring.h gridproblems.h matrix.h quadratic.h toReal.h euclideanDomain.h
modules = ring.cpp gridproblems.cpp matrix.cpp quadratic.cpp toReal.cpp euclideanDomain.cpp
tests = test/ringTest test/gridproblemsTest test/matrixTest test/quadraticTest test/toRealTest test/euclideanDomainTest

all: main tests

runtests: $(tests)
	test/ringTest && test/gridproblemsTest && test/matrixTest && test/quadraticTest && test/toRealTest && test/euclideanDomainTest

tests: $(tests)

main: main.o
main.o: $(headers) $(modules)

test/ringTest: test/ringTest.o
test/gridproblemsTest: test/gridproblemsTest.o
test/matrixTest: test/matrixTest.o
test/quadraticTest: test/quadraticTest.o
test/toRealTest: test/toRealTest.o
test/euclideanDomainTest: test/euclideanDomainTest.o

test/ringTest.o: $(headers) $(modules)
test/gridproblemsTest.o: $(headers) $(modules)
test/matrixTest.o: $(headers) $(modules)
test/quadraticTest.o: $(headers) $(modules)
test/toRealTest.o: $(headers) $(modules)
test/euclideanDomainTest.o: $(headers) $(modules)

clean:
	rm -f *.o test/*.o main $(tests)