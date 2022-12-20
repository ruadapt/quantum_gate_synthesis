CC = g++
CXX = g++
CPPFLAGS= -isystem /usr/local/include -std=c++17 -Werror -Wall -Wextra -Wpedantic -Wconversion
LDLIBS = -lgmp -lgmpxx

headers = utils.h types.h ring.h gridproblems.h matrix.h quadratic.h toReal.h euclideanDomain.h diophantine.h stepComp.h
modules = ring.cpp gridproblems.cpp matrix.cpp quadratic.cpp toReal.cpp euclideanDomain.cpp diophantine.cpp stepComp.cpp
tests = test/ringTest test/gridproblemsTest test/matrixTest test/quadraticTest test/toRealTest test/euclideanDomainTest test/diophantineTest test/stepCompTest

all: main tests

runtests: $(tests)
	test/ringTest && test/gridproblemsTest && test/matrixTest && test/quadraticTest && test/toRealTest && test/euclideanDomainTest && test/diophantineTest && test/stepCompTest

tests: $(tests)

main: main.o
main.o: $(headers) $(modules)

test/ringTest: test/ringTest.o
test/gridproblemsTest: test/gridproblemsTest.o
test/matrixTest: test/matrixTest.o
test/quadraticTest: test/quadraticTest.o
test/toRealTest: test/toRealTest.o
test/euclideanDomainTest: test/euclideanDomainTest.o
test/diophantineTest: test/diophantineTest.o
test/stepCompTest: test/stepCompTest.o

test/ringTest.o: $(headers) $(modules)
test/gridproblemsTest.o: $(headers) $(modules)
test/matrixTest.o: $(headers) $(modules)
test/quadraticTest.o: $(headers) $(modules)
test/toRealTest.o: $(headers) $(modules)
test/euclideanDomainTest.o: $(headers) $(modules)
test/diophantineTest.o: $(headers) $(modules)
test/stepCompTest.o: $(headers) $(modules)

clean:
	rm -f *.o test/*.o main $(tests)