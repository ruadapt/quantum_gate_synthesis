CC = g++
CXX = g++
CPPFLAGS= -isystem /usr/local/include -std=c++17 -Werror -Wall -Wextra -Wpedantic -Wconversion
LDLIBS = -lgmp -lgmpxx

empty = 
space = $(empty) $(empty)

headers = $(wildcard *.h)
modules = $(wildcard *.cpp)
tests_cpp = $(wildcard test/*Test.cpp)
tests_o = $(tests_cpp:%.cpp=%.o)
tests = $(tests_cpp:%.cpp=%)

runtests = $(subst $(space), && ,$(tests))

all: main tests

main: main.o
main.o: $(headers) $(modules)

runtests: $(tests)
	$(runtests)

tests: $(tests)

$(tests): %: %.o
$(tests_o): %.o: $(headers) $(modules) 

clean:
	rm -f *.o test/*.o $(tests)