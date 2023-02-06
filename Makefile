CXX = g++
INCLUDE_DIR = /usr/local/include
CPPFLAGS= -isystem $(INCLUDE_DIR) -std=c++17 -Werror -Wall -Wextra -Wpedantic -Wconversion -DREAL_DIGITS=$(REAL_DIGITS)
LDLIBS = -lgmp -lgmpxx

empty = 
space = $(empty) $(empty)

headers = $(wildcard src/*.h)
modules = $(wildcard src/*.cpp)
tests_with_dir_cpp = $(wildcard test/*.cpp)
tests_cpp = $(tests_with_dir_cpp:test/%=%)
tests_o = $(tests_cpp:%.cpp=build/%.o)
tests = $(tests_cpp:%.cpp=build/%)

runtests = $(subst $(space), && ,$(tests))

REAL_DIGITS ?= 100

.PHONY: all
all: main tests

.PHONY: main
main: build/main

build/main: build/main.o
	$(CXX) $< $(LDLIBS) -o $@

build/main.o: src/main.cpp $(headers) $(modules)
	$(CXX) -c $(CPPFLAGS) $< -o $@

.PHONY: runtests
runtests: $(tests)
	$(runtests)

.PHONY: tests
tests: $(tests)

$(tests): %: %.o
	$(CXX) $< $(LDLIBS) -o $@

$(tests_o): build/%.o: test/%.cpp $(headers) $(modules) 
	$(CXX) -c $(CPPFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -rf build

$(shell mkdir -p build)