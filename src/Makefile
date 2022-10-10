CFLAGS=-g -std=c++17 -Wall -Wextra -Werror 
TESTFLAGS=-lgtest
CC=g++
EXEC_NAME=test.out

SOURCES_LIB=s21_matrix_oop.cc
OBJECTS_LIB=$(SOURCES_LIB:.cc=.o)

all: s21_matrix_oop.a test gcov_report 

s21_matrix_oop.a:
	$(CC) $(CFLAGS) $(SOURCES_LIB) -c
	ar -rcs s21_matrix_oop.a $(OBJECTS_LIB)

leaks: s21_matrix_oop.a
	leaks --atExit -- ./$(EXEC_NAME)

clean:
	rm -rf *.out *.o s21_matrix_oop.a *.gcda *.gcno *.info 
	-rm -rf report

test: s21_matrix_oop.a
	$(CC) $(CFLAGS) s21_matrix_test.cc s21_matrix_oop.a -o test.out -lgtest
	./test.out

gcov_report:
	$(CC) s21_matrix_test.cc -c
	$(CC) --coverage  s21_matrix_oop.cc  s21_matrix_test.o -o test.out $(TESTFLAGS)
	./test.out
	lcov -t "test" -o test.info -c -d ./
	genhtml -o report test.info
	open report/index.html

# support target
clang-format:
	@echo "---clang-format---"
	@clang-format -n *.cc *.h
	@echo "------------------"
	@echo

check: clang-format

.PHONY: all clean test gcov_report