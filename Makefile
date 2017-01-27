
all:
	gcc -c test_c_function.c weak_strong_4d_c.c 
	gcc -c -std=c99 Faddeeva.c 
	gcc -std=c99 -o test_c_function test_c_function.o weak_strong_4d_c.o Faddeeva.o -lm
