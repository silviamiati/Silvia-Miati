CC = g++
CFLAGS = -Wall -O3 --std=c++11

ese02 : ese02.o random.o
	$(CC) random.o ese02.o -o ese02
ese02.o : ese02_2.cxx
	$(CC) -c ese02_2.cxx -o ese02.o $(CFLAGS)
random.o : random.cxx random.h
	$(CC) -c random.cxx -o random.o $(CFLAGS)
clean :
	rm *.o ese02

