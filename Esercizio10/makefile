CC = g++
CFLAGS = -Wall -O3 --std=c++11

ese10 : ese10.o random.o
	$(CC) random.o ese10.o -o ese10
ese10.o : ese10.cpp
	$(CC) -c ese10.cpp -o ese10.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o ese10

