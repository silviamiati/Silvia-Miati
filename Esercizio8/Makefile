CC = g++
CFLAGS = -Wall -O3 --std=c++11

ese08 : ese08.o random.o
	$(CC) random.o ese08.o -o ese08
ese08.o : ese.08.cpp
	$(CC) -c ese.08.cpp -o ese08.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o ese08
