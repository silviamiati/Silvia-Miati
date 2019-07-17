CC = g++
CFLAGS = -g --std=c++11

ese09 : ese09.o random.o
	$(CC) random.o ese09.o -o ese09
ese09.o : ese09.cpp
	$(CC) -c ese09.cpp -o ese09.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o ese09

