CC = g++
CFLAGS = -Wall -O3 --std=c++11

ese07 : ese07.o random.o
	$(CC) random.o ese07.o -o ese07
ese07.o : Monte_Carlo_NVT.cpp
	$(CC) -c Monte_Carlo_NVT.cpp -o ese07.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o ese07
