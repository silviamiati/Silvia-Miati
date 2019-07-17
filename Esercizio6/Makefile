CC = g++
CFLAGS = -Wall -O3 --std=c++11

ese06 : ese06.o random.o
	$(CC) random.o ese06.o -o ese06
ese06.o : Monte_Carlo_ISING_1D.cpp
	$(CC) -c Monte_Carlo_ISING_1D.cpp -o ese06.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o ese06
