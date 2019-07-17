CC = mpicxx
CFLAG = -g --std=c++11

ese10 : ese10_2.o random.o
	$(CC) random.o ese10_2.o -o ese10
ese10_2.o : ese10_2.cpp
	$(CC) -c ese10_2.cpp -o ese10_2.o $(CFLAG)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAG)
clean :
	rm *.o ese10

