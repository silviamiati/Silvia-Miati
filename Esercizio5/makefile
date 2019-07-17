LIBS:=`root-config --libs`
INCS:=`root-config --cflags`
CFLAGS = -Wall -O3 --std=c++11


make: es05
es05: es05.o random.o 
	g++ es05.o random.o -o es05 ${INCS} ${LIBS}
es05.o: es05.cxx random.h 
	g++ -c es05.cxx -o es05.o ${INCS}
random.o: random.cpp random.h
	g++ -c random.cpp -o random.o

