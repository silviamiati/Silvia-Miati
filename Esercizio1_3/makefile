LIBS:=`root-config --libs`
INCS:=`root-config --cflags`
CFLAGS = -Wall -O3 --std=c++11

make: es01
es01: es01.o random.o 
	g++ es01.o  random.o -o es01 ${INCS} ${LIBS}
es01.o: es01.3.cxx random.h
	g++ -c es01.3.cxx -o es01.o ${INCS}
random.o: random.cpp random.h
	g++ -c random.cpp -o random.o

