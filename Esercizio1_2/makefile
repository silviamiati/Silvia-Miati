LIBS:=`root-config --libs`
INCS:=`root-config --cflags`
CFLAGS = -Wall -O3 --std=c++11

make: es01
es01: es01.o random.o funzione.o
	g++ es01.o  random.o funzione.o -o es01 ${INCS} ${LIBS}
es01.o: es01.2.cxx random.h funzione.h
	g++ -c es01.2.cxx -o es01.o ${INCS}
random.o: random.cpp random.h
	g++ -c random.cpp -o random.o
funzione.o: funzione.cxx funzione.h
	g++ -c funzione.cxx -o funzione.o


