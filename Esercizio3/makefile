LIBS:=`root-config --libs`
INCS:=`root-config --cflags`
CFLAGS = -Wall -O3 --std=c++11

make: ese03
ese03: ese03.o funzione.o random.o 
	g++ ese03.o random.o funzione.o -o ese03 ${INCS} ${LIBS}
ese03.o: ese03_1.cxx random.h funzione.h 
	g++ -c ese03_1.cxx -o ese03.o ${INCS}
random.o: random.cxx random.h
	g++ -c random.cxx -o random.o
funzione.o: funzione.cxx funzione.h
	g++ -c funzione.cxx -o funzione.o

