LIBS:=`root-config --libs`
INCS:=`root-config --cflags`
CFLAGS = -Wall -O3 --std=c++11

make: ese02
ese02: ese02.o integra.o random.o
	g++ ese02.o random.o -o ese02 ${INCS} ${LIBS}
ese02.o: ese02_1.cxx random.h 
	g++ -c ese02_1.cxx -o ese02.o ${INCS}
random.o: random.cxx random.h
	g++ -c random.cxx -o random.o

