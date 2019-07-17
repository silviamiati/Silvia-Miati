LIBS:=`root-config --libs`
INCS:=`root-config --cflags`


make: ese04
ese04: ese04.o random.o 
	g++ ese04.o random.o -o ese04 ${INCS} ${LIBS}
ese04.o: MolDyn_NVE1.cxx random.h
	g++ -std=c++11 -c MolDyn_NVE1.cxx -o ese04.o ${INCS} 
random.o: random.cxx random.h
	g++ -c random.cxx -o random.o

