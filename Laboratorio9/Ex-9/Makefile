CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o popolazione.o
	$(CC) random.o popolazione.o main.o -o main.exe
main.o : main.cpp input.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
popolazione.o : popolazione.cpp popolazione.h input.h
	$(CC) -c popolazione.cpp  -o popolazione.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
clean_circ : 
	rm Circolare/*.dat popolaz_iniz.dat
clean_quadr : 
	rm Quadrato/*.dat popolaz_iniz.dat

compila: main.exe

esegui:
	./main.exe


