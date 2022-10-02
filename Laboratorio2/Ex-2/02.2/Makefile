CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.x : main.o random.o RW.o
	$(CC) random.o RW.o main.o -o main.exe
main.o : main.cpp 
	$(CC) -c main.cpp -o main.o $(CFLAGS)
RW.o : RW.cpp RW.h random.h
	$(CC) -c RW.cpp -o RW.o $(CFLAGS)
random.o : random.cpp random.h 
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean:
	rm *.o main.exe 
clean_dat : 
	rm *.dat
compila: 
	main.x

esegui:
	./main.exe


