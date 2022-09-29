CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.x : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out

compila: main.x

esegui:
	./main.x

clean:
	rm *.o

