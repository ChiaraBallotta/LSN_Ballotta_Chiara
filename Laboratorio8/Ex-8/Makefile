CC = g++
CFLAGS = -Wall -O3 --std=c++11

Ex-8.exe : Ex-8.o random.o
	$(CC) random.o Ex-8.o -o Ex-8.exe
Ex-8.o : Ex-8.cpp
	$(CC) -c Ex-8.cpp -o Ex-8.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Ex-8.exe
clean_all : 
	rm 8.1/*.out 8.2/*.out
clean_1 :
	rm 8.1/*.out
clean_2 : 
	rm 8.2/*.out
esegui :
	./Ex-8.exe
