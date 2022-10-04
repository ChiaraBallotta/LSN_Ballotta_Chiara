CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

main : main.o random.o popolazione.o
	$(CC) random.o popolazione.o main.o -o main.exe
main.o : main.cpp input.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
popolazione.o : popolazione.cpp popolazione.h input.h
	$(CC) -c popolazione.cpp  -o popolazione.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe 
clean_1 :
	rm 1/*.dat
clean_2 :
	rm 2/*.dat
clean_3 :
	rm 3/*.dat
clean_all : clean_3 clean_2 clean_1 

esegui: main
	mpiexec -np 1 ./main.exe


