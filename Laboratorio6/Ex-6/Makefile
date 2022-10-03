CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_ISING_1D.exe : Monte_Carlo_ISING_1D.o random.o
	$(CC) random.o Monte_Carlo_ISING_1D.o -o Monte_Carlo_ISING_1D.exe
Monte_Carlo_ISING_1D.o : Monte_Carlo_ISING_1D.cpp
	$(CC) -c Monte_Carlo_ISING_1D.cpp -o Monte_Carlo_ISING_1D.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Monte_Carlo_ISING_1D.exe seed.out
clean_config :
	rm config.final
clean_all :
	rm *.o Monte_Carlo_ISING_1D.exe seed.out config.final  *.0
clean_output :
	rm *.0
clean_dat :
	rm output/Metropolis/*.dat output/Gibbs/*.dat
clean_metr : 
	rm output/Metropolis/*.dat
clean_gibbs :
	rm output/Gibbs/*.dat
esegui :
	./Monte_Carlo_ISING_1D.exe


