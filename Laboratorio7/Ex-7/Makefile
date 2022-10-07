CC = g++
CFLAGS = -Wall -O3 --std=c++11

NVE_NVT.exe : MD_MC.o random.o
	$(CC) random.o MD_MC.o -o NVE_NVT.exe
MD_MC.o : MD_MC.cpp
	$(CC) -c MD_MC.cpp -o MD_MC.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o NVE_NVT.exe
esegui : 
	./NVE_NVT.exe
clean_dat_sol_NVE :
	rm *NVE/output/solid/*.dat
clean_dat_liq_NVE :
	rm *NVE/output/liquid/*.dat
clean_dat_gas_NVE :
	rm *NVE/output/gas/*.dat
clean_all_NVE :
	rm *NVE/output/solid/*.dat *NVE/output/liquid/*.dat *NVE/output/gas/*.dat ./NVE_NVT.exe seed.out
clean_dat_sol_NVT :
	rm *NVT/output/solid/*.dat
clean_dat_liq_NVT :
	rm *NVT/output/liquid/*.dat
clean_dat_gas_NVT :
	rm *NVT/output/gas/*.dat
clean_all_NVT :
	rm *NVT/output/solid/*.dat *NVE/output/liquid/*.dat *NVE/output/gas/*.dat ./NVE_NVT.exe seed.out