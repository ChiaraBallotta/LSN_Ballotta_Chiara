#!/bin/bash
for i in 0.5 0.7 0.9 1.1 1.3 1.5 1.7 1.9
do
    sh ./clean.sh        
    ./Monte_Carlo_ISING_1D.exe ${i} 2   # 2 = GIBBS; 1 = METROPOLIS
done
