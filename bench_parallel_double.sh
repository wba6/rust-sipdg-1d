#!/bin/bash
set -e

P=1.0
Q=0.0
F=1.0
SIGMA=20.0

SIZES=(50 100 250 500 750 1000 2000 3000 5000 10000 20000 100000)

echo "Building..."
cargo build --release --bin sipdg > /dev/null 2>&1

echo "Elements,Serial_R1,Serial_R2,Parallel_R1,Parallel_R2"

for E in "${SIZES[@]}"; do
    echo "p = $P
q = $Q
f = $F
num_elements = $E
sigma_0 = $SIGMA" > parameters.txt
    
    # Serial Run 1
    S1=$(RAYON_NUM_THREADS=1 /usr/bin/time -p ./target/release/sipdg parameters.txt 2>&1 | grep real | awk '{print $2}')
    # Serial Run 2
    S2=$(RAYON_NUM_THREADS=1 /usr/bin/time -p ./target/release/sipdg parameters.txt 2>&1 | grep real | awk '{print $2}')
    
    # Parallel Run 1
    P1=$(RAYON_NUM_THREADS=8 /usr/bin/time -p ./target/release/sipdg parameters.txt 2>&1 | grep real | awk '{print $2}')
    # Parallel Run 2
    P2=$(RAYON_NUM_THREADS=8 /usr/bin/time -p ./target/release/sipdg parameters.txt 2>&1 | grep real | awk '{print $2}')
    
    echo "$E,$S1,$S2,$P1,$P2"
done
