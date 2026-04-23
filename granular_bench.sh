#!/bin/bash
set -e

# Problem parameters
P=1.0
Q=0.0
F=1.0
SIGMA=20.0

# Sizes to test
SIZES=(100 250 500 750 1000 1500 2000 3000 5000 10000 20000)

echo "Building..."
cargo build --release --bin sipdg > /dev/null 2>&1

echo "Elements,Serial(s),Parallel(8_threads)"

for E in "${SIZES[@]}"; do
    # Update parameters.txt
    echo "p = $P
q = $Q
f = $F
num_elements = $E
sigma_0 = $SIGMA" > parameters.txt
    
    # Serial run (1 thread)
    SERIAL_TIME=$(RAYON_NUM_THREADS=1 /usr/bin/time -p ./target/release/sipdg parameters.txt 2>&1 | grep real | awk '{print $2}')
    
    # Parallel run (8 threads)
    PARALLEL_TIME=$(RAYON_NUM_THREADS=8 /usr/bin/time -p ./target/release/sipdg parameters.txt 2>&1 | grep real | awk '{print $2}')
    
    echo "$E,$SERIAL_TIME,$PARALLEL_TIME"
done
