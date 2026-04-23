#!/bin/bash
set -e

# Problem parameters
P=1.0
Q=0.0
F=1.0
SIGMA=20.0

SIZES=(50 100 250 500 750 1000 2000 3000 5000)

echo "Building..."
cargo build --release --bin sipdg > /dev/null 2>&1

echo "Elements,Serial_1_thread(s),Parallel_8_threads(s)"

for E in "${SIZES[@]}"; do
    # Update parameters.txt
    echo "p = $P
q = $Q
f = $F
num_elements = $E
sigma_0 = $SIGMA" > parameters.txt
    
    # Serial run (1 thread)
    SERIAL_TIME=$(RAYON_NUM_THREADS=1 ./target/release/sipdg parameters.txt 2>&1 > /dev/null | /usr/bin/time -p cat 2>&1 | grep real | awk '{print $2}')
    # Wait, the timing above might not be accurate because of 'cat'. 
    # Let's use the same method as before.
    SERIAL_TIME=$(RAYON_NUM_THREADS=1 /usr/bin/time -p ./target/release/sipdg parameters.txt 2>&1 | grep real | awk '{print $2}')
    PARALLEL_TIME=$(RAYON_NUM_THREADS=8 /usr/bin/time -p ./target/release/sipdg parameters.txt 2>&1 | grep real | awk '{print $2}')
    
    echo "$E,$SERIAL_TIME,$PARALLEL_TIME"
done
