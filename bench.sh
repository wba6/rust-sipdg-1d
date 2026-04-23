#!/bin/bash
set -e

# Update parameters.txt for benchmark
echo "p = 1.0
q = 0.0
f = 1.0
num_elements = 100000
sigma_0 = 20.0" > parameters.txt

echo "Building..."
cargo build --release --bin sipdg > /dev/null 2>&1

echo "Running Benchmarks (100,000 elements)..."
echo "----------------------------------------"

for T in 1 2 4 8; do
    echo -n "Threads $T: "
    RAYON_NUM_THREADS=$T /usr/bin/time -p ./target/release/sipdg parameters.txt 2>&1 | grep real | awk '{print $2}'
done
