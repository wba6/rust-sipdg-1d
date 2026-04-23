#!/bin/bash
set -e

SIZES=(5000 10000)

echo "Elements,Main_R1,Main_R2"

for E in "${SIZES[@]}"; do
    sed -i '' "s/num_elements: usize = [0-9]*/num_elements: usize = $E/" src/SIPDG/src/main.rs
    cargo build --release --bin sipdg > /dev/null 2>&1
    
    echo "Running $E..."
    R1=$(/usr/bin/time -p ./target/release/sipdg parameters.txt 2>&1 | grep real | awk '{print $2}')
    R2=$(/usr/bin/time -p ./target/release/sipdg parameters.txt 2>&1 | grep real | awk '{print $2}')
    
    echo "$E,$R1,$R2"
done
