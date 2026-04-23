#!/bin/bash
set -e

# Problem parameters for main
P=1.0
Q=0.0
F=1.0

# Up to 3000 as requested/compared
SIZES=(50 100 250 500 750 1000 3000)

echo "Elements,Main_R1,Main_R2"

for E in "${SIZES[@]}"; do
    # Update main.rs
    sed -i '' "s/num_elements: usize = [0-9]*/num_elements: usize = $E/" src/SIPDG/src/main.rs
    
    cargo build --release --bin SIPDG > /dev/null 2>&1
    
    # Run 1
    R1=$(/usr/bin/time -p ./target/release/SIPDG parameters.txt 2>&1 | grep real | awk '{print $2}')
    # Run 2
    R2=$(/usr/bin/time -p ./target/release/SIPDG parameters.txt 2>&1 | grep real | awk '{print $2}')
    
    echo "$E,$R1,$R2"
done
