#!/bin/bash
set -e

# Problem parameters for main (using the same as main.rs defaults where possible)
P=1.0
Q=0.0
F=1.0

SIZES=(50 100 250 500 750 1000)

echo "Elements,Main_Time(s)"

for E in "${SIZES[@]}"; do
    # Update main.rs with current size
    sed -i '' "s/num_elements: usize = [0-9]*/num_elements: usize = $E/" src/SIPDG/src/main.rs
    
    cargo build --release --bin SIPDG > /dev/null 2>&1
    
    # Run and time
    TIME=$(/usr/bin/time -p ./target/release/SIPDG parameters.txt 2>&1 | grep real | awk '{print $2}')
    
    echo "$E,$TIME"
done
