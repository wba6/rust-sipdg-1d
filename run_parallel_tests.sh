#!/bin/bash
set -e

OUTPUT="parallel_tests.csv"
echo "Version,Elements,Assembly_Time(ms),Solve_Time(ms)" > $OUTPUT

SIZES=(5000 10000 20000 40000)
LARGE_SIZES=(100000)

CURRENT_BRANCH="34-add-more-unit-tests"

run_tests() {
    local version=$1
    local sizes=("${@:2}")
    
    echo "Running tests for $version..."
    
    # Build
    cargo build --release --bin sipdg
    
    for E in "${sizes[@]}"; do
        echo "  Size: $E"
        echo "p = 1.0
q = 1.0
f = 1.0
num_elements = $E
sigma_0 = 20.0" > parameters.txt
        
        # Warmup
        ./target/release/sipdg parameters.txt > /dev/null 2>&1
        
        # Run and extract CSV_DATA
        ./target/release/sipdg parameters.txt | grep "CSV_DATA:" | sed 's/CSV_DATA://' >> $OUTPUT
    done
}

# 1. Run Parallel (Current Branch)
git checkout $CURRENT_BRANCH > /dev/null 2>&1
run_tests "Parallel" "${SIZES[@]}" "${LARGE_SIZES[@]}"

# 2. Run Serial (Old Branch)
git checkout old-serial > /dev/null 2>&1
run_tests "Serial" "${SIZES[@]}"

# Cleanup: Switch back to original branch
git checkout $CURRENT_BRANCH > /dev/null 2>&1
echo "Benchmarks completed. Results saved to $OUTPUT"
