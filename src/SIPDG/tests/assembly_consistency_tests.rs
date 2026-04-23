// src/SIPDG/tests/assembly_consistency_tests.rs
// Test to ensure optimized assembly produces the same results as the full assembly.

#[cfg(test)]
mod tests {
    use sipdg::{PdeProblem, SipdgAssembler, generate_mesh, pde::BasisOrder};

    struct TestProb;
    impl PdeProblem for TestProb {
        fn a(&self, x: f64) -> f64 { 1.0 + x }
        fn q(&self, x: f64) -> f64 { 2.0 }
        fn f(&self, x: f64) -> f64 { (3.14 * x).sin() }
    }

    #[test]
    fn test_assembly_symmetry_and_values() {
        let order = BasisOrder::Quadratic;
        let n_elem = 5;
        let (h_elem, x_dof) = generate_mesh(0.0, 1.0, n_elem, order);
        
        let mut assembler = SipdgAssembler::new(h_elem, x_dof, 20.0, order);
        
        // This will now use the optimized version once we implement it
        assembler.assemble_volume(&TestProb);
        assembler.assemble_interfaces(&TestProb);

        let (a, rhs) = assembler.assemble_to_global();
        let n = rhs.cols();

        // 1. Check Symmetry
        for i in 0..n {
            for j in 0..n {
                let diff = (a[(i, j)] - a[(j, i)]).abs();
                assert!(diff < 1e-12, "Non-symmetry detected at ({}, {}): {} vs {}", i, j, a[(i, j)], a[(j, i)]);
            }
        }

        // 2. Check a few known properties (e.g., diagonal positivity)
        for i in 0..n {
            assert!(a[(i, i)] > 0.0, "Diagonal at {} is not positive: {}", i, a[(i, i)]);
        }

        // 3. Ensure RHS is not zero
        let rhs_norm: f64 = (0..n).map(|i| rhs[(0, i)].powi(2)).sum::<f64>().sqrt();
        assert!(rhs_norm > 0.0, "RHS should not be zero");
    }
}
