// src/SIPDG/tests/symmetry_tests.rs
// Undergraduate-style tests to verify the properties of the assembled system matrix.
// For the SIPDG method, the matrix should be symmetric if we use the symmetric version of the flux.

#[cfg(test)]
mod tests {
    use sipdg::{PdeProblem, SipdgAssembler, generate_mesh, pde::BasisOrder};

    struct SimpleProb;
    impl PdeProblem for SimpleProb {
        fn a(&self, _x: f64) -> f64 { 1.0 }
        fn q(&self, _x: f64) -> f64 { 1.0 }
        fn f(&self, _x: f64) -> f64 { 1.0 }
    }

    #[test]
    fn test_matrix_symmetry() {
        let order = BasisOrder::Linear;
        let n_elem = 4;
        let (h_elem, x_dof) = generate_mesh(0.0, 1.0, n_elem, order);
        
        // Assemble with a standard penalty
        let mut assembler = SipdgAssembler::new(h_elem, x_dof, 20.0, order);
        assembler.assemble_volume(&SimpleProb);
        assembler.assemble_interfaces(&SimpleProb);

        let (a, rhs) = assembler.assemble_to_global();

        // Check symmetry: A[i, j] == A[j, i]
        let n = rhs.cols();
        let tol = 1e-12;
        for i in 0..n {
            for j in 0..n {
                let diff = (a[(i, j)] - a[(j, i)]).abs();
                assert!(diff < tol, "Matrix is not symmetric at ({}, {}): Aij={}, Aji={}", i, j, a[(i, j)], a[(j, i)]);
            }
        }
    }

    #[test]
    fn test_diagonal_positivity() {
        // For elliptic problems like this, the diagonal entries of the stiffness matrix should be positive.
        let order = BasisOrder::Quadratic;
        let n_elem = 3;
        let (h_elem, x_dof) = generate_mesh(0.0, 1.0, n_elem, order);
        
        let mut assembler = SipdgAssembler::new(h_elem, x_dof, 10.0, order);
        assembler.assemble_volume(&SimpleProb);
        assembler.assemble_interfaces(&SimpleProb);

        let (a, rhs) = assembler.assemble_to_global();

        let n = rhs.cols();
        for i in 0..n {
            assert!(a[(i, i)] > 0.0, "Diagonal entry A[{}, {}] is not positive: {}", i, i, a[(i, i)]);
        }
    }

    #[test]
    fn test_load_vector_constant() {
        // For f(x) = 1.0, the load vector entries on an element of length h should be:
        // Linear: [h/2, h/2]
        // Quadratic: [h/6, 4h/6, h/6]
        let n_elem = 2;
        let h = 0.5;
        
        // Test Linear
        {
            let (h_elem, x_dof) = generate_mesh(0.0, 1.0, n_elem, BasisOrder::Linear);
            let mut assembler = SipdgAssembler::new(h_elem, x_dof, 20.0, BasisOrder::Linear);
            assembler.assemble_volume(&SimpleProb);
            let (_, rhs) = assembler.assemble_to_global();
            
            for i in 0..rhs.cols() {
                assert!((rhs[(0, i)] - h/2.0).abs() < 1e-12, "Linear RHS mismatch at {}: expected {}, got {}", i, h/2.0, rhs[(0, i)]);
            }
        }

        // Test Quadratic
        {
            let (h_elem, x_dof) = generate_mesh(0.0, 1.0, n_elem, BasisOrder::Quadratic);
            let mut assembler = SipdgAssembler::new(h_elem, x_dof, 20.0, BasisOrder::Quadratic);
            assembler.assemble_volume(&SimpleProb);
            let (_, rhs) = assembler.assemble_to_global();
            
            let expected_h6 = h / 6.0;
            let expected_4h6 = 4.0 * h / 6.0;
            
            for e in 0..n_elem {
                assert!((rhs[(0, 3*e)] - expected_h6).abs() < 1e-12);
                assert!((rhs[(0, 3*e + 1)] - expected_4h6).abs() < 1e-12);
                assert!((rhs[(0, 3*e + 2)] - expected_h6).abs() < 1e-12);
            }
        }
    }
}
