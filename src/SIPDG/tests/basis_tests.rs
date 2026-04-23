// src/SIPDG/tests/basis_tests.rs
// Undergraduate-style tests for Finite Element basis functions.
// We verify properties like Kronecker delta and Partition of Unity.

#[cfg(test)]
mod tests {
    use sipdg::pde::{Element, BasisOrder};

    /// Tolerance for floating point comparisons.
    const TOL: f64 = 1e-12;

    #[test]
    fn test_linear_kronecker_delta() {
        // Linear basis on [-1, 1] (standard reference element)
        // phi_0(-1) = 1, phi_0(1) = 0
        // phi_1(-1) = 0, phi_1(1) = 1
        let elem = Element::new(vec![0, 1], -1.0, 1.0, BasisOrder::Linear);

        assert!((elem.phi(0, -1.0) - 1.0).abs() < TOL);
        assert!((elem.phi(0,  1.0) - 0.0).abs() < TOL);
        assert!((elem.phi(1, -1.0) - 0.0).abs() < TOL);
        assert!((elem.phi(1,  1.0) - 1.0).abs() < TOL);
    }

    #[test]
    fn test_linear_partition_of_unity() {
        // The sum of all basis functions should be exactly 1 everywhere.
        let elem = Element::new(vec![0, 1], -1.0, 1.0, BasisOrder::Linear);
        
        let test_points = vec![-1.0, -0.5, 0.0, 0.5, 1.0];
        for xi in test_points {
            let sum = elem.phi(0, xi) + elem.phi(1, xi);
            assert!((sum - 1.0).abs() < TOL, "Partition of unity failed at xi={}", xi);
        }
    }

    #[test]
    fn test_quadratic_kronecker_delta() {
        // Quadratic basis on [-1, 1]
        // Nodes are at xi = -1, 0, 1
        let elem = Element::new(vec![0, 1, 2], -1.0, 1.0, BasisOrder::Quadratic);

        // phi_0
        assert!((elem.phi(0, -1.0) - 1.0).abs() < TOL);
        assert!((elem.phi(0,  0.0) - 0.0).abs() < TOL);
        assert!((elem.phi(0,  1.0) - 0.0).abs() < TOL);

        // phi_1
        assert!((elem.phi(1, -1.0) - 0.0).abs() < TOL);
        assert!((elem.phi(1,  0.0) - 1.0).abs() < TOL);
        assert!((elem.phi(1,  1.0) - 0.0).abs() < TOL);

        // phi_2
        assert!((elem.phi(2, -1.0) - 0.0).abs() < TOL);
        assert!((elem.phi(2,  0.0) - 0.0).abs() < TOL);
        assert!((elem.phi(2,  1.0) - 1.0).abs() < TOL);
    }

    #[test]
    fn test_quadratic_partition_of_unity() {
        let elem = Element::new(vec![0, 1, 2], -1.0, 1.0, BasisOrder::Quadratic);
        
        for i in 0..11 {
            let xi = -1.0 + 0.2 * (i as f64);
            let sum = elem.phi(0, xi) + elem.phi(1, xi) + elem.phi(2, xi);
            assert!((sum - 1.0).abs() < TOL, "Partition of unity failed at xi={}", xi);
        }
    }

    #[test]
    fn test_linear_derivatives() {
        // For phi_0 = 0.5 * (1 - xi), dphi/dxi = -0.5
        // If element is [0, 2], h=2, J=1, then dphi/dx = dphi/dxi * (1/J) = -0.5
        let elem = Element::new(vec![0, 1], 0.0, 2.0, BasisOrder::Linear);
        
        assert!((elem.dphi_dx(0, 0.0) - (-0.5)).abs() < TOL);
        assert!((elem.dphi_dx(1, 0.0) - ( 0.5)).abs() < TOL);
    }

    #[test]
    fn test_quadratic_derivatives() {
        // phi_1 = 1 - xi^2, dphi/dxi = -2*xi
        // If element is [0, 2], J=1. At xi=0.5, dphi/dx = -1.0
        let elem = Element::new(vec![0, 1, 2], 0.0, 2.0, BasisOrder::Quadratic);
        
        assert!((elem.dphi_dx(1, 0.5) - (-1.0)).abs() < TOL);
        assert!((elem.dphi_dx(1, 0.0) - ( 0.0)).abs() < TOL);
        assert!((elem.dphi_dx(1, -0.5) - ( 1.0)).abs() < TOL);
    }
}
