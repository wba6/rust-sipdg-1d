// src/SIPDG/tests/mesh_tests.rs
// We verify that the mesh correctly places nodes and elements.

#[cfg(test)]
mod tests {
    use sipdg::{generate_mesh, pde::BasisOrder};

    #[test]
    fn test_mesh_linear_counts() {
        let n_elem = 5;
        let order = BasisOrder::Linear;
        let (h_elem, x_dof) = generate_mesh(0.0, 10.0, n_elem, order);

        // For linear (P1), each element has 2 nodes.
        // In DG, nodes are not shared, so we expect n_elem * 2 DoFs.
        assert_eq!(h_elem.len(), n_elem);
        assert_eq!(x_dof.len(), n_elem * 2);

        // Each element should have length 10.0 / 5 = 2.0
        for &h in &h_elem {
            assert!((h - 2.0).abs() < 1e-12);
        }
    }

    #[test]
    fn test_mesh_quadratic_counts() {
        let n_elem = 3;
        let order = BasisOrder::Quadratic;
        let (h_elem, x_dof) = generate_mesh(0.0, 1.0, n_elem, order);

        // For quadratic (P2), each element has 3 nodes.
        assert_eq!(h_elem.len(), n_elem);
        assert_eq!(x_dof.len(), n_elem * 3);
    }

    #[test]
    fn test_mesh_linear_positions() {
        // Domain [0, 1] with 2 elements.
        // E0: [0.0, 0.5], E1: [0.5, 1.0]
        // Nodes should be [0.0, 0.5, 0.5, 1.0] (discontinuous at 0.5)
        let n_elem = 2;
        let order = BasisOrder::Linear;
        let (_, x_dof) = generate_mesh(0.0, 1.0, n_elem, order);

        let expected = vec![0.0, 0.5, 0.5, 1.0];
        for i in 0..x_dof.len() {
            assert!((x_dof[i] - expected[i]).abs() < 1e-12, "Mismatch at index {}", i);
        }
    }

    #[test]
    fn test_mesh_quadratic_positions() {
        // Domain [0, 1] with 1 element.
        // Nodes: Left end, Midpoint, Right end.
        // Expected: [0.0, 0.5, 1.0]
        let n_elem = 1;
        let order = BasisOrder::Quadratic;
        let (_, x_dof) = generate_mesh(0.0, 1.0, n_elem, order);

        let expected = vec![0.0, 0.5, 1.0];
        for i in 0..x_dof.len() {
            assert!((x_dof[i] - expected[i]).abs() < 1e-12, "Mismatch at index {}", i);
        }
    }

    #[test]
    fn test_mesh_negative_domain() {
        let (h_elem, x_dof) = generate_mesh(-5.0, 5.0, 1, BasisOrder::Linear);
        assert_eq!(h_elem[0], 10.0);
        assert_eq!(x_dof[0], -5.0);
        assert_eq!(x_dof[1], 5.0);
    }
}
