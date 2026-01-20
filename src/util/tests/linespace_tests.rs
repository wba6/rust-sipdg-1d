#[cfg(test)]
mod tests {
    use util::linespace::linespace;

    // Helper to compare vectors with a small tolerance for floating point errors
    fn assert_vec_approx_eq(a: &Vec<f64>, b: &Vec<f64>) {
        let epsilon = 1e-10;
        assert_eq!(a.len(), b.len(), "Vector lengths do not match. Got {}, expected {}", a.len(), b.len());
        
        for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
            let diff = (x - y).abs();
            assert!(
                diff < epsilon,
                "Mismatch at index {}: {} != {} (diff: {})", i, x, y, diff
            );
        }
    }

    #[test]
    fn test_standard_ascending() {
        // 0 to 10 with 5 points -> Intervals of 2.5
        // Expected: [0.0, 2.5, 5.0, 7.5, 10.0]
        let result = linespace(0.0, 10.0, 5);
        let expected = vec![0.0, 2.5, 5.0, 7.5, 10.0];
        assert_vec_approx_eq(&result, &expected);
    }

    #[test]
    fn test_standard_descending() {
        // 10 to 0 with 5 points -> Intervals of -2.5
        // Expected: [10.0, 7.5, 5.0, 2.5, 0.0]
        let result = linespace(10.0, 0.0, 5);
        let expected = vec![10.0, 7.5, 5.0, 2.5, 0.0];
        assert_vec_approx_eq(&result, &expected);
    }

    #[test]
    fn test_negative_range() {
        // -2 to 2 with 5 points -> Intervals of 1.0
        // Expected: [-2.0, -1.0, 0.0, 1.0, 2.0]
        let result = linespace(-2.0, 2.0, 5);
        let expected = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
        assert_vec_approx_eq(&result, &expected);
    }

    #[test]
    fn test_single_element() {
        // Edge case: Requesting exactly 1 element should return [start]
        let result = linespace(5.5, 10.0, 1);
        let expected = vec![5.5];
        assert_vec_approx_eq(&result, &expected);
    }

    #[test]
    fn test_two_elements() {
        // Edge case: 2 elements is exactly [start, end]
        let result = linespace(1.0, 5.0, 2);
        let expected = vec![1.0, 5.0];
        assert_vec_approx_eq(&result, &expected);
    }

    #[test]
    fn test_start_equals_end() {
        // If start == end, all points should be the same
        let result = linespace(3.0, 3.0, 4);
        let expected = vec![3.0, 3.0, 3.0, 3.0];
        assert_vec_approx_eq(&result, &expected);
    }

    #[test]
    fn test_fractional_steps() {
        // 0 to 1 with 11 points -> 0.1 steps
        let result = linespace(0.0, 1.0, 11);
        let expected = vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
        assert_vec_approx_eq(&result, &expected);
    }

    #[test]
    #[should_panic(expected = "num_elements must be greater than 0")]
    fn test_panic_zero_elements() {
        // Should trigger the assertion
        linespace(0.0, 10.0, 0);
    }
}
