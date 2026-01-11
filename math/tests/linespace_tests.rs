#[cfg(test)]
pub mod tests {
    use math::linespace::linespace;

    // Helper for float comparisons
    fn assert_vec_approx_eq(got: &[f64], expected: &[f64], eps: f64) {
        assert_eq!(got.len(), expected.len(), "length mismatch: got {}, expected {}", got.len(), expected.len());
        for (i, (g, e)) in got.iter().zip(expected.iter()).enumerate() {
            assert!(
                (g - e).abs() <= eps,
                "index {i}: got {g}, expected {e}, diff {} > eps {eps}",
                (g - e).abs()
            );
        }
    }

    #[test]
    fn length_is_num_elements_plus_one() {
        let v = linespace(0.0, 10.0, 5);
        assert_eq!(v.len(), 6);
    }

    #[test]
    fn num_elements_one_gives_two_points() {
        // spacing = |end-start|/1 = 10
        // output = [0, 10]
        let v = linespace(0.0, 10.0, 1);
        assert_vec_approx_eq(&v, &[0.0, 10.0], 1e-12);
    }

    #[test]
    fn evenly_spaced_basic_case() {
        // spacing = |10-0|/4 = 2.5
        // output = [0, 2.5, 5.0, 7.5, 10.0]
        let v = linespace(0.0, 10.0, 4);
        assert_vec_approx_eq(&v, &[0.0, 2.5, 5.0, 7.5, 10.0], 1e-12);
    }

    #[test]
    fn start_is_ignored_in_current_implementation() {
        // spacing = |11-5|/3 = 2
        // output uses index*spacing only => [0,2,4,6], NOT [5,7,9,11]
        let v = linespace(5.0, 11.0, 3);
        assert_vec_approx_eq(&v, &[0.0, 2.0, 4.0, 6.0], 1e-12);
    }

    #[test]
    fn reversed_endpoints_still_produce_increasing_values_due_to_abs() {
        // spacing = |0-10|/4 = 2.5
        // output = [0, 2.5, 5.0, 7.5, 10.0]
        let v = linespace(10.0, 0.0, 4);
        assert_vec_approx_eq(&v, &[0.0, 2.5, 5.0, 7.5, 10.0], 1e-12);
    }

    #[test]
    fn start_equal_end_returns_all_zeros() {
        // spacing = 0, so everything is 0
        let v = linespace(3.14, 3.14, 10);
        assert!(v.iter().all(|&x| x == 0.0));
    }

    #[test]
    fn consecutive_differences_are_constant_spacing() {
        let start = 0.0;
        let end = 1.0;
        let n = 10;
        let v = linespace(start, end, n);

        let spacing = (end - start).abs() / n as f64;

        for i in 0..n {
            let diff = v[i + 1] - v[i];
            assert!(
                (diff - spacing).abs() <= 1e-12,
                "diff at i={i} was {diff}, expected {spacing}"
            );
        }
    }

    #[test]
    fn last_element_equals_abs_end_minus_start() {
        let start = -2.0;
        let end = 3.5;
        let n = 7;
        let v = linespace(start, end, n);

        let expected_last = (end - start).abs();
        let got_last = *v.last().unwrap();

        assert!(
            (got_last - expected_last).abs() <= 1e-12,
            "got last {got_last}, expected {expected_last}"
        );
    }

    #[test]
    #[should_panic]
    fn panics_when_num_elements_is_zero() {
        let _ = linespace(0.0, 1.0, 0);
    }

    #[test]
    fn works_with_negative_endpoints() {
        // spacing = |-1 - (-5)|/4 = 1
        // output = [0,1,2,3,4]
        let v = linespace(-5.0, -1.0, 4);
        assert_vec_approx_eq(&v, &[0.0, 1.0, 2.0, 3.0, 4.0], 1e-12);
    }

    #[test]
    fn floating_point_non_terminating_spacing_is_close() {
        // spacing = |1-0|/3 = 0.333...
        let v = linespace(0.0, 1.0, 3);
        // expected: [0, 1/3, 2/3, 1]
        let expected = [0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0];
        assert_vec_approx_eq(&v, &expected, 1e-12);
    }
}
