
#[cfg(test)]
mod tests {
    use math::gauss_pp::guass_pp;
    use math::matrix::Matrix;

    fn assert_vec_approx_eq(a: &[f64], b: &[f64], tol: f64) {
        assert_eq!(a.len(), b.len(), "length mismatch");
        for i in 0..a.len() {
            let diff = (a[i] - b[i]).abs();
            assert!(
                diff <= tol,
                "index {}: a={} b={} |diff|={} > tol={}",
                i, a[i], b[i], diff, tol
            );
        }
    }

    #[test]
    fn gauss_pp_solves_2x2_with_row_rhs() {
        // A = [2 1; 5 7], b = [11; 13]
        // Solution: x = 64/9, y = -29/9
        let a = Matrix::from_vec(2, 2, vec![2.0, 1.0,
                                            5.0, 7.0]);
        let f = Matrix::from_vec(1, 2, vec![11.0, 13.0]); // 1×n

        let x = guass_pp(a, f);

        assert_vec_approx_eq(&x, &[64.0 / 9.0, -29.0 / 9.0], 1e-12);
    }

    #[test]
    fn gauss_pp_solves_2x2_with_column_rhs() {
        // Same system as above, but b is n×1
        let a = Matrix::from_vec(2, 2, vec![2.0, 1.0,
                                            5.0, 7.0]);
        let f = Matrix::from_vec(2, 1, vec![11.0, 13.0]); // n×1

        let x = guass_pp(a, f);

        assert_vec_approx_eq(&x, &[64.0 / 9.0, -29.0 / 9.0], 1e-12);
    }

    #[test]
    fn gauss_pp_requires_pivoting_and_gets_correct_answer() {
        // This system will be numerically nasty without pivoting because A(0,0) is tiny.
        // A = [1e-12  1;
        //      1      1], b = [1; 2]
        // Exact-ish solution: x = 1/(1-1e-12), y = 1 - 1e-12*x
        let a = Matrix::from_vec(2, 2, vec![1e-12, 1.0,
                                            1.0,   1.0]);
        let f = Matrix::from_vec(2, 1, vec![1.0, 2.0]);

        let x = guass_pp(a, f);

        let x_expected = 1.0 / (1.0 - 1e-12);
        let y_expected = 1.0 - 1e-12 * x_expected;

        assert_vec_approx_eq(&x, &[x_expected, y_expected], 1e-9);
    }

    #[test]
    fn gauss_pp_solves_3x3_known_solution() {
        // Choose x_true and build b = A*x_true
        let a = Matrix::from_vec(3, 3, vec![
            3.0,  2.0, -1.0,
            2.0, -2.0,  4.0,
           -1.0,  0.5, -1.0
        ]);

        let x_true = vec![1.0, -2.0, -2.0];

        // b = A*x_true
        let mut b = vec![0.0; 3];
        for i in 0..3 {
            let mut s = 0.0;
            for j in 0..3 {
                s += a[(i, j)] * x_true[j];
            }
            b[i] = s;
        }

        let f = Matrix::from_vec(3, 1, b);
        let x = guass_pp(a, f);

        assert_vec_approx_eq(&x, &x_true, 1e-12);
    }

    #[test]
    #[should_panic(expected = "A must be square")]
    fn gauss_pp_panics_if_a_not_square() {
        let a = Matrix::from_vec(2, 3, vec![1.0, 2.0, 3.0,
                                            4.0, 5.0, 6.0]);
        let f = Matrix::from_vec(1, 2, vec![1.0, 1.0]);
        let _ = guass_pp(a, f);
    }

    #[test]
    #[should_panic] // message includes dimensions; keep generic
    fn gauss_pp_panics_if_f_wrong_shape() {
        let a = Matrix::from_vec(2, 2, vec![1.0, 0.0,
                                            0.0, 1.0]);
        // not 1×n and not n×1
        let f = Matrix::from_vec(2, 2, vec![1.0, 2.0,
                                            3.0, 4.0]);
        let _ = guass_pp(a, f);
    }

    #[test]
    #[should_panic(expected = "singular")]
    fn gauss_pp_panics_on_singular_matrix() {
        // Singular: rows are multiples
        let a = Matrix::from_vec(2, 2, vec![1.0, 2.0,
                                            2.0, 4.0]);
        let f = Matrix::from_vec(2, 1, vec![3.0, 6.0]);
        let _ = guass_pp(a, f);
    }

    #[test]
    fn gauss_pp_solves_identity_returns_b() {
        let n = 5;
        let mut data = vec![0.0; n * n];
        for i in 0..n {
            data[i * n + i] = 1.0;
        }
        let a = Matrix::from_vec(n, n, data);

        let b = vec![0.5, -1.0, 2.0, 3.5, -4.0];
        let f = Matrix::from_vec(1, n, b.clone());

        let x = guass_pp(a, f);
        assert_vec_approx_eq(&x, &b, 1e-12);
    }
}

