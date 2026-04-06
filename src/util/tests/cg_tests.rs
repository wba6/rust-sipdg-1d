
#[cfg(test)]
mod tests {
    use util::cg::cg;
    use util::matrix::Matrix;

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
    fn cg_solves_2x2_spd() {
        // A = [4 1; 1 3], b = [1; 2]
        // SPD matrix.
        // Solution: x = [1/11, 7/11] = [0.090909, 0.636363]
        let a = Matrix::from_vec(2, 2, vec![4.0, 1.0,
                                            1.0, 3.0]);
        let f = Matrix::from_vec(2, 1, vec![1.0, 2.0]);

        let x = cg(&a, &f, 1e-12, 100);

        assert_vec_approx_eq(&x, &[1.0/11.0, 7.0/11.0], 1e-10);
    }

    #[test]
    fn cg_solves_3x3_spd() {
        // A = [2 -1 0; -1 2 -1; 0 -1 2], b = [1; 0; 1]
        // SPD (Tridiagonal matrix with 2 on diagonal and -1 on off-diagonals)
        // Solution: x = [1, 1, 1]
        let a = Matrix::from_vec(3, 3, vec![
            2.0, -1.0,  0.0,
           -1.0,  2.0, -1.0,
            0.0, -1.0,  2.0
        ]);
        let f = Matrix::from_vec(3, 1, vec![1.0, 0.0, 1.0]);

        let x = cg(&a, &f, 1e-12, 100);

        assert_vec_approx_eq(&x, &[1.0, 1.0, 1.0], 1e-10);
    }

    #[test]
    fn cg_solves_identity() {
        let n = 5;
        let mut data = vec![0.0; n * n];
        for i in 0..n {
            data[i * n + i] = 1.0;
        }
        let a = Matrix::from_vec(n, n, data);

        let b = vec![0.5, -1.0, 2.0, 3.5, -4.0];
        let f = Matrix::from_vec(1, n, b.clone());

        let x = cg(&a, &f, 1e-12, 100);
        assert_vec_approx_eq(&x, &b, 1e-10);
    }
}
