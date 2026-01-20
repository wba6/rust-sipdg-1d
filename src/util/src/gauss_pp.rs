use crate::matrix::Matrix;



pub fn guass_pp(mut a: Matrix<f64>, f: Matrix<f64>) -> Vec<f64> {
    // basic checks
    let n = a.rows();
    assert_eq!(a.cols(), n, "A must be square (n x n)");

    // Convert F (either 1xn or nx1) into a vector b of length n
    let mut b: Vec<f64> = vec![0.0; n];
    if f.rows() == 1 && f.cols() == n {
        for i in 0..n {
            b[i] = f[(0, i)];
        }
    } else if f.cols() == 1 && f.rows() == n {
        for i in 0..n {
            b[i] = f[(i, 0)];
        }
    } else {
        panic!("F must be 1×n or n×1 (got {}×{})", f.rows(), f.cols());
    }

    //  forward elimination with partial pivoting 
    let eps = 1e-14;

    for k in 0..(n.saturating_sub(1)) {
        // find pivot row p with max absolute value in column k
        let mut p = k;
        let mut max_val = a[(k, k)].abs();
        for i in (k + 1)..n {
            let val = a[(i, k)].abs();
            if val > max_val {
                max_val = val;
                p = i;
            }
        }

        // swap rows in A and b if needed
        if p != k {
            // swap b entries
            b.swap(k, p);

            // swap the full rows of A
            for j in 0..n {
                let tmp = a[(k, j)];
                a[(k, j)] = a[(p, j)];
                a[(p, j)] = tmp;
            }
        }

        // check pivot
        let pivot = a[(k, k)];
        if pivot.abs() < eps {
            panic!("Matrix is singular or nearly singular at pivot k={}", k);
        }

        // eliminate below pivot
        for i in (k + 1)..n {
            let m = a[(i, k)] / pivot;

            // update row i from column k..n-1
            for j in k..n {
                a[(i, j)] -= m * a[(k, j)];
            }

            // update RHS
            b[i] -= m * b[k];

            // optional cleanup
            a[(i, k)] = 0.0;
        }
    }

    // final pivot check
    if n > 0 && a[(n - 1, n - 1)].abs() < eps {
        panic!("Matrix is singular or nearly singular at final pivot");
    }

    // back substitution 
    let mut x = vec![0.0; n];

    for i in (0..n).rev() {
        let mut sum = 0.0;
        for j in (i + 1)..n {
            sum += a[(i, j)] * x[j];
        }
        x[i] = (b[i] - sum) / a[(i, i)];
    }

    x

}
