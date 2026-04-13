use crate::matrix::Matrix;

pub trait LinearOperator {
    fn columns(&self) -> usize;
    fn multiply_vec(&self, v: &[f64]) -> Vec<f64>;
}

impl LinearOperator for Matrix<f64> {
    fn columns(&self) -> usize {
        self.cols()
    }
    fn multiply_vec(&self, v: &[f64]) -> Vec<f64> {
        self.multiply_vec(v)
    }
}

/// Conjugate Gradient solver for Ax = b.
/// Assumes A is symmetric and positive-definite.
pub fn cg<A: LinearOperator>(a: &A, f: &Matrix<f64>, tol: f64, max_iter: usize) -> Vec<f64> {
    let n = a.columns();

    // Convert F (either 1xn or nx1) into a vector b of length n
    let b: Vec<f64> = if f.rows() == 1 && f.cols() == n {
        (0..n).map(|i| f[(0, i)]).collect()
    } else if f.cols() == 1 && f.rows() == n {
        (0..n).map(|i| f[(i, 0)]).collect()
    } else {
        panic!("F must be 1×n or n×1 (got {}×{})", f.rows(), f.cols());
    };

    // Initial guess x0 = 0
    let mut x = vec![0.0; n];
    
    // Initial residual r0 = b - A*x0 = b
    let mut r = b.clone();
    
    // Initial search direction p0 = r0
    let mut p = r.clone();
    
    // r^T * r
    let mut r_dot_r = dot(&r, &r);
    
    let tol_sq = tol * tol;
    if r_dot_r < tol_sq {
        return x;
    }

    for iter in 0..max_iter {
        // ap = A * p
        let ap = a.multiply_vec(&p);
        
        // p_ap = p^T * A * p
        let p_ap = dot(&p, &ap);
        
        // Check for breakdown or convergence
        if p_ap.abs() < 1e-20 {
            eprintln!("CG breakdown: p_ap too small at iteration {}", iter);
            break;
        }

        // alpha = (r^T * r) / (p^T * A * p)
        let alpha = r_dot_r / p_ap;
        
        // x_{k+1} = x_k + alpha * p
        for i in 0..n {
            x[i] += alpha * p[i];
        }
        
        // r_{k+1} = r_k - alpha * A * p
        for i in 0..n {
            r[i] -= alpha * ap[i];
        }
        
        let r_dot_r_new = dot(&r, &r);
        
        if r_dot_r_new < tol_sq {
            return x;
        }
        
        // beta = (r_{k+1}^T * r_{k+1}) / (r_k^T * r_k)
        let beta = r_dot_r_new / r_dot_r;
        
        // p_{k+1} = r_{k+1} + beta * p_k
        for i in 0..n {
            p[i] = r[i] + beta * p[i];
        }
        
        r_dot_r = r_dot_r_new;
    }

    x
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}
