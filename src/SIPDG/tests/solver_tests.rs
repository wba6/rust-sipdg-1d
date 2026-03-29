// src/SIPDG/tests/solver_tests.rs

use SIPDG::{PdeProblem, SipdgAssembler, DirichletBC, generate_mesh};
use util::gauss_pp::gauss_pp;
use std::f64::consts::PI;

/// A flexible problem struct that lets us define p, q, f, and exact_u using closures.
struct TestProblem<P, Q, F, U>
where
    P: Fn(f64) -> f64,
    Q: Fn(f64) -> f64,
    F: Fn(f64) -> f64,
    U: Fn(f64) -> f64,
{
    p_fn: P,
    q_fn: Q,
    f_fn: F,
    u_exact: U,
}

impl<P, Q, F, U> PdeProblem for TestProblem<P, Q, F, U>
where
    P: Fn(f64) -> f64,
    Q: Fn(f64) -> f64,
    F: Fn(f64) -> f64,
    U: Fn(f64) -> f64,
{
    fn p(&self, x: f64) -> f64 { (self.p_fn)(x) }
    fn q(&self, x: f64) -> f64 { (self.q_fn)(x) }
    fn f(&self, x: f64) -> f64 { (self.f_fn)(x) }
}

/// Helper function to run the full solver pipeline and return the L2 error.
/// Now accepts `penalty` to allow tuning for convergence tests.
fn run_solver_and_compute_error(
    prob: &impl PdeProblem, 
    exact_soln: impl Fn(f64) -> f64, 
    num_elements: usize,
    penalty: f64
) -> f64 {
    // Generate Mesh
    let (h_elem, x_dof) = generate_mesh(0.0, 1.0, num_elements);
    
    // Assemble (using the provided penalty parameter)
    let mut assembler = SipdgAssembler::new(h_elem.clone(), x_dof.clone(), penalty);
    assembler.assemble_volume(prob);
    assembler.assemble_interfaces(prob);

    let (mut a, mut rhs) = assembler.assemble_to_global();

    // Apply BCs (assuming homogeneous Dirichlet u(0)=u(1)=0 for these tests)
    let bc = DirichletBC { value: 0.0 };
    assembler.apply_boundaries(prob, &bc, &bc, &mut a, &mut rhs);

    // Solve
    let u = gauss_pp(a, rhs);

    // Compute L2 Error (Trapezoidal rule)
    let mut l2_err_sq = 0.0;
    let npts = 10;

    for k in 0..num_elements {
        let idx0 = 2 * k;
        let idx1 = 2 * k + 1;
        let x_l = x_dof[idx0];
        let x_r = x_dof[idx1];
        let h = h_elem[k];
        let ul = u[idx0];
        let ur = u[idx1];

        let step = (x_r - x_l) / ((npts - 1) as f64);
        let mut prev_x = x_l;
        
        // Helper to compute error squared at a point
        let get_err2 = |x: f64| {
            let phi1 = (x_r - x) / h;
            let phi2 = (x - x_l) / h;
            let uh = ul * phi1 + ur * phi2;
            let uex = exact_soln(x);
            (uh - uex).powi(2)
        };

        let mut prev_err2 = get_err2(prev_x);

        for i in 1..npts {
            let x = x_l + step * (i as f64);
            let err2 = get_err2(x);
            
            // Trapezoid area
            l2_err_sq += 0.5 * (x - prev_x) * (prev_err2 + err2);
            
            prev_x = x;
            prev_err2 = err2;
        }
    }

    l2_err_sq.sqrt()
}

#[test]
fn test_exact_linear_solution() {
    // Problem: -u'' = 1, u(0)=0, u(1)=0 
    // Exact: u(x) = x * (1.0 - x) / 2.0 (Quadratic)
    // NOTE: Linear elements cannot capture quadratic exactly, but error should be small.
    
    let prob = TestProblem {
        p_fn: |_| 1.0,
        q_fn: |_| 0.0,
        f_fn: |_| 1.0, 
        u_exact: |x| x * (1.0 - x) / 2.0,
    };

    // Use standard penalty of 80.0
    let error = run_solver_and_compute_error(&prob, &prob.u_exact, 10, 80.0);
    
    // Error for N=10 should be around ~1e-3
    assert!(error < 5e-3, "Error was too high for simple Poisson: {}", error);
}

#[test]
fn test_convergence_poisson_sine() {
    // Problem: -u'' = f
    // Target: u(x) = sin(pi * x)
    // f(x) = pi^2 * sin(pi * x)
    
    let prob = TestProblem {
        p_fn: |_| 1.0,
        q_fn: |_| 0.0,
        f_fn: |x| PI.powi(2) * (PI * x).sin(),
        u_exact: |x| (PI * x).sin(),
    };

    let penalty = 10.0;

    let err_coarse = run_solver_and_compute_error(&prob, &prob.u_exact, 10, penalty);
    let err_fine = run_solver_and_compute_error(&prob, &prob.u_exact, 20, penalty);

    println!("Sine Test: Coarse Err = {:.4e}, Fine Err = {:.4e}", err_coarse, err_fine);

    // For linear elements, L2 error is O(h^2). Doubling elements = 1/4 error.
    let rate = err_coarse / err_fine;
    assert!(rate > 3.0, "Convergence rate too slow! Expected ~4.0, got {:.2}", rate);
}

#[test]
fn test_reaction_term() {
    // Problem: -u'' + u = f
    // Target: u(x) = sin(pi * x)
    // f(x) = (pi^2 + 1) * sin(pi * x)
    
    let prob = TestProblem {
        p_fn: |_| 1.0,
        q_fn: |_| 1.0,
        f_fn: |x| (PI.powi(2) + 1.0) * (PI * x).sin(),
        u_exact: |x| (PI * x).sin(),
    };

    let err = run_solver_and_compute_error(&prob, &prob.u_exact, 20, 80.0);
    
    assert!(err < 1e-2, "Reaction-Diffusion error too high: {}", err);
}

#[test]
fn test_variable_coefficient_p() {
    // Problem: -( (1+x) u' )' = f
    // Target: u(x) = x * (1.0 - x)
    // f(x) = 1 + 4x
    
    let prob = TestProblem {
        p_fn: |x| 1.0 + x,
        q_fn: |_| 0.0,
        f_fn: |x| 1.0 + 4.0 * x,
        u_exact: |x| x * (1.0 - x),
    };

    // approximating a quadratic function + variable coefficients.
    let err = run_solver_and_compute_error(&prob, &prob.u_exact, 20, 80.0);

    assert!(err < 5e-3, "Variable coefficient error too high: {}", err);
}
