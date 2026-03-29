// src/SIPDG/tests/solver_tests.rs

use SIPDG::{PdeProblem, SipdgAssembler, DirichletBC, generate_mesh};
use util::gauss_pp::gauss_pp;
use std::f64::consts::PI;

/// A flexible problem struct that lets us define a, q, f, and exact_p using closures.
struct TestProblem<A, Q, F, P>
where
    A: Fn(f64) -> f64,
    Q: Fn(f64) -> f64,
    F: Fn(f64) -> f64,
    P: Fn(f64) -> f64,
{
    a_fn: A,
    q_fn: Q,
    f_fn: F,
    p_exact: P,
}

impl<A, Q, F, P> PdeProblem for TestProblem<A, Q, F, P>
where
    A: Fn(f64) -> f64,
    Q: Fn(f64) -> f64,
    F: Fn(f64) -> f64,
    P: Fn(f64) -> f64,
{
    fn a(&self, x: f64) -> f64 { (self.a_fn)(x) }
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

    // Apply BCs (assuming homogeneous Dirichlet p(0)=p(1)=0 for these tests)
    let bc = DirichletBC { value: 0.0 };
    assembler.apply_boundaries(prob, &bc, &bc, &mut a, &mut rhs);

    // Solve
    let p = gauss_pp(a, rhs);

    // Compute L2 Error (Trapezoidal rule)
    let mut l2_err_sq = 0.0;
    let npts = 10;

    for k in 0..num_elements {
        let idx0 = 2 * k;
        let idx1 = 2 * k + 1;
        let x_l = x_dof[idx0];
        let x_r = x_dof[idx1];
        let h = h_elem[k];
        let pl = p[idx0];
        let pr = p[idx1];

        let step = (x_r - x_l) / ((npts - 1) as f64);
        let mut prev_x = x_l;
        
        // Helper to compute error squared at a point
        let get_err2 = |x: f64| {
            let phi1 = (x_r - x) / h;
            let phi2 = (x - x_l) / h;
            let ph = pl * phi1 + pr * phi2;
            let pex = exact_soln(x);
            (ph - pex).powi(2)
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
    // Problem: -p'' = 1, p(0)=0, p(1)=0 
    // Exact: p(x) = x * (1.0 - x) / 2.0 (Quadratic)
    // NOTE: Linear elements cannot capture quadratic exactly, but error should be small.
    
    let prob = TestProblem {
        a_fn: |_| 1.0,
        q_fn: |_| 0.0,
        f_fn: |_| 1.0, 
        p_exact: |x| x * (1.0 - x) / 2.0,
    };

    // Use standard penalty of 80.0
    let error = run_solver_and_compute_error(&prob, &prob.p_exact, 10, 80.0);
    
    // Error for N=10 should be around ~1e-3
    assert!(error < 5e-3, "Error was too high for simple Poisson: {}", error);
}

#[test]
fn test_convergence_poisson_sine() {
    // Problem: -p'' = f
    // Target: p(x) = sin(pi * x)
    // f(x) = pi^2 * sin(pi * x)
    
    let prob = TestProblem {
        a_fn: |_| 1.0,
        q_fn: |_| 0.0,
        f_fn: |x| PI.powi(2) * (PI * x).sin(),
        p_exact: |x| (PI * x).sin(),
    };

    let penalty = 10.0;

    let err_coarse = run_solver_and_compute_error(&prob, &prob.p_exact, 10, penalty);
    let err_fine = run_solver_and_compute_error(&prob, &prob.p_exact, 20, penalty);

    println!("Sine Test: Coarse Err = {:.4e}, Fine Err = {:.4e}", err_coarse, err_fine);

    // For linear elements, L2 error is O(h^2). Doubling elements = 1/4 error.
    let rate = err_coarse / err_fine;
    assert!(rate > 3.0, "Convergence rate too slow! Expected ~4.0, got {:.2}", rate);
}

#[test]
fn test_reaction_term() {
    // Problem: -p'' + p = f
    // Target: p(x) = sin(pi * x)
    // f(x) = (pi^2 + 1) * sin(pi * x)
    
    let prob = TestProblem {
        a_fn: |_| 1.0,
        q_fn: |_| 1.0,
        f_fn: |x| (PI.powi(2) + 1.0) * (PI * x).sin(),
        p_exact: |x| (PI * x).sin(),
    };

    let err = run_solver_and_compute_error(&prob, &prob.p_exact, 20, 80.0);
    
    assert!(err < 1e-2, "Reaction-Diffusion error too high: {}", err);
}

#[test]
fn test_variable_coefficient_p() {
    // Problem: -( (1+x) p' )' = f
    // Target: p(x) = x * (1.0 - x)
    // f(x) = 1 + 4x
    
    let prob = TestProblem {
        a_fn: |x| 1.0 + x,
        q_fn: |_| 0.0,
        f_fn: |x| 1.0 + 4.0 * x,
        p_exact: |x| x * (1.0 - x),
    };

    // approximating a quadratic function + variable coefficients.
    let err = run_solver_and_compute_error(&prob, &prob.p_exact, 20, 80.0);

    assert!(err < 5e-3, "Variable coefficient error too high: {}", err);
}
