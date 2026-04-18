// src/SIPDG/tests/solver_tests.rs

use sipdg::{PdeProblem, SipdgAssembler, DirichletBC, generate_mesh, pde::BasisOrder};
use util::cg::cg;
use std::f64::consts::PI;

/// A flexible problem struct that lets us define a, q, f, and exact_p using closures.
struct TestProblem<A, Q, F, P>
where
    A: Fn(f64) -> f64 + Sync,
    Q: Fn(f64) -> f64 + Sync,
    F: Fn(f64) -> f64 + Sync,
    P: Fn(f64) -> f64 + Sync,
{
    a_fn: A,
    q_fn: Q,
    f_fn: F,
    p_exact: P,
}

impl<A, Q, F, P> PdeProblem for TestProblem<A, Q, F, P>
where
    A: Fn(f64) -> f64 + Sync,
    Q: Fn(f64) -> f64 + Sync,
    F: Fn(f64) -> f64 + Sync,
    P: Fn(f64) -> f64 + Sync,
{
    fn a(&self, x: f64) -> f64 { (self.a_fn)(x) }
    fn q(&self, x: f64) -> f64 { (self.q_fn)(x) }
    fn f(&self, x: f64) -> f64 { (self.f_fn)(x) }
}

fn run_solver_and_compute_error(
    prob: &impl PdeProblem, 
    exact_soln: impl Fn(f64) -> f64, 
    num_elements: usize,
    penalty: f64
) -> f64 {
    let order = BasisOrder::Linear;
    // Generate Mesh
    let (h_elem, x_dof) = generate_mesh(0.0, 1.0, num_elements, order);
    
    // Assemble (using the provided penalty parameter)
    let mut assembler = SipdgAssembler::new(h_elem.clone(), x_dof.clone(), penalty, order);
    assembler.assemble_volume(prob);
    assembler.assemble_interfaces(prob);

    let (mut a, mut rhs) = assembler.assemble_to_global();

    // Apply BCs (assuming homogeneous Dirichlet p(0)=p(1)=0 for these tests)
    let bc = DirichletBC { value: 0.0 };
    assembler.apply_boundaries(prob, &bc, &bc, &mut a, &mut rhs);

    // Solve (Matrix-Free)
    let op = assembler.matrix_free_op(prob, &bc, &bc);
    let p = cg(&op, &rhs, 1e-10, 10000);

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
fn test_convergence_rates() {
    let prob = TestProblem {
        a_fn: |_| 1.0,
        q_fn: |_| 1.0,
        f_fn: |x| (PI.powi(2) + 1.0) * (PI * x).sin(),
        p_exact: |x| (PI * x).sin(),
    };

    let penalty = 20.0;
    let sizes = [10, 20, 40];

    // Linear
    let mut errors_p1 = Vec::new();
    for &n in &sizes {
        errors_p1.push(run_solver_and_compute_error_with_order(&prob, &prob.p_exact, n, penalty, BasisOrder::Linear));
    }
    
    for i in 1..errors_p1.len() {
        let rate = (errors_p1[i-1] / errors_p1[i]).log2();
        assert!(rate > 1.9, "Linear convergence rate too low: {:.2}", rate);
    }

    // Quadratic
    let mut errors_p2 = Vec::new();
    for &n in &sizes {
        errors_p2.push(run_solver_and_compute_error_with_order(&prob, &prob.p_exact, n, penalty, BasisOrder::Quadratic));
    }
    
    for i in 1..errors_p2.len() {
        let rate = (errors_p2[i-1] / errors_p2[i]).log2();
        assert!(rate > 2.8, "Quadratic convergence rate too low: {:.2}", rate);
    }
}

fn run_solver_and_compute_error_with_order(
    prob: &impl PdeProblem, 
    exact_soln: impl Fn(f64) -> f64, 
    num_elements: usize,
    penalty: f64,
    order: BasisOrder
) -> f64 {
    let (h_elem, x_dof) = generate_mesh(0.0, 1.0, num_elements, order);
    let mut assembler = SipdgAssembler::new(h_elem.clone(), x_dof.clone(), penalty, order);
    assembler.assemble_volume(prob);
    assembler.assemble_interfaces(prob);

    let (mut a, mut rhs) = assembler.assemble_to_global();
    let bc = DirichletBC { value: 0.0 };
    assembler.apply_boundaries(prob, &bc, &bc, &mut a, &mut rhs);

    let op = assembler.matrix_free_op(prob, &bc, &bc);
    let p = cg(&op, &rhs, 1e-12, 10000);

    // Compute L2 Error (Trapezoidal rule used in tests for simplicity, 
    // but with enough points it should be fine for rate verification)
    let mut l2_err_sq = 0.0;
    let npts = 20;

    for k in 0..num_elements {
        let n_local = order.num_nodes();
        let x_l = assembler.elements[k].x_l;
        let x_r = assembler.elements[k].x_r;
        let h = x_r - x_l;

        let step = h / ((npts - 1) as f64);
        
        let get_ph = |x: f64| {
            let xi = (x - assembler.elements[k].x_mid) / (h / 2.0);
            let mut val = 0.0;
            for i in 0..n_local {
                val += p[assembler.elements[k].nodes[i]] * assembler.elements[k].phi(i, xi);
            }
            val
        };

        let mut prev_x = x_l;
        let mut prev_err2 = (get_ph(prev_x) - exact_soln(prev_x)).powi(2);

        for i in 1..npts {
            let x = x_l + step * (i as f64);
            let err2 = (get_ph(x) - exact_soln(x)).powi(2);
            l2_err_sq += 0.5 * (x - prev_x) * (prev_err2 + err2);
            prev_x = x;
            prev_err2 = err2;
        }
    }

    l2_err_sq.sqrt()
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
