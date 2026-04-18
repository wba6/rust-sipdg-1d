use util::matrix::Matrix;
use util::quadrature::get_gauss_legendre_3pts;
use util::cg::LinearOperator;
use rayon::prelude::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BasisOrder {
    Linear,
    Quadratic,
}

impl BasisOrder {
    pub fn num_nodes(&self) -> usize {
        match self {
            BasisOrder::Linear => 2,
            BasisOrder::Quadratic => 3,
        }
    }
}

pub struct Element {
    // Global DoF indices for the nodes
    pub nodes: Vec<usize>,
    // Local stiffness/mass matrix
    pub local_matrix: Matrix<f64>,
    // Local right-hand side
    pub local_rhs: Matrix<f64>,
    // Basis order
    pub order: BasisOrder,
    // Affine translation data: stores element width (h_k) and midpoint (x_mid)
    pub h_k: f64,
    pub x_mid: f64,
    // Physical coordinates
    pub x_l: f64,
    pub x_r: f64,
}

impl Element {
    pub fn new(nodes: Vec<usize>, x_l: f64, x_r: f64, order: BasisOrder) -> Self {
        let h_k = x_r - x_l;
        let n_nodes = order.num_nodes();
        Self {
            nodes,
            local_matrix: Matrix::new(n_nodes, n_nodes, 0.0),
            local_rhs: Matrix::new(1, n_nodes, 0.0),
            order,
            h_k,
            x_mid: (x_l + x_r) / 2.0,
            x_l,
            x_r,
        }
    }

    /// Maps a reference coordinate xi in [-1, 1] to the physical coordinate x in [x_l, x_r]
    pub fn map_to_physical(&self, xi: f64) -> f64 {
        self.x_mid + (self.h_k / 2.0) * xi
    }

    /// Returns the Jacobian J = dx/dxi = h_k / 2
    pub fn jacobian(&self) -> f64 {
        self.h_k / 2.0
    }

    /// Evaluates the i-th basis function at reference coordinate xi
    pub fn phi(&self, i: usize, xi: f64) -> f64 {
        match self.order {
            BasisOrder::Linear => match i {
                0 => 0.5 * (1.0 - xi),
                1 => 0.5 * (1.0 + xi),
                _ => panic!("Invalid linear basis function index"),
            },
            BasisOrder::Quadratic => match i {
                0 => 0.5 * xi * (xi - 1.0),
                1 => 1.0 - xi * xi,
                2 => 0.5 * xi * (xi + 1.0),
                _ => panic!("Invalid quadratic basis function index"),
            },
        }
    }

    /// Evaluates the physical derivative of the i-th basis function dphi_i/dx
    pub fn dphi_dx(&self, i: usize, xi: f64) -> f64 {
        let dphi_dxi = match self.order {
            BasisOrder::Linear => match i {
                0 => -0.5,
                1 => 0.5,
                _ => panic!("Invalid linear basis function index"),
            },
            BasisOrder::Quadratic => match i {
                0 => xi - 0.5,
                1 => -2.0 * xi,
                2 => xi + 0.5,
                _ => panic!("Invalid quadratic basis function index"),
            },
        };
        dphi_dxi / self.jacobian()
    }
}

pub struct Interface {
    // Physical x coordinate of the interface x_n
    pub x_coord: f64,
    // Index of the element to the left (K-)
    pub k_minus: Option<usize>,
    // Index of the element to the right (K+)
    pub k_plus: Option<usize>,
    // Interaction matrix for the interface
    pub matrix: Matrix<f64>,
}

impl Interface {
    pub fn new(x_coord: f64, k_minus: Option<usize>, k_plus: Option<usize>, order: BasisOrder) -> Self {
        let n_side = order.num_nodes();
        Self {
            x_coord,
            k_minus,
            k_plus,
            matrix: Matrix::new(2 * n_side, 2 * n_side, 0.0),
        }
    }
}

/// PDE Problem definition: - (a(x) p'(x))' + q(x) p(x) = f(x)
pub trait PdeProblem: Sync {
    fn a(&self, x: f64) -> f64;
    fn q(&self, x: f64) -> f64;
    fn f(&self, x: f64) -> f64;
}

pub struct ConfigurableProblem {
    pub a_val: f64,
    pub q_val: f64,
    pub f_val: f64,
    pub num_elements: usize,
    pub sigma_0: f64,
}

impl Default for ConfigurableProblem {
    fn default() -> Self {
        Self { 
            a_val: 1.0, 
            q_val: 1.0, 
            f_val: 1.0, // Use 1.0 as a flag for default non-constant f(x) or just a constant
            num_elements: 3000,
            sigma_0: 20.0,
        }
    }
}

impl PdeProblem for ConfigurableProblem {
    fn a(&self, _x: f64) -> f64 { self.a_val }
    fn q(&self, _x: f64) -> f64 { self.q_val }
    fn f(&self, x: f64) -> f64 { 
        // If q=1, a=1 and f=1, we assume the user wants the "full" SL default: (pi^2 + 1) * sin(pi * x)
        if self.a_val == 1.0 && self.q_val == 1.0 && self.f_val == 1.0 {
            use std::f64::consts::PI;
            (PI.powi(2) + 1.0) * (PI * x).sin()
        } else {
            self.f_val 
        }
    }
}

pub struct SipdgAssembler {
    pub elements: Vec<Element>,
    pub interfaces: Vec<Interface>,
    pub sigma_0: f64, // Penalty parameter sigma^0 in PDF
    pub x_dof: Vec<f64>,
    pub order: BasisOrder,
}

pub enum Side { Left, Right }

pub trait BoundaryCondition {
    /// Modifies the Global Matrix A and RHS F based on the BC type
    fn apply(
        &self,
        side: Side,
        elem_idx: usize,
        assembler: &SipdgAssembler,
        prob: &impl PdeProblem,
        a_global: &mut Matrix<f64>,
        f_global: &mut Matrix<f64>,
    );

    /// Matrix-free contribution to A*v
    fn apply_matrix_free(
        &self,
        side: Side,
        elem_idx: usize,
        assembler: &SipdgAssembler,
        prob: &impl PdeProblem,
        v_in: &[f64],
        v_out: &mut [f64],
    );
}

/// Dirichlet Boundary Condition: p(x) = value
/// Formulation follows PDF Equation (18)
pub struct DirichletBC { pub value: f64 }

impl BoundaryCondition for DirichletBC {
    fn apply(
        &self,
        side: Side,
        elem_idx: usize,
        assembler: &SipdgAssembler,
        prob: &impl PdeProblem,
        a_global: &mut Matrix<f64>,
        f_global: &mut Matrix<f64>,
    ) {
        let elem = &assembler.elements[elem_idx];
        let h = elem.h_k;
        let x = match side {
            Side::Left => elem.x_l,
            Side::Right => elem.x_r,
        };
        let xi = match side {
            Side::Left => -1.0,
            Side::Right => 1.0,
        };
        let n = match side {
            Side::Left => -1.0,
            Side::Right => 1.0,
        };

        let a_val = prob.a(x);
        let pen = assembler.sigma_0 * (a_val / h);
        let n_nodes = elem.order.num_nodes();

        // Matrix A contributions (LHS)
        for i in 0..n_nodes {
            let row_idx = elem.nodes[i];
            let vi = elem.phi(i, xi);
            let dvi_dx = elem.dphi_dx(i, xi);

            for j in 0..n_nodes {
                let col_idx = elem.nodes[j];
                let pj = elem.phi(j, xi);
                let dpj_dx = elem.dphi_dx(j, xi);

                // Terms: pen * p * v - a * grad_p * n * v - a * grad_v * n * p
                a_global[(row_idx, col_idx)] += pen * pj * vi 
                                              - a_val * dpj_dx * n * vi 
                                              - a_val * dvi_dx * n * pj;
            }

            // RHS F contributions
            match side {
                Side::Left => {
                    // + v'(x0) g1 - (sigma/h) g1 v(x0)
                    f_global[(0, row_idx)] += self.value * (a_val * dvi_dx - pen * vi);
                }
                Side::Right => {
                    // - v'(xN) g2 + (sigma/h) g2 v(xN)
                    f_global[(0, row_idx)] += self.value * (-a_val * dvi_dx + pen * vi);
                }
            }
        }
    }

    fn apply_matrix_free(
        &self,
        side: Side,
        elem_idx: usize,
        assembler: &SipdgAssembler,
        prob: &impl PdeProblem,
        v_in: &[f64],
        v_out: &mut [f64],
    ) {
        let elem = &assembler.elements[elem_idx];
        let h = elem.h_k;
        let x = match side {
            Side::Left => elem.x_l,
            Side::Right => elem.x_r,
        };
        let xi = match side {
            Side::Left => -1.0,
            Side::Right => 1.0,
        };
        let n = match side {
            Side::Left => -1.0,
            Side::Right => 1.0,
        };

        let a_val = prob.a(x);
        let pen = assembler.sigma_0 * (a_val / h);
        let n_nodes = elem.order.num_nodes();

        for i in 0..n_nodes {
            let row_idx = elem.nodes[i];
            let vi = elem.phi(i, xi);
            let dvi_dx = elem.dphi_dx(i, xi);

            let mut sum = 0.0;
            for j in 0..n_nodes {
                let col_idx = elem.nodes[j];
                let pj = elem.phi(j, xi);
                let dpj_dx = elem.dphi_dx(j, xi);

                sum += (pen * pj * vi 
                      - a_val * dpj_dx * n * vi 
                      - a_val * dvi_dx * n * pj) * v_in[col_idx];
            }
            v_out[row_idx] += sum;
        }
    }
}

/// Neumann Boundary Condition: a(x) p'(x) = value
/// Formulation follows PDF Equation (27)
#[allow(dead_code)]
pub struct NeumannBC { pub value: f64 }

impl BoundaryCondition for NeumannBC {
    fn apply(
        &self,
        side: Side,
        elem_idx: usize,
        assembler: &SipdgAssembler,
        _prob: &impl PdeProblem,
        _a_global: &mut Matrix<f64>,
        f_global: &mut Matrix<f64>,
    ) {
        let elem = &assembler.elements[elem_idx];
        let xi = match side {
            Side::Left => -1.0,
            Side::Right => 1.0,
        };
        let n_nodes = elem.order.num_nodes();

        for i in 0..n_nodes {
            let idx = elem.nodes[i];
            let vi = elem.phi(i, xi);
            // Neumann terms on RHS: - p'(x0)v(x0) at left, + p'(xN)v(xN) at right
            // Since a(x)p'(x) = value, p'(x) = value / a(x)
            match side {
                Side::Left => {
                    f_global[(0, idx)] -= self.value * vi;
                }
                Side::Right => {
                    f_global[(0, idx)] += self.value * vi;
                }
            }
        }
    }

    fn apply_matrix_free(
        &self,
        _side: Side,
        _elem_idx: usize,
        _assembler: &SipdgAssembler,
        _prob: &impl PdeProblem,
        _v_in: &[f64],
        _v_out: &mut [f64],
    ) {
        // Neumann BC only affects RHS, not the operator A
    }
}

impl SipdgAssembler {
    pub fn new(h_elem: Vec<f64>, x_dof: Vec<f64>, sigma_0: f64, order: BasisOrder) -> Self {
        let num_elements = h_elem.len();
        let n_nodes = order.num_nodes();
        let mut elements = Vec::with_capacity(num_elements);
        for i in 0..num_elements {
            let mut nodes = Vec::with_capacity(n_nodes);
            for j in 0..n_nodes {
                nodes.push(i * n_nodes + j);
            }
            elements.push(Element::new(nodes, x_dof[i * n_nodes], x_dof[(i + 1) * n_nodes - 1], order));
        }

        let mut interfaces = Vec::with_capacity(num_elements - 1);
        for i in 0..num_elements - 1 {
            let x_int = x_dof[(i + 1) * n_nodes - 1];
            interfaces.push(Interface::new(x_int, Some(i), Some(i + 1), order));
        }

        Self {
            elements,
            interfaces,
            sigma_0,
            x_dof,
            order,
        }
    }

    /// Jump: [v] = v^- - v^+
    fn jump(v_minus: f64, v_plus: f64) -> f64 {
        v_minus - v_plus
    }

    /// Average: {v} = 0.5 * (v^- + v^+)
    fn average(v_minus: f64, v_plus: f64) -> f64 {
        0.5 * (v_minus + v_plus)
    }

    /// Assembles (Int a p' v') and (Int q p v) and (Int f v) using the reference element [-1, 1]
    pub fn assemble_volume(&mut self, prob: &impl PdeProblem) {
        let quad_pts = get_gauss_legendre_3pts();
        let n_nodes = self.order.num_nodes();

        self.elements.par_iter_mut().for_each(|elem| {
            let det_j = elem.jacobian();

            for i in 0..n_nodes {
                for j in 0..n_nodes {
                    let mut stiffness = 0.0;
                    let mut mass = 0.0;
                    
                    for pt in &quad_pts {
                        let x = elem.map_to_physical(pt.xi);
                        let w = pt.weight;

                        let dphi_i_dx = elem.dphi_dx(i, pt.xi);
                        let dphi_j_dx = elem.dphi_dx(j, pt.xi);

                        // a(x) p' v'
                        stiffness += w * prob.a(x) * dphi_i_dx * dphi_j_dx;

                        // q(x) p v
                        mass += w * prob.q(x) * elem.phi(i, pt.xi) * elem.phi(j, pt.xi);
                    }

                    elem.local_matrix[(i, j)] = (stiffness + mass) * det_j;
                }

                // f(x) v
                let mut load = 0.0;
                for pt in &quad_pts {
                    let x = elem.map_to_physical(pt.xi);
                    load += pt.weight * prob.f(x) * elem.phi(i, pt.xi);
                }
                elem.local_rhs[(0, i)] = load * det_j;
            }
        });
    }

    pub fn assemble_interfaces(&mut self, prob: &impl PdeProblem) {
        let n_side = self.order.num_nodes();
        let elements = &self.elements;
        self.interfaces.par_iter_mut().for_each(|interface| {
            let k_minus_idx = interface.k_minus.expect("Interface missing left element");
            let k_plus_idx = interface.k_plus.expect("Interface missing right element");

            let x_int = interface.x_coord;
            let a_val = prob.a(x_int);

            let elem_m = &elements[k_minus_idx];
            let elem_p = &elements[k_plus_idx];

            // Penality coefficient
            let h_avg = Self::average(elem_m.h_k, elem_p.h_k);
            let pen = self.sigma_0 * (a_val / h_avg);

            // Gradients and values at interface
            // K- at its right end (xi=1)
            let mut v_m = Vec::with_capacity(n_side);
            let mut g_m = Vec::with_capacity(n_side);
            for i in 0..n_side {
                v_m.push(elem_m.phi(i, 1.0));
                g_m.push(elem_m.dphi_dx(i, 1.0));
            }

            // K+ at its left end (xi=-1)
            let mut v_p = Vec::with_capacity(n_side);
            let mut g_p = Vec::with_capacity(n_side);
            for i in 0..n_side {
                v_p.push(elem_p.phi(i, -1.0));
                g_p.push(elem_p.dphi_dx(i, -1.0));
            }

            for i in 0..2 * n_side {
                for j in 0..2 * n_side {
                    // i and j are indices into [nodes_m, nodes_p]
                    let (vi_m, vi_p) = if i < n_side { (v_m[i], 0.0) } else { (0.0, v_p[i - n_side]) };
                    let (dvi_m, dvi_p) = if i < n_side { (g_m[i], 0.0) } else { (0.0, g_p[i - n_side]) };

                    let (pj_m, pj_p) = if j < n_side { (v_m[j], 0.0) } else { (0.0, v_p[j - n_side]) };
                    let (dpj_m, dpj_p) = if j < n_side { (g_m[j], 0.0) } else { (0.0, g_p[j - n_side]) };

                    // Terms: pen [p] [v] - {a p'} [v] - {a v'} [p]
                    let jump_v = Self::jump(vi_m, vi_p);
                    let jump_p = Self::jump(pj_m, pj_p);
                    let avg_dp = Self::average(dpj_m, dpj_p);
                    let avg_dv = Self::average(dvi_m, dvi_p);

                    interface.matrix[(i, j)] += pen * jump_p * jump_v
                                              - a_val * avg_dp * jump_v
                                              - a_val * avg_dv * jump_p;
                }
            }
        });
    }

    pub fn assemble_to_global(&self) -> (Matrix<f64>, Matrix<f64>) {
        let n_dof = self.x_dof.len();
        let mut a = Matrix::new(n_dof, n_dof, 0.0);
        let mut rhs = Matrix::new(1, n_dof, 0.0);
        let n_nodes = self.order.num_nodes();

        // Volume terms
        for elem in &self.elements {
            for i in 0..n_nodes {
                for j in 0..n_nodes {
                    a[(elem.nodes[i], elem.nodes[j])] += elem.local_matrix[(i, j)];
                }
                rhs[(0, elem.nodes[i])] += elem.local_rhs[(0, i)];
            }
        }

        // Interface terms
        for interface in &self.interfaces {
            let k_minus_idx = interface.k_minus.expect("Interface missing left element");
            let k_plus_idx = interface.k_plus.expect("Interface missing right element");
            
            let mut nodes = Vec::<usize>::with_capacity(2 * n_nodes);
            nodes.extend(&self.elements[k_minus_idx].nodes);
            nodes.extend(&self.elements[k_plus_idx].nodes);

            for i in 0..2 * n_nodes {
                for j in 0..2 * n_nodes {
                    a[(nodes[i], nodes[j])] += interface.matrix[(i, j)];
                }
            }
        }

        (a, rhs)
    }

    pub fn apply_boundaries(
        &self,
        prob: &impl PdeProblem,
        left_bc: &impl BoundaryCondition,
        right_bc: &impl BoundaryCondition,
        a_global: &mut Matrix<f64>,
        f_global: &mut Matrix<f64>,
    ) {
        // Apply to Left Boundary (Element 0, Left side)
        left_bc.apply(
            Side::Left,
            0,
            self,
            prob,
            a_global,
            f_global
        );

        // Apply to Right Boundary (Last Element, Right side)
        right_bc.apply(
            Side::Right,
            self.elements.len() - 1,
            self,
            prob,
            a_global,
            f_global
        );
    }

    pub fn matrix_free_op<'a, P, B1, B2>(
        &'a self,
        prob: &'a P,
        left_bc: &'a B1,
        right_bc: &'a B2,
    ) -> SipdgOperator<'a, P, B1, B2>
    where
        P: PdeProblem,
        B1: BoundaryCondition,
        B2: BoundaryCondition,
    {
        SipdgOperator {
            assembler: self,
            prob,
            left_bc,
            right_bc,
        }
    }
}

pub struct SipdgOperator<'a, P, B1, B2>
where
    P: PdeProblem,
    B1: BoundaryCondition,
    B2: BoundaryCondition,
{
    assembler: &'a SipdgAssembler,
    prob: &'a P,
    left_bc: &'a B1,
    right_bc: &'a B2,
}

impl<'a, P, B1, B2> LinearOperator for SipdgOperator<'a, P, B1, B2>
where
    P: PdeProblem,
    B1: BoundaryCondition,
    B2: BoundaryCondition,
{
    fn columns(&self) -> usize {
        self.assembler.x_dof.len()
    }

    fn multiply_vec(&self, v: &[f64]) -> Vec<f64> {
        use std::sync::atomic::{AtomicU64, Ordering};

        let n = self.columns();
        let out_atomic: Vec<AtomicU64> = (0..n)
            .map(|_| AtomicU64::new(0.0f64.to_bits()))
            .collect();
        let n_nodes = self.assembler.order.num_nodes();

        // Helper to add f64 to AtomicU64
        let add_to_atomic = |idx: usize, val: f64| {
            let mut current_bits = out_atomic[idx].load(Ordering::Relaxed);
            loop {
                let current_val = f64::from_bits(current_bits);
                let new_val = current_val + val;
                let new_bits = new_val.to_bits();
                match out_atomic[idx].compare_exchange_weak(
                    current_bits,
                    new_bits,
                    Ordering::SeqCst,
                    Ordering::Relaxed,
                ) {
                    Ok(_) => break,
                    Err(actual_bits) => current_bits = actual_bits,
                }
            }
        };

        // Volume contributions
        self.assembler.elements.par_iter().for_each(|elem| {
            for i in 0..n_nodes {
                let row_idx = elem.nodes[i];
                let mut sum = 0.0;
                for j in 0..n_nodes {
                    let col_idx = elem.nodes[j];
                    sum += elem.local_matrix[(i, j)] * v[col_idx];
                }
                add_to_atomic(row_idx, sum);
            }
        });

        // Interface contributions
        self.assembler.interfaces.par_iter().for_each(|interface| {
            let k_m = interface.k_minus.expect("Interface missing left element");
            let k_p = interface.k_plus.expect("Interface missing right element");

            let mut nodes = Vec::<usize>::with_capacity(2 * n_nodes);
            nodes.extend(&self.assembler.elements[k_m].nodes);
            nodes.extend(&self.assembler.elements[k_p].nodes);

            for i in 0..2 * n_nodes {
                let row_idx = nodes[i];
                let mut sum = 0.0;
                for j in 0..2 * n_nodes {
                    let col_idx = nodes[j];
                    sum += interface.matrix[(i, j)] * v[col_idx];
                }
                add_to_atomic(row_idx, sum);
            }
        });

        // Boundary contributions (typically small, can stay serial or use atomics)
        let mut b_out = vec![0.0; n];
        self.left_bc.apply_matrix_free(
            Side::Left,
            0,
            self.assembler,
            self.prob,
            v,
            &mut b_out,
        );
        self.right_bc.apply_matrix_free(
            Side::Right,
            self.assembler.elements.len() - 1,
            self.assembler,
            self.prob,
            v,
            &mut b_out,
        );

        // Combine atomics and boundary
        out_atomic
            .into_par_iter()
            .enumerate()
            .map(|(i, val)| f64::from_bits(val.load(Ordering::SeqCst)) + b_out[i])
            .collect()
    }
}
