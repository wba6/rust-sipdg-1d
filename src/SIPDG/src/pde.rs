use util::matrix::Matrix;
use util::quadrature::get_gauss_legendre_3pts;

pub struct Element {
    // Global DoF index for the left node
    pub n_left: usize,
    // Global DoF index for the right node
    pub n_right: usize,
    // 2x2 Local stiffness/mass matrix
    pub local_matrix: Matrix<f64>,
    // 1x2 Local right-hand side
    pub local_rhs: Matrix<f64>,
    // Local solution p (populated after global solve)
    pub local_p: Vec<f64>,
    // Affine translation data: stores element width (h_k) and midpoint (x_mid)
    pub h_k: f64,
    pub x_mid: f64,
    // Physical coordinates
    pub x_l: f64,
    pub x_r: f64,
}

impl Element {
    pub fn new(n_left: usize, n_right: usize, x_l: f64, x_r: f64) -> Self {
        let h_k = x_r - x_l;
        Self {
            n_left,
            n_right,
            local_matrix: Matrix::new(2, 2, 0.0),
            local_rhs: Matrix::new(1, 2, 0.0),
            local_p: Vec::new(),
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
        match i {
            0 => 0.5 * (1.0 - xi),
            1 => 0.5 * (1.0 + xi),
            _ => panic!("Invalid basis function index"),
        }
    }

    /// Evaluates the physical derivative of the i-th basis function dphi_i/dx
    pub fn dphi_dx(&self, i: usize) -> f64 {
        let dphi_dxi = match i {
            0 => -0.5,
            1 => 0.5,
            _ => panic!("Invalid basis function index"),
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
    // 4x4 matrix containing contributions to:
    // [L-, R-, L+, R+] x [L-, R-, L+, R+]
    pub matrix: Matrix<f64>,
}

impl Interface {
    pub fn new(x_coord: f64, k_minus: Option<usize>, k_plus: Option<usize>) -> Self {
        Self {
            x_coord,
            k_minus,
            k_plus,
            matrix: Matrix::new(4, 4, 0.0),
        }
    }
}

/// PDE Problem definition: - (a(x) p'(x))' + q(x) p(x) = f(x)
pub trait PdeProblem {
    fn a(&self, x: f64) -> f64;
    fn q(&self, x: f64) -> f64;
    fn f(&self, x: f64) -> f64;
}

pub struct ConfigurableProblem {
    pub a_val: f64,
    pub q_val: f64,
    pub f_val: f64,
}

impl Default for ConfigurableProblem {
    fn default() -> Self {
        Self { a_val: 1.0, q_val: 0.0, f_val: 1.0 }
    }
}

impl PdeProblem for ConfigurableProblem {
    fn a(&self, _x: f64) -> f64 { self.a_val }
    fn q(&self, _x: f64) -> f64 { self.q_val }
    fn f(&self, _x: f64) -> f64 { self.f_val }
}

pub struct SipdgAssembler {
    pub elements: Vec<Element>,
    pub interfaces: Vec<Interface>,
    pub sigma_0: f64, // Penalty parameter sigma^0 in PDF
    pub x_dof: Vec<f64>,
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
        let idx = match side {
            Side::Left => elem.n_left,
            Side::Right => elem.n_right,
        };
        let other_idx = match side {
            Side::Left => elem.n_right,
            Side::Right => elem.n_left,
        };

        let a_val = prob.a(x);
        let pen = assembler.sigma_0 * (a_val / h);

        // Basis function gradients at the boundary
        let grad_phi_i = match side { Side::Left => -1.0 / h, Side::Right => 1.0 / h };
        let grad_phi_j = match side { Side::Left => 1.0 / h, Side::Right => -1.0 / h };

        // For Dirichlet, [p] = -p^+ at x_0 and [p] = p^- at x_N (using PDF Eq 5 signs)
        // However, we handle signs consistently: 
        // LHS term: - {a p'} [v] - {a v'} [p] + (sigma/h) [p] [v]
        // Left boundary (n=0): [v] = -v^+, {v'} = v'^+, [p] = -p^+, {p'} = p'^+
        // Term 1: - p'^+ (-v^+) = p'^+ v^+
        // Term 2: - v'^+ (-p^+) = v'^+ p^+
        // Term 3: (sigma/h) (-p^+) (-v^+) = (sigma/h) p^+ v^+
        
        let n = match side { Side::Left => -1.0, Side::Right => 1.0 };
        // Consistency & Symmetry: - a * grad_p * n * v - a * grad_v * n * p
        // Penalty: (sigma/h) * p * v
        
        // Matrix A contributions (LHS)
        // (p_i, v_i)
        a_global[(idx, idx)] += pen - 2.0 * a_val * grad_phi_i * n;
        // (p_j, v_i)
        a_global[(other_idx, idx)] -= a_val * grad_phi_j * n;
        // (p_i, v_j)
        a_global[(idx, other_idx)] -= a_val * grad_phi_j * n;

        // RHS F contributions (following Eq 18 signs)
        // Left (n=0): + v'(x0) g1 - (sigma/h) g1 v(x0)
        // Right (n=N): - v'(xN) g2 + (sigma/h) g2 v(xN)
        match side {
            Side::Left => {
                f_global[(0, idx)] += self.value * (a_val * grad_phi_i - pen);
                f_global[(0, other_idx)] += self.value * (a_val * grad_phi_j);
            }
            Side::Right => {
                f_global[(0, idx)] += self.value * (-a_val * grad_phi_i + pen);
                f_global[(0, other_idx)] += self.value * (-a_val * grad_phi_j);
            }
        }
    }
}

/// Neumann Boundary Condition: a(x) p'(x) = value
/// Formulation follows PDF Equation (27)
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
        let idx = match side {
            Side::Left => elem.n_left,
            Side::Right => elem.n_right,
        };

        // Neumann terms on RHS: - p'(x0)v(x0) at left, + p'(xN)v(xN) at right
        // Since a(x)p'(x) = value, p'(x) = value / a(x)
        match side {
            Side::Left => {
                f_global[(0, idx)] -= self.value;
            }
            Side::Right => {
                f_global[(0, idx)] += self.value;
            }
        }
    }
}

impl SipdgAssembler {
    pub fn new(h_elem: Vec<f64>, x_dof: Vec<f64>, sigma_0: f64) -> Self {
        let num_elements = h_elem.len();
        let mut elements = Vec::with_capacity(num_elements);
        for i in 0..num_elements {
            elements.push(Element::new(2 * i, 2 * i + 1, x_dof[2 * i], x_dof[2 * i + 1]));
        }

        let mut interfaces = Vec::with_capacity(num_elements - 1);
        for i in 0..num_elements - 1 {
            let x_int = x_dof[2 * i + 1];
            interfaces.push(Interface::new(x_int, Some(i), Some(i + 1)));
        }

        Self {
            elements,
            interfaces,
            sigma_0,
            x_dof,
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

        for elem in &mut self.elements {
            let det_j = elem.jacobian();

            for i in 0..2 {
                for j in 0..2 {
                    let mut stiffness = 0.0;
                    let mut mass = 0.0;
                    
                    let dphi_i_dx = elem.dphi_dx(i);
                    let dphi_j_dx = elem.dphi_dx(j);

                    for pt in &quad_pts {
                        let x = elem.map_to_physical(pt.xi);
                        let w = pt.weight;

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
        }
    }

    pub fn assemble_interfaces(&mut self, prob: &impl PdeProblem) {
        for interface in &mut self.interfaces {
            let k_minus_idx = interface.k_minus.expect("Interface missing left element");
            let k_plus_idx = interface.k_plus.expect("Interface missing right element");

            let x_int = interface.x_coord;
            let a_val = prob.a(x_int);

            let elem_m = &self.elements[k_minus_idx];
            let elem_p = &self.elements[k_plus_idx];

            // Penality coefficient
            let h_avg = Self::average(elem_m.h_k, elem_p.h_k);
            let pen = self.sigma_0 * (a_val / h_avg);

            // Gradients at interface (evaluated on physical elements)
            let g_m = [elem_m.dphi_dx(0), elem_m.dphi_dx(1)];
            let g_p = [elem_p.dphi_dx(0), elem_p.dphi_dx(1)];

            // Basis values at the interface:
            // K- at its right end (xi=1): phi_0=0, phi_1=1
            // K+ at its left end (xi=-1): phi_0=1, phi_1=0
            let v_m = [0.0, 1.0];
            let v_p = [1.0, 0.0];

            for i in 0..4 {
                for j in 0..4 {
                    // Extract basis values and derivatives for trial/test functions
                    let (vi_m, vi_p) = match i {
                        0 => (v_m[0], 0.0), 1 => (v_m[1], 0.0),
                        2 => (0.0, v_p[0]), 3 => (0.0, v_p[1]),
                        _ => unreachable!(),
                    };
                    let (dvi_m, dvi_p) = match i {
                        0 => (g_m[0], 0.0), 1 => (g_m[1], 0.0),
                        2 => (0.0, g_p[0]), 3 => (0.0, g_p[1]),
                        _ => unreachable!(),
                    };

                    let (pj_m, pj_p) = match j {
                        0 => (v_m[0], 0.0), 1 => (v_m[1], 0.0),
                        2 => (0.0, v_p[0]), 3 => (0.0, v_p[1]),
                        _ => unreachable!(),
                    };
                    let (dpj_m, dpj_p) = match j {
                        0 => (g_m[0], 0.0), 1 => (g_m[1], 0.0),
                        2 => (0.0, g_p[0]), 3 => (0.0, g_p[1]),
                        _ => unreachable!(),
                    };

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
        }
    }

    pub fn assemble_to_global(&self) -> (Matrix<f64>, Matrix<f64>) {
        let n_dof = self.elements.len() * 2;
        let mut a = Matrix::new(n_dof, n_dof, 0.0);
        let mut rhs = Matrix::new(1, n_dof, 0.0);

        // Volume terms
        for elem in &self.elements {
            let indices = [elem.n_left, elem.n_right];
            for i in 0..2 {
                for j in 0..2 {
                    a[(indices[i], indices[j])] += elem.local_matrix[(i, j)];
                }
                rhs[(0, indices[i])] += elem.local_rhs[(0, i)];
            }
        }

        // Interface terms
        for interface in &self.interfaces {
            let k_minus = interface.k_minus.expect("Interface missing left element");
            let k_plus = interface.k_plus.expect("Interface missing right element");
            
            let nodes = [
                self.elements[k_minus].n_left,
                self.elements[k_minus].n_right,
                self.elements[k_plus].n_left,
                self.elements[k_plus].n_right,
            ];

            for i in 0..4 {
                for j in 0..4 {
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
}
