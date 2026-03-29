use util::matrix::Matrix;
use util::quadrature::{get_gauss_legendre_3pts, integrate};

pub struct Element {
    // Global DoF index for the left node
    pub n_left: usize,
    // Global DoF index for the right node
    pub n_right: usize,
    // 2x2 Local stiffness/mass matrix
    pub local_matrix: Matrix<f64>,
    // 1x2 Local right-hand side
    pub local_rhs: Matrix<f64>,
    // Local solution (populated after global solve)
    pub local_soln: Vec<f64>,
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
            local_soln: Vec::new(),
            h_k,
            x_mid: (x_l + x_r) / 2.0,
            x_l,
            x_r,
        }
    }
}

pub struct Interface {
    // Physical x coordinate of the interface
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

pub trait PdeProblem {
    fn p(&self, x: f64) -> f64;
    fn q(&self, x: f64) -> f64;
    fn f(&self, x: f64) -> f64;
}

pub struct ConfigurableProblem {
    pub p_val: f64,
    pub q_val: f64,
    pub f_val: f64,
}

impl Default for ConfigurableProblem {
    fn default() -> Self {
        Self { p_val: 1.0, q_val: 0.0, f_val: 1.0 }
    }
}

impl PdeProblem for ConfigurableProblem {
    fn p(&self, _x: f64) -> f64 { self.p_val }
    fn q(&self, _x: f64) -> f64 { self.q_val }
    fn f(&self, _x: f64) -> f64 { self.f_val }
}

pub struct SipdgAssembler {
    pub elements: Vec<Element>,
    pub interfaces: Vec<Interface>,
    pub penalty_param: f64,
    pub x_dof: Vec<f64>,
}

pub enum Side { Left, Right }

pub trait BoundaryCondition {
    /// Modifies the Global Matrix A and RHS F based on the BC type
    fn apply(
        &self,
        side: Side,
        elem_idx: usize,
        assembler: &mut SipdgAssembler,
        prob: &impl PdeProblem,
        a_global: &mut Matrix<f64>,
        f_global: &mut Matrix<f64>,
    );
}

pub struct DirichletBC { pub value: f64 }

impl BoundaryCondition for DirichletBC {
    fn apply(
        &self,
        side: Side,
        elem_idx: usize,
        assembler: &mut SipdgAssembler,
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

        let n = match side { Side::Left => -1.0, Side::Right => 1.0 };
        // Grad phi for the boundary node
        let grad_phi = match side { Side::Left => -1.0 / h, Side::Right => 1.0 / h };
        // Grad phi for the OTHER node in the same element
        let grad_phi_other = match side { Side::Left => 1.0 / h, Side::Right => -1.0 / h };
        let other_idx = match side { Side::Left => elem.n_right, Side::Right => elem.n_left };

        let p_val = prob.p(x);
        let pen = assembler.penalty_param * (p_val / h);

        // Matrix A contributions (Stability + Consistency + Symmetry)
        // (u_idx, v_idx)
        a_global[(idx, idx)] += pen - (2.0 * p_val * grad_phi * n);
        
        // (u_other, v_idx)
        a_global[(other_idx, idx)] -= p_val * grad_phi_other * n;
        
        // (u_idx, v_other)
        a_global[(idx, other_idx)] -= p_val * grad_phi_other * n;

        // RHS F contributions (Nitsche enforcement)
        f_global[(0, idx)] += self.value * (pen - (p_val * grad_phi * n));
        f_global[(0, other_idx)] -= self.value * (p_val * grad_phi_other * n);
    }
}

impl SipdgAssembler {
    pub fn new(h_elem: Vec<f64>, x_dof: Vec<f64>, penalty: f64) -> Self {
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
            penalty_param: penalty,
            x_dof,
        }
    }

    /// Assembles (Int p u' v') and (Int q u v) and (Int f v) using Gauss-Legendre quadrature
    pub fn assemble_volume(&mut self, prob: &impl PdeProblem) {
        let quad_pts = get_gauss_legendre_3pts();

        for elem in &mut self.elements {
            let h = elem.h_k;
            let x_l = elem.x_l;
            let x_r = elem.x_r;

            // Assemble local matrix (2x2)
            for i in 0..2 {
                for j in 0..2 {
                    let dphi_i = if i == 0 { -1.0 / h } else { 1.0 / h };
                    let dphi_j = if j == 0 { -1.0 / h } else { 1.0 / h };

                    // p(x) u' v' term
                    let stiff = integrate(|x| prob.p(x) * dphi_i * dphi_j, x_l, x_r, &quad_pts);

                    // q(x) u v term
                    let mass = integrate(
                        |x| {
                            let phi_i = if i == 0 { (x_r - x) / h } else { (x - x_l) / h };
                            let phi_j = if j == 0 { (x_r - x) / h } else { (x - x_l) / h };
                            prob.q(x) * phi_i * phi_j
                        },
                        x_l,
                        x_r,
                        &quad_pts,
                    );

                    elem.local_matrix[(i, j)] = stiff + mass;
                }

                // f(x) v term
                let load = integrate(
                    |x| {
                        let phi_i = if i == 0 { (x_r - x) / h } else { (x - x_l) / h };
                        prob.f(x) * phi_i
                    },
                    x_l,
                    x_r,
                    &quad_pts,
                );
                elem.local_rhs[(0, i)] = load;
            }
        }
    }

    pub fn assemble_interfaces(&mut self, prob: &impl PdeProblem) {
        for interface in &mut self.interfaces {
            let k_minus_idx = interface.k_minus.expect("Interface missing left element");
            let k_plus_idx = interface.k_plus.expect("Interface missing right element");

            let x_int = interface.x_coord;
            let p_val = prob.p(x_int);

            let h_l = self.elements[k_minus_idx].h_k;
            let h_r = self.elements[k_plus_idx].h_k;
            let h_avg = 0.5 * (h_l + h_r);

            let pen = self.penalty_param * (p_val / h_avg);

            // Gradients at interface:
            // K-: phi_L' = -1/h_l, phi_R' = 1/h_l
            // K+: phi_L' = -1/h_r, phi_R' = 1/h_r
            let g_l0 = -1.0 / h_l;
            let g_l1 = 1.0 / h_l;
            let g_r0 = -1.0 / h_r;
            let g_r1 = 1.0 / h_r;

            // [u] = u_R^- - u_L^+ = u_1 - u_2  (using indices 0..3 for [L-, R-, L+, R+])
            // {p u'} = 0.5 * p * ( (u_0*g_l0 + u_1*g_l1) + (u_2*g_r0 + u_3*g_r1) )
            
            // Term: - {p u'} [v] - {p v'} [u] + pen [u] [v]
            // where [v] = v_1 - v_2
            
            let c = 0.5 * p_val;

            // Mapping: 0->L-, 1->R-, 2->L+, 3->R+
            
            // Penalty: pen * (u1 - u2) * (v1 - v2)
            interface.matrix[(1, 1)] += pen;
            interface.matrix[(1, 2)] -= pen;
            interface.matrix[(2, 1)] -= pen;
            interface.matrix[(2, 2)] += pen;

            // Consistency: - 0.5 * p * (u0*g_l0 + u1*g_l1 + u2*g_r0 + u3*g_r1) * (v1 - v2)
            for u_node in 0..4 {
                let grad = match u_node {
                    0 => g_l0,
                    1 => g_l1,
                    2 => g_r0,
                    3 => g_r1,
                    _ => unreachable!(),
                };
                interface.matrix[(1, u_node)] -= c * grad; // v1 part
                interface.matrix[(2, u_node)] += c * grad; // v2 part
            }

            // Symmetry: - 0.5 * p * (v0*g_l0 + v1*g_l1 + v2*g_r0 + v3*g_r1) * (u1 - u2)
            for v_node in 0..4 {
                let grad = match v_node {
                    0 => g_l0,
                    1 => g_l1,
                    2 => g_r0,
                    3 => g_r1,
                    _ => unreachable!(),
                };
                interface.matrix[(v_node, 1)] -= c * grad; // u1 part
                interface.matrix[(v_node, 2)] += c * grad; // u2 part
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
        &mut self,
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
