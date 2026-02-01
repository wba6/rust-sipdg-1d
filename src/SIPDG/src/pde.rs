use util::matrix::Matrix;

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
    pub a: Matrix<f64>,
    pub rhs: Matrix<f64>,
    pub h_elem: Vec<f64>,
    pub x_dof: Vec<f64>,
    pub penalty_param: f64,
}


pub enum Side { Left, Right }

pub trait BoundaryCondition {
    /// Modifies the Global Matrix A and RHS F based on the BC type
    fn apply(
        &self,
        side: Side,
        node_idx: usize,
        h_elem: f64,
        prob: &impl PdeProblem,
        penalty_param: f64,
        a: &mut Matrix<f64>,
        f: &mut Matrix<f64>,
    );
}

pub struct DirichletBC { pub value: f64 }

impl BoundaryCondition for DirichletBC {
    fn apply(&self, side: Side, idx: usize, h: f64, prob: &impl PdeProblem, pen_param: f64, a: &mut Matrix<f64>, f: &mut Matrix<f64>) {
        let n = match side { Side::Left => -1.0, Side::Right => 1.0 };
        let grad_phi = match side { Side::Left => -1.0 / h, Side::Right => 1.0 / h };
        let p_val = prob.p(0.0); // Simplified: should use actual x coord

        let pen = pen_param * (p_val / h);

        // Matrix A contributions (Stability + Consistency + Symmetry)
        a[(idx, idx)] += pen + (2.0 * p_val * grad_phi * n);

        // RHS F contributions (Nitsche enforcement)
        // F += g * [Penalty - (p * grad_u * n)]
        f[(0, idx)] += self.value * (pen - (p_val * grad_phi * n));
    }
}

impl SipdgAssembler {
    pub fn new(h_elem: Vec<f64>, x_dof: Vec<f64>, penalty: f64) -> Self {
        let n_dof = x_dof.len();
        Self {
            a: Matrix::new(n_dof, n_dof, 0.0),
            rhs: Matrix::new(1, n_dof, 0.0),
            h_elem,
            x_dof,
            penalty_param: penalty,
        }
    }

    /// Assembles (Int p u' v') and (Int q u v)
    pub fn assemble_volume(&mut self, prob: &impl PdeProblem) {
        for i in 0..self.h_elem.len() {
            let h_k = self.h_elem[i];
            let idx = [2 * i, 2 * i + 1];
            let xc = (self.x_dof[idx[0]] + self.x_dof[idx[1]]) / 2.0;

            let p_k = prob.p(xc);
            let q_k = prob.q(xc);
            let f_k = prob.f(xc);

            // Local contributions
            let stiff = p_k / h_k;
            let mass = q_k * h_k / 6.0;

            // Apply to global A
            self.a[(idx[0], idx[0])] += stiff + 2.0 * mass;
            self.a[(idx[0], idx[1])] += -stiff + 1.0 * mass;
            self.a[(idx[1], idx[0])] += -stiff + 1.0 * mass;
            self.a[(idx[1], idx[1])] += stiff + 2.0 * mass;

            // Apply to global RHS
            let load = f_k * h_k / 2.0;
            self.rhs[(0, idx[0])] += load;
            self.rhs[(0, idx[1])] += load;
        }
    }

    pub fn assemble_interfaces(&mut self, prob: &impl PdeProblem) {
        for i in 0..self.h_elem.len() - 1 {
            let (idx_l, idx_r) = (2 * i + 1, 2 * i + 2);
            let x_int = self.x_dof[idx_l];
            let p_val = prob.p(x_int);
            let h_avg = 0.5 * (self.h_elem[i] + self.h_elem[i+1]);

            // Penalty
            let pen = self.penalty_param * (p_val / h_avg);
            self.a[(idx_l, idx_l)] += pen;
            self.a[(idx_l, idx_r)] -= pen;
            self.a[(idx_r, idx_l)] -= pen;
            self.a[(idx_r, idx_r)] += pen;

            // Continuity & Symmetry (Simplified logic)
            let g_l = 1.0 / self.h_elem[i];
            let g_r = -1.0 / self.h_elem[i+1];

            self.a[(idx_l, idx_l)] -= p_val * g_l;
            self.a[(idx_l, idx_r)] -= 0.5 * p_val * (g_r - g_l);
            self.a[(idx_r, idx_l)] += 0.5 * p_val * (g_l - g_r);
            self.a[(idx_r, idx_r)] += p_val * g_r;
        }
    }

    pub fn apply_boundaries(
        &mut self,
        prob: &impl PdeProblem,
        left_bc: &impl BoundaryCondition,
        right_bc: &impl BoundaryCondition,
    ) {
        let n_dof = self.x_dof.len();
        
        // Apply to Left Boundary (Node 0)
        left_bc.apply(
            Side::Left, 
            0, 
            self.h_elem[0], 
            prob, 
            self.penalty_param, 
            &mut self.a, 
            &mut self.rhs
        );

        // Apply to Right Boundary (Node n-1)
        right_bc.apply(
            Side::Right, 
            n_dof - 1, 
            *self.h_elem.last().unwrap(), 
            prob, 
            self.penalty_param, 
            &mut self.a, 
            &mut self.rhs
        );
    }
}

