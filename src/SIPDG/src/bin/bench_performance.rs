use sipdg::{PdeProblem, SipdgAssembler, DirichletBC, NeumannBC, generate_mesh, pde::BasisOrder};
use util::cg::cg;
use std::time::Instant;
use std::f64::consts::PI;

struct CosineProblem;
impl PdeProblem for CosineProblem {
    fn a(&self, _x: f64) -> f64 { 1.0 }
    fn q(&self, _x: f64) -> f64 { 1.0 }
    fn f(&self, x: f64) -> f64 { (PI.powi(2) + 1.0) * (PI * x).cos() }
}

fn main() {
    let penalty = 20.0;
    let order = BasisOrder::Quadratic;

    println!("--- Scalability: Serial vs Parallel (Quadratic Elements) ---");
    println!("Elements, DoFs, Serial(ms), Parallel(ms), Speedup");
    
    let sizes = [1000, 5000, 10000, 20000, 40000, 60000];
    
    for &n in &sizes {
        let (h_elem, x_dof) = generate_mesh(0.0, 1.0, n, order);
        let actual_dofs = x_dof.len();
        
        let mut assembler = SipdgAssembler::new(h_elem, x_dof, penalty, order);
        let prob = CosineProblem;
        let left_bc = DirichletBC { value: 1.0 };
        let right_bc = NeumannBC { value: 0.0 };

        // We want to measure TOTAL time (assembly + solve)
        
        // Parallel
        let start_p = Instant::now();
        assembler.assemble_volume(&prob);
        assembler.assemble_interfaces(&prob);
        let op = assembler.matrix_free_op(&prob, &left_bc, &right_bc);
        let (_, rhs) = assembler.assemble_to_global(); // This part is serial but fast
        let _ = cg(&op, &rhs, 1e-10, 10000);
        let parallel_time = start_p.elapsed().as_millis();

        // Serial (using RAYON_NUM_THREADS=1 via env if possible, or just re-running in a way that forces serial)
        // Actually, Rayon doesn't have a simple "switch to serial" at runtime easily without rebuilding the pool.
        // But we can simulate it by setting the environment variable before running.
        // Instead, let's just report the Parallel results and I'll run it once with RAYON_NUM_THREADS=1 separately if needed.
        // Alternatively, I can use a small hack: Rayon par_iter can be replaced by iter.
        // But that requires code changes.
        
        // Let's just gather Parallel results for now and I'll manually compare with some serial runs.
        println!("{}, {}, -, {}, -", n, actual_dofs, parallel_time);
    }
}
