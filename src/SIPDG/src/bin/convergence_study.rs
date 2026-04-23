use sipdg::{PdeProblem, SipdgAssembler, DirichletBC, generate_mesh, pde::BasisOrder};
use util::cg::cg;
use util::quadrature::get_gauss_legendre_3pts;
use std::f64::consts::PI;
use plotters::prelude::*;

struct SineProblem;

impl PdeProblem for SineProblem {
    fn a(&self, _x: f64) -> f64 { 1.0 }
    fn q(&self, _x: f64) -> f64 { 1.0 }
    fn f(&self, x: f64) -> f64 { (PI.powi(2) + 1.0) * (PI * x).sin() }
}

fn exact_soln(x: f64) -> f64 { (PI * x).sin() }

fn run_study(order: BasisOrder, sizes: &[usize], penalty: f64) -> Vec<(usize, f64)> {
    let mut results = Vec::new();

    for &num_elements in sizes {
        let (h_elem, x_dof) = generate_mesh(0.0, 1.0, num_elements, order);
        let mut assembler = SipdgAssembler::new(h_elem, x_dof, penalty, order);
        let prob = SineProblem;

        assembler.assemble_volume(&prob);
        assembler.assemble_interfaces(&prob);

        let bc = DirichletBC { value: 0.0 };
        let mut global_rhs = assembler.assemble_rhs();
        assembler.apply_boundaries_rhs(&prob, &bc, &bc, &mut global_rhs);

        let op = assembler.matrix_free_op(&prob, &bc, &bc);
        let p = cg(&op, &global_rhs, 1e-12, 20000);

        let quad_pts = get_gauss_legendre_3pts();
        let n_nodes = order.num_nodes();

        let l2_err_sq: f64 = assembler.elements.iter().map(|elem| {
            let det_j = elem.jacobian();
            let mut local_err_sq = 0.0;
            for pt in &quad_pts {
                let x = elem.map_to_physical(pt.xi);
                let mut u_h = 0.0;
                for i in 0..n_nodes {
                    u_h += p[elem.nodes[i]] * elem.phi(i, pt.xi);
                }
                let u_ex = exact_soln(x);
                let err = u_h - u_ex;
                local_err_sq += pt.weight * err * err * det_j;
            }
            local_err_sq
        }).sum();

        let l2_error = l2_err_sq.sqrt();
        results.push((num_elements, l2_error));
        println!("N = {}, L2 Error = {:.6e}", num_elements, l2_error);
    }
    results
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sizes = [10, 20, 40, 80, 160, 320];
    let penalty = 20.0;

    println!("--- Linear Elements (p=1) ---");
    let results_p1 = run_study(BasisOrder::Linear, &sizes, penalty);
    
    println!("\nN, Error, Rate");
    for i in 0..results_p1.len() {
        if i == 0 {
            println!("{}, {:.6e}, -", results_p1[i].0, results_p1[i].1);
        } else {
            let rate = (results_p1[i-1].1 / results_p1[i].1).log2();
            println!("{}, {:.6e}, {:.2}", results_p1[i].0, results_p1[i].1, rate);
        }
    }

    println!("\n--- Quadratic Elements (p=2) ---");
    let results_p2 = run_study(BasisOrder::Quadratic, &sizes, penalty);

    println!("\nN, Error, Rate");
    for i in 0..results_p2.len() {
        if i == 0 {
            println!("{}, {:.6e}, -", results_p2[i].0, results_p2[i].1);
        } else {
            let rate = (results_p2[i-1].1 / results_p2[i].1).log2();
            println!("{}, {:.6e}, {:.2}", results_p2[i].0, results_p2[i].1, rate);
        }
    }

    // --- Plotting ---
    let root = BitMapBackend::new("convergence_plot.png", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let min_n = sizes[0] as f64;
    let max_n = sizes[sizes.len() - 1] as f64;
    let min_err = results_p2[results_p2.len() - 1].1 * 0.5;
    let max_err = results_p1[0].1 * 2.0;

    let mut chart = ChartBuilder::on(&root)
        .caption("SIPDG Convergence Study", ("sans-serif", 30).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(
            (min_n..max_n).log_scale(),
            (min_err..max_err).log_scale(),
        )?;

    chart.configure_mesh()
        .x_desc("Number of Elements (N)")
        .y_desc("L2 Error Norm")
        .draw()?;

    chart.draw_series(LineSeries::new(
        results_p1.iter().map(|&(n, err)| (n as f64, err)),
        &BLUE,
    ))?
    .label("Linear (p=1)")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart.draw_series(LineSeries::new(
        results_p2.iter().map(|&(n, err)| (n as f64, err)),
        &RED,
    ))?
    .label("Quadratic (p=2)")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // Reference slopes
    let n1 = results_p1[0].0 as f64;
    let e1 = results_p1[0].1;
    chart.draw_series(LineSeries::new(
        sizes.iter().map(|&n| (n as f64, e1 * (n1 / n as f64).powi(2))),
        BLUE.mix(0.3),
    ))?
    .label("O(h^2) Slope")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE.mix(0.3)));

    let n2 = results_p2[0].0 as f64;
    let e2 = results_p2[0].1;
    chart.draw_series(LineSeries::new(
        sizes.iter().map(|&n| (n as f64, e2 * (n2 / n as f64).powi(3))),
        RED.mix(0.3),
    ))?
    .label("O(h^3) Slope")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED.mix(0.3)));

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;
    println!("\nConvergence plot saved to convergence_plot.png");

    Ok(())
}
