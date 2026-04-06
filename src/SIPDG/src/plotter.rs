use plotters::prelude::*;
use crate::pde::BasisOrder;

fn evaluate_uh(xi: f64, local_u: &[f64], order: BasisOrder) -> f64 {
    match order {
        BasisOrder::Linear => {
            let phi0 = 0.5 * (1.0 - xi);
            let phi1 = 0.5 * (1.0 + xi);
            local_u[0] * phi0 + local_u[1] * phi1
        }
        BasisOrder::Quadratic => {
            let phi0 = 0.5 * xi * (xi - 1.0);
            let phi1 = 1.0 - xi * xi;
            let phi2 = 0.5 * xi * (xi + 1.0);
            local_u[0] * phi0 + local_u[1] * phi1 + local_u[2] * phi2
        }
    }
}

/// Plots the SIPDG numerical solution against the exact solution.
///
/// # Arguments
/// * `x_dof` - Global coordinates of the degrees of freedom.
/// * `u` - Numerical solution vector.
/// * `soln_function` - Closure representing the exact solution u(x).
/// * `num_elements` - The number of elements in the mesh.
/// * `order` - The basis order used.
pub fn plot_results(
    x_dof: &[f64],
    u: &[f64],
    soln_function: impl Fn(f64) -> f64,
    num_elements: usize,
    order: BasisOrder,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new("solution_plot.png", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    // Calculate scaling limits based on both exact and numerical data
    let exact_max = (0..100)
        .map(|i| soln_function(i as f64 / 100.0))
        .fold(0.0, f64::max);
    
    let numerical_max = u.iter()
        .filter(|&&val| !val.is_nan() && !val.is_infinite())
        .fold(0.0, |a: f64, &b| a.max(b));

    // Use a buffer for the Y-axis height
    let y_limit = exact_max.max(numerical_max) * 1.2;

    let mut chart = ChartBuilder::on(&root)
        .caption("1D SIPDG Solver Results", ("sans-serif", 30).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(0.0..1.0, 0.0..y_limit)?;

    chart.configure_mesh()
        .x_desc("x")
        .y_desc("u(x)")
        .draw()?;

    // Plot Exact Solution as a smooth continuous line
    chart.draw_series(LineSeries::new(
        (0..100).map(|i| {
            let x = i as f64 / 100.0;
            (x, soln_function(x))
        }),
        &RED,
    ))?
    .label("Exact Solution")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // Plot SIPDG Numerical Solution
    // In DG, we plot each element as a separate segment because nodes are not shared.
    let n_nodes = order.num_nodes();
    let n_samples = 10;
    for i in 0..num_elements {
        let start_idx = i * n_nodes;
        let local_u = &u[start_idx..start_idx + n_nodes];
        let x_l = x_dof[start_idx];
        let x_r = x_dof[start_idx + n_nodes - 1];

        chart.draw_series(LineSeries::new(
            (0..n_samples).map(|j| {
                let xi = -1.0 + 2.0 * (j as f64) / ((n_samples - 1) as f64);
                let x = x_l + (xi + 1.0) * (x_r - x_l) / 2.0;
                (x, evaluate_uh(xi, local_u, order))
            }),
            &BLUE,
        ))?;
    }

    // Add a dummy entry to the legend for the Blue SIPDG lines
    chart.draw_series(LineSeries::new(vec![], &BLUE))?
        .label("SIPDG (Numerical)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    println!("Plot successfully saved to solution_plot.png");
    Ok(())
}
