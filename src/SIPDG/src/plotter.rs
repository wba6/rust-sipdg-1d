use plotters::prelude::*;

/// Plots the SIPDG numerical solution against the exact solution.
///
/// # Arguments
/// * `x_dof` - Global coordinates of the degrees of freedom.
/// * `u` - Numerical solution vector.
/// * `soln_function` - Closure representing the exact solution u(x).
/// * `num_elements` - The number of elements in the mesh.
pub fn plot_results(
    x_dof: &[f64],
    u: &[f64],
    soln_function: impl Fn(f64) -> f64,
    num_elements: usize,
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

    // 1. Plot Exact Solution as a smooth continuous line
    chart.draw_series(LineSeries::new(
        (0..100).map(|i| {
            let x = i as f64 / 100.0;
            (x, soln_function(x))
        }),
        &RED,
    ))?
    .label("Exact Solution")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // 2. Plot SIPDG Numerical Solution
    // In DG, we plot each element as a separate segment because nodes are not shared.
    for i in 0..num_elements {
        let idx_l = 2 * i;
        let idx_r = 2 * i + 1;

        // Safety check for n_dof bounds (len is 2 * num_elements)
        if idx_r < u.len() {
            chart.draw_series(LineSeries::new(
                vec![(x_dof[idx_l], u[idx_l]), (x_dof[idx_r], u[idx_r])],
                &BLUE,
            ))?;
        }
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
