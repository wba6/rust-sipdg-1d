use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use crate::pde::ConfigurableProblem;

pub fn load_problem_from_file(path: &Path) -> ConfigurableProblem {
    let mut prob = ConfigurableProblem::default();

    // Open the file or panic with a descriptive message
    let file = File::open(path).expect("Unable to open config file");
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line.unwrap_or_default();
        if line.is_empty() || line.starts_with('#') { continue; } // Skip empty/comments

        let parts: Vec<&str> = line.split('=').map(|s| s.trim()).collect();
        
        if parts.len() == 2 {
            if let Ok(val) = parts[1].parse::<f64>() {
                match parts[0] {
                    "p" => prob.p_val = val,
                    "q" => prob.q_val = val,
                    "f" => prob.f_val = val,
                    _ => println!("Warning: Unknown parameter {}", parts[0]),
                }
            }
        }
    }
    
    prob
}
