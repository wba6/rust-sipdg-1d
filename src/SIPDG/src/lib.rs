// src/SIPDG/src/lib.rs

pub mod mesh;
pub mod pde;
pub mod problem_loader;
pub mod plotter;

// Re-export commonly used items for convenience
pub use pde::{PdeProblem, SipdgAssembler, DirichletBC, NeumannBC, BoundaryCondition, Side};
pub use mesh::generate_mesh;
