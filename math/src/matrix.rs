use std::ops::{Index, IndexMut};

#[derive(Debug, Clone, PartialEq)]
pub struct Matrix {
    pub data: Vec<Vec<f64>>
}

impl Matrix {
    // Creates a new matrix initialized to zero
    pub fn new(x_size: usize, y_size: usize) -> Self {
       Matrix {
           data: vec![vec![0.0; y_size]; x_size]
       }
    }
}

impl Index<usize> for Matrix {
    type Output = Vec<f64>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl IndexMut<usize> for Matrix {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

