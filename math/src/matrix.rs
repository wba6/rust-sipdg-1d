use std::ops::{Index, IndexMut, Mul};

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

impl Mul<f64> for Matrix {
    type Output = Matrix;

    fn mul(self, scalar: f64) -> Self::Output {
        let mut result = Matrix::new(self.data.len(), self.data[0].len());
        for (index, row) in self.data.iter().enumerate() {
            result.data[index] = row.iter().map(|&x| x*scalar).collect();
        }
        result
    }
}

impl Mul<Matrix> for f64 {
    type Output = Matrix;

    fn mul(self, matrix: Matrix) -> Self::Output {
        let mut result = Matrix::new(matrix.data.len(), matrix.data[0].len());
        for (index, row) in matrix.data.iter().enumerate() {
            result.data[index] = row.iter().map(|&x| x*self).collect();
        }
        result
    }
}

