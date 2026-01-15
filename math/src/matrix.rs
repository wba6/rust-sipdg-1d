use std::{ops::{Add, Index, IndexMut, Mul}, process::Output};

#[derive(Debug, Clone, PartialEq)]
pub struct Matrix<T> {
    rows: usize,
    cols: usize,
    data: Vec<T>
}

impl<T> Matrix<T> {
    // Creates a new matrix initialized to value 
    pub fn new(rows: usize, cols: usize, value: T) -> Self
    where 
        T: Clone,
    {
        let len = rows * cols;
        Self {
            rows,
            cols,
            data: vec![value; len],
        }
    }

    pub fn from_vec(rows: usize, cols: usize, data: Vec<T>) -> Self 
    where 
        T: Clone,
    {
        assert_eq!(rows * cols, data.len(), "Vector is not the correct size (rows * cols)");
        Self { rows, cols, data}
    }

    pub fn rows(&self) -> usize { self.rows }
    pub fn cols(&self) -> usize { self.cols }
    pub fn len(&self) -> usize { self.data.len() }
    pub fn is_empty(&self) -> bool { self.data.is_empty() }

    fn flat_index(&self, r: usize, c: usize) -> usize {
        debug_assert!(r < self.rows && c < self.cols, "Index out of bounds");
        r * self.cols + c
    }

    pub fn get(&self, r: usize, c: usize) -> Option<&T> {
        if r < self.rows && c < self.cols {
            Some(&self.data[r * self.cols + c])
        } else {
            None
        }
    }

    pub fn get_mut(&mut self, r: usize, c: usize) -> Option<&mut T> {
        if r < self.rows && c < self.cols {
            let i = r * self.cols + c;
            Some(&mut self.data[i])
        } else {
            None
        }
    }

}

impl Matrix<f64> {
    /// Convenience: zero matrix for f64.
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self::new(rows, cols, 0.0)
    }
}

impl<T> Index<(usize, usize)> for Matrix<T> {
    type Output = T;

    fn index(&self, (r, c): (usize, usize)) -> &Self::Output {
        let i = self.flat_index(r, c);
        &self.data[i]
    }
}

impl<T> IndexMut<(usize, usize)> for Matrix<T> {
    fn index_mut(&mut self, (r, c): (usize, usize)) -> &mut Self::Output {
        let i = self.flat_index(r, c);
        &mut self.data[i]
    }
}

impl<T> Mul<&T> for &Matrix<T>
where
    T: Clone + Default + Mul<Output = T>,
{
    type Output = Matrix<T>;

    fn mul(self, scalar: &T) -> Self::Output {
        let mut result = Matrix::<T>::new(self.rows, self.cols, T::default());

        for i in 0..self.data.len() {
            result.data[i] = self.data[i].clone() * scalar.clone();
        }

        result
    }
}
