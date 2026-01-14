#[cfg(test)]
pub mod tests {
    use math::matrix::Matrix;

    #[test]
    fn new_creates_correct_dimensions() {
        let m = Matrix::new(3, 5);
        assert_eq!(m.data.len(), 3);
        assert_eq!(m.data[0].len(), 5);
        assert_eq!(m.data[1].len(), 5);
        assert_eq!(m.data[2].len(), 5);
    }

    #[test]
    fn new_initializes_all_entries_to_zero() {
        let m = Matrix::new(4, 6);
        for row in 0..4 {
            for col in 0..6 {
                assert_eq!(m[row][col], 0.0);
            }
        }
    }

    #[test]
    fn index_returns_row_reference() {
        let m = Matrix::new(2, 3);
        // Index gives you a &Vec<f64> (row)
        assert_eq!(m[0].len(), 3);
        assert_eq!(m[1].len(), 3);
    }

    #[test]
    fn index_mut_allows_writing_values() {
        let mut m = Matrix::new(2, 2);

        m[0][0] = 1.0;
        m[0][1] = 2.5;
        m[1][0] = -3.0;
        m[1][1] = 4.25;

        assert_eq!(m[0][0], 1.0);
        assert_eq!(m[0][1], 2.5);
        assert_eq!(m[1][0], -3.0);
        assert_eq!(m[1][1], 4.25);
    }

    #[test]
    fn clone_and_partial_eq_work() {
        let mut a = Matrix::new(2, 3);
        a[1][2] = 9.0;

        let b = a.clone();
        assert_eq!(a, b);

        let mut c = b.clone();
        c[0][0] = 1.0;
        assert_ne!(a, c);
    }

    #[test]
    fn rows_are_independent_not_aliased() {
        let mut m = Matrix::new(3, 3);
        m[0][0] = 7.0;

        // If rows were accidentally aliased, this would also change.
        assert_eq!(m[1][0], 0.0);
        assert_eq!(m[2][0], 0.0);
    }

    #[test]
    #[should_panic]
    fn index_out_of_bounds_panics() {
        let m = Matrix::new(2, 2);
        let _ = &m[2]; // valid indices: 0,1
    }

    #[test]
    #[should_panic]
    fn inner_index_out_of_bounds_panics() {
        let m = Matrix::new(2, 2);
        let _ = m[0][2]; // valid indices: 0,1
    }

    #[test]
    fn zero_sized_matrix_is_allowed() {
        let m0 = Matrix::new(0, 5);
        assert_eq!(m0.data.len(), 0);

        let m1 = Matrix::new(3, 0);
        assert_eq!(m1.data.len(), 3);
        assert_eq!(m1.data[0].len(), 0);
        assert_eq!(m1.data[1].len(), 0);
        assert_eq!(m1.data[2].len(), 0);
    }

    #[test]
    fn scalar_mul_scales_all_entries() {
        let mut m = Matrix::new(2, 3);
        m[0][0] = 1.0;
        m[0][1] = -2.0;
        m[0][2] = 0.5;
        m[1][0] = 3.0;
        m[1][1] = 4.0;
        m[1][2] = -1.5;

        let r = m.clone() * 2.0;

        assert_eq!(r[0][0], 2.0);
        assert_eq!(r[0][1], -4.0);
        assert_eq!(r[0][2], 1.0);
        assert_eq!(r[1][0], 6.0);
        assert_eq!(r[1][1], 8.0);
        assert_eq!(r[1][2], -3.0);

        // original unchanged (since we used clone())
        assert_eq!(m[0][0], 1.0);
        assert_eq!(m[1][2], -1.5);
    }

   #[test]
   fn scalar_mul_by_zero_gives_zero_matrix() {
       let mut m = Matrix::new(2, 2);
       m[0][0] = 1.0;
       m[0][1] = 2.0;
       m[1][0] = -3.0;
       m[1][1] = 4.0;

       let r = m * 0.0;

       assert_eq!(r[0][0], 0.0);
       assert_eq!(r[0][1], 0.0);
       assert_eq!(r[1][0], 0.0);
       assert_eq!(r[1][1], 0.0);
   }

   #[test]
   fn scalar_mul_preserves_dimensions() {
       let m = Matrix::new(4, 7);
       let r = m * 3.0;

       assert_eq!(r[0].len(), 7);
       assert_eq!(r[3].len(), 7);
       assert_eq!(r.data.len(), 4);
   }

   #[test]
   fn scalar_mul_by_zero_gives_zero_matrix_rhs_matrix() {
       let mut m = Matrix::new(2, 2);
       m[0][0] = 1.0;
       m[0][1] = 2.0;
       m[1][0] = -3.0;
       m[1][1] = 4.0;

       let r = 0.0 * m;

       assert_eq!(r[0][0], 0.0);
       assert_eq!(r[0][1], 0.0);
       assert_eq!(r[1][0], 0.0);
       assert_eq!(r[1][1], 0.0);
   }

   #[test]
   fn scalar_mul_preserves_dimensions_rhs_matrix() {
       let m = Matrix::new(4, 7);
       let r = 3.0 * m;

       assert_eq!(r[0].len(), 7);
       assert_eq!(r[3].len(), 7);
       assert_eq!(r.data.len(), 4);
   }

}

