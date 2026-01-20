#[cfg(test)]
pub mod tests {
    use util::matrix::Matrix;

    #[test]
    fn new_fills_with_value_row_major_len() {
        let m = Matrix::new(2, 3, 7i32);
        assert_eq!(m.rows(), 2);
        assert_eq!(m.cols(), 3);
        assert_eq!(m.len(), 6);
        assert!(!m.is_empty());

        // all values are 7
        for r in 0..2 {
            for c in 0..3 {
                assert_eq!(m[(r, c)], 7);
            }
        }
    }

    #[test]
    fn new_zero_by_zero_is_empty() {
        let m: Matrix<i32> = Matrix::new(0, 0, 123);
        assert_eq!(m.rows(), 0);
        assert_eq!(m.cols(), 0);
        assert_eq!(m.len(), 0);
        assert!(m.is_empty());
    }

    #[test]
    fn from_vec_success_and_row_major_mapping() {
        // row-major: [ (0,0)=1, (0,1)=2, (0,2)=3, (1,0)=4, (1,1)=5, (1,2)=6 ]
        let m = Matrix::from_vec(2, 3, vec![1, 2, 3, 4, 5, 6]);
        assert_eq!(m[(0, 0)], 1);
        assert_eq!(m[(0, 1)], 2);
        assert_eq!(m[(0, 2)], 3);
        assert_eq!(m[(1, 0)], 4);
        assert_eq!(m[(1, 1)], 5);
        assert_eq!(m[(1, 2)], 6);
    }

    #[test]
    #[should_panic(expected = "Vector is not the correct size")]
    fn from_vec_panics_on_wrong_length() {
        let _ = Matrix::from_vec(2, 3, vec![1, 2, 3]); // should be 6 entries
    }

    #[test]
    fn get_returns_some_in_bounds_and_none_out_of_bounds() {
        let m = Matrix::from_vec(2, 2, vec![10, 20, 30, 40]);

        assert_eq!(m.get(0, 0), Some(&10));
        assert_eq!(m.get(0, 1), Some(&20));
        assert_eq!(m.get(1, 0), Some(&30));
        assert_eq!(m.get(1, 1), Some(&40));

        // out of bounds
        assert_eq!(m.get(2, 0), None);
        assert_eq!(m.get(0, 2), None);
        assert_eq!(m.get(99, 99), None);
    }

    #[test]
    fn get_mut_allows_mutation_in_bounds_and_none_out_of_bounds() {
        let mut m = Matrix::from_vec(2, 2, vec![1, 2, 3, 4]);

        // mutate (1,0)
        if let Some(x) = m.get_mut(1, 0) {
            *x = 99;
        } else {
            panic!("expected Some for in-bounds get_mut");
        }

        assert_eq!(m[(1, 0)], 99);

        // out of bounds
        assert!(m.get_mut(2, 0).is_none());
        assert!(m.get_mut(0, 2).is_none());
    }


    #[test]
    fn index_and_index_mut_work() {
        let mut m = Matrix::new(2, 2, 0i32);
        m[(0, 0)] = 5;
        m[(0, 1)] = 6;
        m[(1, 0)] = 7;
        m[(1, 1)] = 8;

        assert_eq!(m[(0, 0)], 5);
        assert_eq!(m[(0, 1)], 6);
        assert_eq!(m[(1, 0)], 7);
        assert_eq!(m[(1, 1)], 8);
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "Index out of bounds")]
    fn index_panics_out_of_bounds_in_debug() {
        let m = Matrix::new(1, 1, 0i32);
        let _ = m[(1, 0)];
    }

    #[test]
    fn zeros_makes_all_zeros() {
        let m = Matrix::<f64>::zeros(2, 3);
        for r in 0..2 {
            for c in 0..3 {
                assert_eq!(m[(r, c)], 0.0);
            }
        }
    }

    #[test]
    fn scalar_mul_f64_basic() {
        let m = Matrix::from_vec(2, 2, vec![1.0, -2.0, 3.0, 0.5]);
        let s = 2.0;

        let out = &m * &s;

        assert_eq!(out.rows(), 2);
        assert_eq!(out.cols(), 2);
        assert_eq!(out[(0, 0)], 2.0);
        assert_eq!(out[(0, 1)], -4.0);
        assert_eq!(out[(1, 0)], 6.0);
        assert_eq!(out[(1, 1)], 1.0);
    }

    #[test]
    fn scalar_mul_by_zero_gives_zero_matrix() {
        let m = Matrix::from_vec(2, 3, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let s = 0.0;

        let out = &m * &s;

        for r in 0..2 {
            for c in 0..3 {
                assert_eq!(out[(r, c)], 0.0);
            }
        }
    }

    #[test]
    fn scalar_mul_by_one_returns_same_values() {
        let m = Matrix::from_vec(2, 3, vec![1.0, -2.0, 3.5, 4.0, 0.0, 6.25]);
        let s = 1.0;

        let out = &m * &s;

        assert_eq!(out, m);
    }

    #[test]
    fn scalar_mul_works_for_non_f64_types_with_default() {
        // i32: Default = 0, Mul output is i32
        let m = Matrix::from_vec(2, 2, vec![1i32, 2, 3, 4]);
        let s = 3i32;

        let out = &m * &s;

        assert_eq!(out[(0, 0)], 3);
        assert_eq!(out[(0, 1)], 6);
        assert_eq!(out[(1, 0)], 9);
        assert_eq!(out[(1, 1)], 12);
    }


    #[test]
    fn multiplying_does_not_modify_original() {
        let m = Matrix::from_vec(2, 2, vec![1.0, 2.0, 3.0, 4.0]);
        let s = 10.0;

        let _out = &m * &s;

        // original unchanged
        assert_eq!(m[(0, 0)], 1.0);
        assert_eq!(m[(0, 1)], 2.0);
        assert_eq!(m[(1, 0)], 3.0);
        assert_eq!(m[(1, 1)], 4.0);
    }

}

