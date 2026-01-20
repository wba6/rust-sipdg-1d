
/// Generates a vector of `num_elements` linearly (evenly) spaced points between two endpoints.
///
/// # Arguments
///
/// * `start` - The starting value of the sequence (first element of the returned vector).
/// * `end` - The ending value of the sequence (last element of the returned vector).
/// * `num_elements` - The total number of points to generate. Must be greater than 0.
///
/// # Returns
///
/// A `Vec<f64>` containing `num_elements` points linearly spaced from `start` to `end`,
/// inclusive. The first element is exactly `start`, and the last element is forced to be
/// exactly `end` to avoid floating-point accumulation error.
///
/// # Edge cases
///
/// * If `num_elements == 1`, the function returns a single-element vector containing only
///   `start`. In this case, `end` is ignored.
/// * If `start > end`, the function still produces linearly spaced values, but the step
///   between elements is negative.
///
/// # Panics
///
/// Panics if `num_elements == 0`.
///
/// # Examples
///
/// Basic usage:
///
/// 
/// # use math::linespace;
/// let xs = linespace(0.0, 1.0, 5);
/// assert_eq!(xs, vec![0.0, 0.25, 0.5, 0.75, 1.0]);
/// 
///
/// Reverse range:
///
///
/// # use math::linespace;
/// let xs = linespace(5.0, 1.0, 5);
/// assert_eq!(xs, vec![5.0, 4.0, 3.0, 2.0, 1.0]);
/// 
pub fn linespace(start: f64, end: f64, num_elements: usize) -> Vec<f64> {
    assert!(num_elements > 0, "num_elements must be greater than 0");

    // Case for exactly one point (avoids division by zero)
    if num_elements == 1 {
        return vec![start];
    }

    let mut result: Vec<f64> = Vec::with_capacity(num_elements);
    
    let step = (end - start) / ((num_elements - 1) as f64);

    for i in 0..num_elements {
        result.push(start + (i as f64 * step));
    }
    
    if let Some(last) = result.last_mut() {
        *last = end;
    }

    return result;
}

