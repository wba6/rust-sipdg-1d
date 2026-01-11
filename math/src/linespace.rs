
/// generates a row vector of a specified number of linearly (evenly) spaced points between two endpoints
pub fn linespace(start: f64, end: f64, num_elements: usize) -> Vec<f64> {
    assert!(num_elements > 0);
    // Create a vector or num_elements + 1 due to this being inclusive of the domain
    let mut result: Vec<f64> = vec![0.0; num_elements + 1];
    let spacing: f64 = (end - start).abs()/num_elements as f64;
    for (index, value) in result.iter_mut().enumerate() {
        *value = start + index as f64 * spacing;
    }
    return result;
}

