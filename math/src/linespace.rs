
/// generates a row vector of a specified number of linearly (evenly) spaced points between two endpoints
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

