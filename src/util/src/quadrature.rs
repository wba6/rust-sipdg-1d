pub struct QuadraturePoint {
    pub xi: f64,
    pub weight: f64,
}

pub fn get_gauss_legendre_2pts() -> Vec<QuadraturePoint> {
    let xi = 1.0 / 3.0f64.sqrt();
    vec![
        QuadraturePoint { xi: -xi, weight: 1.0 },
        QuadraturePoint { xi: xi, weight: 1.0 },
    ]
}

pub fn get_gauss_legendre_3pts() -> Vec<QuadraturePoint> {
    let xi = 0.6f64.sqrt();
    vec![
        QuadraturePoint { xi: -xi, weight: 5.0 / 9.0 },
        QuadraturePoint { xi: 0.0, weight: 8.0 / 9.0 },
        QuadraturePoint { xi: xi, weight: 5.0 / 9.0 },
    ]
}

pub fn integrate<F>(f: F, domain_a: f64, domain_b: f64, points: &[QuadraturePoint]) -> f64
where
    F: Fn(f64) -> f64,
{
    let h = (domain_b - domain_a) / 2.0;
    let mid = (domain_a + domain_b) / 2.0;
    let mut sum = 0.0;
    for pt in points {
        let x = mid + h * pt.xi;
        sum += pt.weight * f(x);
    }
    sum * h
}
