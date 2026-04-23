// src/util/tests/quadrature_tests.rs
// We check if the Gauss-Legendre rules can integrate polynomials exactly.

#[cfg(test)]
mod tests {
    use util::quadrature::{get_gauss_legendre_2pts, get_gauss_legendre_3pts, integrate};

    /// A small tolerance for comparing floating point results.
    const TOL: f64 = 1e-12;

    #[test]
    fn test_gauss_2pt_constant() {
        // Integral of f(x) = 5 from 0 to 2 is 10.
        let pts = get_gauss_legendre_2pts();
        let result = integrate(|_x| 5.0, 0.0, 2.0, &pts);
        assert!((result - 10.0).abs() < TOL, "2-point rule should be exact for constants. Got {}", result);
    }

    #[test]
    fn test_gauss_2pt_linear() {
        // Integral of f(x) = x from 0 to 1 is 0.5.
        let pts = get_gauss_legendre_2pts();
        let result = integrate(|x| x, 0.0, 1.0, &pts);
        assert!((result - 0.5).abs() < TOL, "2-point rule should be exact for linear functions. Got {}", result);
    }

    #[test]
    fn test_gauss_2pt_quadratic() {
        // Integral of f(x) = x^2 from 0 to 1 is 1/3.
        let pts = get_gauss_legendre_2pts();
        let result = integrate(|x| x * x, 0.0, 1.0, &pts);
        assert!((result - 1.0/3.0).abs() < TOL, "2-point rule should be exact for quadratics. Got {}", result);
    }

    #[test]
    fn test_gauss_2pt_cubic() {
        // Integral of f(x) = x^3 from 0 to 1 is 1/4.
        // Gauss-Legendre with N points is exact for polynomials of degree 2N-1.
        // For N=2, it should be exact for degree 3.
        let pts = get_gauss_legendre_2pts();
        let result = integrate(|x| x * x * x, 0.0, 1.0, &pts);
        assert!((result - 0.25).abs() < TOL, "2-point rule should be exact for cubics. Got {}", result);
    }

    #[test]
    fn test_gauss_3pt_quartic() {
        // Integral of f(x) = x^4 from 0 to 1 is 1/5 = 0.2.
        // For N=3, it should be exact for degree 2*3 - 1 = 5.
        let pts = get_gauss_legendre_3pts();
        let result = integrate(|x| x.powi(4), 0.0, 1.0, &pts);
        assert!((result - 0.2).abs() < TOL, "3-point rule should be exact for quartics. Got {}", result);
    }

    #[test]
    fn test_gauss_3pt_quintic() {
        // Integral of f(x) = x^5 from 0 to 1 is 1/6.
        let pts = get_gauss_legendre_3pts();
        let result = integrate(|x| x.powi(5), 0.0, 1.0, &pts);
        assert!((result - 1.0/6.0).abs() < TOL, "3-point rule should be exact for quintics. Got {}", result);
    }
}
