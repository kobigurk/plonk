//! A polynomial represented in evaluations form over a domain of size 2^n.

use super::domain::EvaluationDomain;
use super::polynomial::Polynomial;
use core::ops::{Add, AddAssign, DivAssign, Index, Mul, MulAssign, Sub, SubAssign};
use algebra::PrimeField;

/// Stores a polynomial in evaluation form.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Evaluations<F: PrimeField> {
    /// The evaluations of a polynomial over the domain `D`
    pub evals: Vec<F>,
    #[doc(hidden)]
    domain: EvaluationDomain<F>,
}

impl<F: PrimeField> Evaluations<F> {
    /// Construct `Self` from evaluations and a domain.
    pub fn from_vec_and_domain(evals: Vec<F>, domain: EvaluationDomain<F>) -> Self {
        Self { evals, domain }
    }

    /// Interpolate a polynomial from a list of evaluations
    pub fn interpolate_by_ref(&self) -> Polynomial<F> {
        Polynomial::from_coefficients_vec(self.domain.ifft(&self.evals))
    }

    /// Interpolate a polynomial from a list of evaluations
    pub fn interpolate(self) -> Polynomial<F> {
        let Self { mut evals, domain } = self;
        domain.ifft_in_place(&mut evals);
        Polynomial::from_coefficients_vec(evals)
    }
}

impl<F: PrimeField> Index<usize> for Evaluations<F> {
    type Output = F;

    fn index(&self, index: usize) -> &F {
        &self.evals[index]
    }
}

impl<'a, 'b, F: PrimeField> Mul<&'a Evaluations<F>> for &'b Evaluations<F> {
    type Output = Evaluations<F>;

    #[inline]
    fn mul(self, other: &'a Evaluations<F>) -> Evaluations<F> {
        let mut result = self.clone();
        result *= other;
        result
    }
}

impl<'a, F: PrimeField> MulAssign<&'a Evaluations<F>> for Evaluations<F> {
    #[inline]
    fn mul_assign(&mut self, other: &'a Evaluations<F>) {
        assert_eq!(self.domain, other.domain, "domains are unequal");
        self.evals
            .iter_mut()
            .zip(&other.evals)
            .for_each(|(a, b)| *a *= b);
    }
}

impl<'a, 'b, F: PrimeField> Add<&'a Evaluations<F>> for &'b Evaluations<F> {
    type Output = Evaluations<F>;

    #[inline]
    fn add(self, other: &'a Evaluations<F>) -> Evaluations<F> {
        let mut result = self.clone();
        result += other;
        result
    }
}

impl<'a, F: PrimeField> AddAssign<&'a Evaluations<F>> for Evaluations<F> {
    #[inline]
    fn add_assign(&mut self, other: &'a Evaluations<F>) {
        assert_eq!(self.domain, other.domain, "domains are unequal");
        self.evals
            .iter_mut()
            .zip(&other.evals)
            .for_each(|(a, b)| *a += b);
    }
}

impl<'a, 'b, F: PrimeField> Sub<&'a Evaluations<F>> for &'b Evaluations<F> {
    type Output = Evaluations<F>;

    #[inline]
    fn sub(self, other: &'a Evaluations<F>) -> Evaluations<F> {
        let mut result = self.clone();
        result -= other;
        result
    }
}

impl<'a, F: PrimeField> SubAssign<&'a Evaluations<F>> for Evaluations<F> {
    #[inline]
    fn sub_assign(&mut self, other: &'a Evaluations<F>) {
        assert_eq!(self.domain, other.domain, "domains are unequal");
        self.evals
            .iter_mut()
            .zip(&other.evals)
            .for_each(|(a, b)| *a -= b);
    }
}

impl<'a, F: PrimeField> DivAssign<&'a Evaluations<F>> for Evaluations<F> {
    #[inline]
    fn div_assign(&mut self, other: &'a Evaluations<F>) {
        assert_eq!(self.domain, other.domain, "domains are unequal");
        self.evals
            .iter_mut()
            .zip(&other.evals)
            .for_each(|(a, b)| *a *= b.inverse().unwrap());
    }
}
