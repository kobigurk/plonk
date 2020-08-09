//! This module contains an implementation of a polynomial in coefficient form
//! Where each coefficient is represented using a position in the underlying vector.
use super::{EvaluationDomain, Evaluations};
use crate::util;
use rand::Rng;
use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelIterator,
};

use std::ops::{Add, AddAssign, Deref, DerefMut, Mul, Neg, Sub, SubAssign};

#[derive(Debug, Eq, PartialEq, Clone)]
/// Polynomial represents a polynomial in coeffiient form.
pub struct Polynomial<F: PrimeField> {
    /// The coefficient of `x^i` is stored at location `i` in `self.coeffs`.
    pub coeffs: Vec<F>,
}

impl<F: PrimeField> Deref for Polynomial<F> {
    type Target = [F];

    fn deref(&self) -> &[F] {
        &self.coeffs
    }
}

impl<F: PrimeField> DerefMut for Polynomial<F> {
    fn deref_mut(&mut self) -> &mut [F] {
        &mut self.coeffs
    }
}

impl<F: PrimeField> Polynomial<F> {
    /// Returns the zero polynomial.
    pub fn zero() -> Self {
        Self { coeffs: Vec::new() }
    }

    /// Checks if the given polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty() || self.coeffs.par_iter().all(|coeff| coeff == &F::zero())
    }

    /// Constructs a new polynomial from a list of coefficients.
    pub fn from_coefficients_slice(coeffs: &[F]) -> Self {
        Self::from_coefficients_vec(coeffs.to_vec())
    }

    /// Constructs a new polynomial from a list of coefficients.
    ///
    /// # Panics
    /// When the length of the coeffs is zero.
    pub fn from_coefficients_vec(coeffs: Vec<F>) -> Self {
        let mut result = Self { coeffs };
        // While there are zeros at the end of the coefficient vector, pop them off.
        result.truncate_leading_zeros();
        // Check that either the coefficients vec is empty or that the last coeff is
        // non-zero.
        assert!(result
            .coeffs
            .last()
            .map_or(true, |coeff| coeff != &F::zero()));

        result
    }

    /// Returns the degree of the polynomial.
    pub fn degree(&self) -> usize {
        if self.is_zero() {
            return 0;
        }
        assert!(self
            .coeffs
            .last()
            .map_or(false, |coeff| coeff != &F::zero()));
        self.coeffs.len() - 1
    }

    fn truncate_leading_zeros(&mut self) {
        while self.coeffs.last().map_or(false, |c| c == &F::zero()) {
            self.coeffs.pop();
        }
    }
    /// Evaluates `self` at the given `point` in the field.
    pub fn evaluate(&self, point: &F) -> F {
        if self.is_zero() {
            return F::zero();
        }

        // Compute powers of points
        let powers = util::powers_of(point, self.len());

        let p_evals: Vec<_> = self
            .par_iter()
            .zip(powers.into_par_iter())
            .map(|(c, p)| p * c)
            .collect();
        let mut sum = F::zero();
        for eval in p_evals.into_iter() {
            sum += &eval;
        }
        sum
    }

    /// Outputs a polynomial of degree `d` where each coefficient is sampled
    /// uniformly at random from the field `F`.
    pub fn rand<R: Rng>(d: usize, mut rng: &mut R) -> Self {
        let mut random_coeffs = Vec::with_capacity(d + 1);
        for _ in 0..=d {
            random_coeffs.push(util::random_scalar(&mut rng));
        }
        Self::from_coefficients_vec(random_coeffs)
    }
}

use std::iter::Sum;
use algebra::PrimeField;

impl<F: PrimeField> Sum for Polynomial<F> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        let sum: Polynomial<F> = iter.fold(Polynomial::zero(), |mut res, val| {
            res = &res + &val;
            res
        });
        sum
    }
}

impl<'a, 'b, F: PrimeField> Add<&'a Polynomial<F>> for &'b Polynomial<F> {
    type Output = Polynomial<F>;

    fn add(self, other: &'a Polynomial<F>) -> Polynomial<F> {
        let mut result = if self.is_zero() {
            other.clone()
        } else if other.is_zero() {
            self.clone()
        } else if self.degree() >= other.degree() {
            let mut result = self.clone();
            for (a, b) in result.coeffs.iter_mut().zip(&other.coeffs) {
                *a += b
            }
            result
        } else {
            let mut result = other.clone();
            for (a, b) in result.coeffs.iter_mut().zip(&self.coeffs) {
                *a += b
            }
            result
        };
        result.truncate_leading_zeros();
        result
    }
}

#[allow(dead_code)]
// Addition method tht uses iterators to add polynomials
// Benchmark this against the original method
fn iter_add<F: PrimeField>(poly_a: &Polynomial<F>, poly_b: &Polynomial<F>) -> Polynomial<F> {
    if poly_a.len() == 0 {
        return poly_b.clone();
    }
    if poly_b.len() == 0 {
        return poly_a.clone();
    }

    let max_len = std::cmp::max(poly_a.len(), poly_b.len());
    let min_len = std::cmp::min(poly_a.len(), poly_b.len());
    let mut data = Vec::with_capacity(max_len);
    let (mut poly_a_iter, mut poly_b_iter) = (poly_a.iter(), poly_b.iter());

    let partial_addition = poly_a_iter
        .by_ref()
        .zip(poly_b_iter.by_ref())
        .map(|(&a, &b)| a + b)
        .take(min_len);

    data.extend(partial_addition);
    data.extend(poly_a_iter);
    data.extend(poly_b_iter);

    assert_eq!(data.len(), std::cmp::max(poly_a.len(), poly_b.len()));

    let mut result = Polynomial::from_coefficients_vec(data);
    result.truncate_leading_zeros();
    result
}

impl<'a, 'b, F: PrimeField> AddAssign<&'a Polynomial<F>> for Polynomial<F> {
    fn add_assign(&mut self, other: &'a Polynomial<F>) {
        if self.is_zero() {
            self.coeffs.truncate(0);
            self.coeffs.extend_from_slice(&other.coeffs);
        } else if other.is_zero() {
        } else if self.degree() >= other.degree() {
            for (a, b) in self.coeffs.iter_mut().zip(&other.coeffs) {
                *a += b
            }
        } else {
            // Add the necessary number of zero coefficients.
            self.coeffs.resize(other.coeffs.len(), F::zero());
            for (a, b) in self.coeffs.iter_mut().zip(&other.coeffs) {
                *a += b
            }
            self.truncate_leading_zeros();
        }
    }
}

impl<'a, 'b, F: PrimeField> AddAssign<(F, &'a Polynomial<F>)> for Polynomial<F> {
    fn add_assign(&mut self, (f, other): (F, &'a Polynomial<F>)) {
        if self.is_zero() {
            self.coeffs.truncate(0);
            self.coeffs.extend_from_slice(&other.coeffs);
            self.coeffs.iter_mut().for_each(|c| *c *= &f);
        } else if other.is_zero() {
        } else if self.degree() >= other.degree() {
            for (a, b) in self.coeffs.iter_mut().zip(&other.coeffs) {
                *a += &(f * b);
            }
        } else {
            // Add the necessary number of zero coefficients.
            self.coeffs.resize(other.coeffs.len(), F::zero());
            for (a, b) in self.coeffs.iter_mut().zip(&other.coeffs) {
                *a += &(f * b);
            }
            self.truncate_leading_zeros();
        }
    }
}

impl<F: PrimeField> Neg for Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn neg(mut self) -> Polynomial<F> {
        for coeff in &mut self.coeffs {
            *coeff = -*coeff;
        }
        self
    }
}

impl<'a, 'b, F: PrimeField> Sub<&'a Polynomial<F>> for &'b Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn sub(self, other: &'a Polynomial<F>) -> Polynomial<F> {
        let mut result = if self.is_zero() {
            let mut result = other.clone();
            for coeff in &mut result.coeffs {
                *coeff = -(*coeff);
            }
            result
        } else if other.is_zero() {
            self.clone()
        } else if self.degree() >= other.degree() {
            let mut result = self.clone();
            for (a, b) in result.coeffs.iter_mut().zip(&other.coeffs) {
                *a -= b
            }
            result
        } else {
            let mut result = self.clone();
            result.coeffs.resize(other.coeffs.len(), F::zero());
            for (a, b) in result.coeffs.iter_mut().zip(&other.coeffs) {
                *a -= b;
            }
            result
        };
        result.truncate_leading_zeros();
        result
    }
}

impl<'a, 'b, F: PrimeField> SubAssign<&'a Polynomial<F>> for Polynomial<F> {
    #[inline]
    fn sub_assign(&mut self, other: &'a Polynomial<F>) {
        if self.is_zero() {
            self.coeffs.resize(other.coeffs.len(), F::zero());
            for (i, coeff) in other.coeffs.iter().enumerate() {
                self.coeffs[i] -= coeff;
            }
        } else if other.is_zero() {
        } else if self.degree() >= other.degree() {
            for (a, b) in self.coeffs.iter_mut().zip(&other.coeffs) {
                *a -= b
            }
        } else {
            // Add the necessary number of zero coefficients.
            self.coeffs.resize(other.coeffs.len(), F::zero());
            for (a, b) in self.coeffs.iter_mut().zip(&other.coeffs) {
                *a -= b
            }
            // If the leading coefficient ends up being zero, pop it off.
            self.truncate_leading_zeros();
        }
    }
}

impl<F: PrimeField> Polynomial<F> {
    #[allow(dead_code)]
    #[inline]
    fn leading_coefficient(&self) -> Option<&F> {
        self.last()
    }

    #[allow(dead_code)]
    #[inline]
    fn iter_with_index(&self) -> Vec<(usize, F)> {
        self.iter().cloned().enumerate().collect()
    }

    /// Divides `self` by x-z using Ruffinis method
    pub fn ruffini(&self, z: F) -> Polynomial<F> {
        let mut quotient: Vec<F> = Vec::with_capacity(self.degree());
        let mut k = F::zero();

        // Reverse the results and use Ruffini's method to compute the quotient
        // The coefficients must be reversed as Ruffini's method
        // starts with the leading coefficient, while Polynomials
        // are stored in increasing order i.e. the leading coefficient is the last element
        for coeff in self.coeffs.iter().rev() {
            let t = *coeff + &k;
            quotient.push(t);
            k = z * &t;
        }

        // Pop off the last element, it is the remainder term
        // For PLONK, we only care about perfect factors
        quotient.pop();

        // Reverse the results for storage in the Polynomial struct
        quotient.reverse();
        Polynomial::from_coefficients_vec(quotient)
    }
}

/// Performs O(nlogn) multiplication of polynomials if F is smooth.
impl<'a, 'b, F: PrimeField> Mul<&'a Polynomial<F>> for &'b Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn mul(self, other: &'a Polynomial<F>) -> Polynomial<F> {
        if self.is_zero() || other.is_zero() {
            Polynomial::zero()
        } else {
            let domain = EvaluationDomain::new(self.coeffs.len() + other.coeffs.len())
                .expect("field is not smooth enough to construct domain");
            let mut self_evals = Evaluations::from_vec_and_domain(domain.fft(&self.coeffs), domain);
            let other_evals = Evaluations::from_vec_and_domain(domain.fft(&other.coeffs), domain);
            self_evals *= &other_evals;
            self_evals.interpolate()
        }
    }
}
/// Convenience Trait to multiply a scalar and polynomial
impl<'a, 'b, F: PrimeField> Mul<&'a F> for &'b Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn mul(self, constant: &'a F) -> Polynomial<F> {
        if self.is_zero() || (constant == &F::zero()) {
            return Polynomial::zero();
        }
        let scaled_coeffs: Vec<_> = self
            .coeffs
            .par_iter()
            .map(|coeff| *coeff * constant)
            .collect();
        Polynomial::from_coefficients_vec(scaled_coeffs)
    }
}
/// Convenience Trait to sub a scalar and polynomial
impl<'a, 'b, F: PrimeField> Add<&'a F> for &'b Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn add(self, constant: &'a F) -> Polynomial<F> {
        if self.is_zero() {
            return Polynomial::from_coefficients_vec(vec![*constant]);
        }
        let mut result = self.clone();
        if constant == &F::zero() {
            return result;
        }

        result[0] += constant;
        result
    }
}
impl<'a, 'b, F: PrimeField> Sub<&'a F> for &'b Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn sub(self, constant: &'a F) -> Polynomial<F> {
        let negated_constant = constant.neg();
        self + &negated_constant
    }
}
#[cfg(test)]
mod test {
    use super::*;
    use algebra::{bls12_381::Fr, Zero, One};
    #[test]
    fn test_ruffini() {
        // X^2 + 4X + 4
        let quadratic = Polynomial::<Fr>::from_coefficients_vec(vec![
            <Fr as From<u64>>::from(4),
            <Fr as From<u64>>::from(4),
            Fr::one(),
        ]);
        // Divides X^2 + 4X + 4 by X+2
        let quotient = quadratic.ruffini(-<Fr as From<u64>>::from(2));
        // X+2
        let expected_quotient =
            Polynomial::from_coefficients_vec(vec![<Fr as From<u64>>::from(2), Fr::one()]);
        assert_eq!(quotient, expected_quotient);
    }
    #[test]
    fn test_ruffini_zero() {
        // Tests the two situations where zero can be added to Ruffini:
        // (1) Zero polynomial in the divided
        // (2) Zero as the constant term for the polynomial you are dividing by
        // In the first case, we should get zero as the quotient
        // In the second case, this is the same as dividing by X
        // (1)
        //
        // Zero polynomial
        let zero = Polynomial::zero();
        // Quotient is invariant under any argument we pass
        let quotient = zero.ruffini(-<Fr as From<u64>>::from(2));
        assert_eq!(quotient, Polynomial::zero());
        // (2)
        //
        // X^2 + X
        let p =
            Polynomial::from_coefficients_vec(vec![Fr::zero(), Fr::one(), Fr::one()]);
        // Divides X^2 + X by X
        let quotient = p.ruffini(Fr::zero());
        // X + 1
        let expected_quotient =
            Polynomial::from_coefficients_vec(vec![Fr::one(), Fr::one()]);
        assert_eq!(quotient, expected_quotient);
    }
}
