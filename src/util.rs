use rand_core::RngCore;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

/// Returns a vector of Frs of increasing powers of x from x^0 to x^d.
pub(crate) fn powers_of<F: PrimeField>(scalar: &F, max_degree: usize) -> Vec<F> {
    let mut powers = Vec::with_capacity(max_degree + 1);
    powers.push(F::one());
    for i in 1..=max_degree {
        powers.push(powers[i - 1] * scalar);
    }
    powers
}

/// Generates a random Fr using a RNG seed.
pub(crate) fn random_scalar<F: PrimeField, R: RngCore>(rng: &mut R) -> F {
    F::rand(rng)
}

/// Generates a random G1 Point using an RNG seed.
pub(crate) fn random_g1_point<E: PairingEngine, R: RngCore>(rng: &mut R) -> E::G1Projective {
    E::G1Projective::rand(rng)
}
/// Generates a random G2 point using an RNG seed.
pub(crate) fn random_g2_point<E: PairingEngine, R: RngCore>(rng: &mut R) -> E::G2Projective {
    E::G2Projective::rand(rng)
}

/// This function is only used to generate the SRS.
/// The intention is just to compute the resulting points
/// of the operation `a*P, b*P, c*P ... (n-1)*P` into a `Vec`.
pub(crate) fn slow_multiscalar_mul_single_base<E: PairingEngine>(
    scalars: &[E::Fr],
    base: E::G1Projective,
) -> Vec<E::G1Projective> {
    scalars.par_iter().map(|s| base.mul(s.into_repr())).collect()
}

// while we do not have batch inversion for scalars
use std::ops::MulAssign;
use algebra::{PrimeField, Field, Zero, One, PairingEngine, ProjectiveCurve, UniformRand};

pub fn batch_inversion<F: PrimeField>(v: &mut [F]) {
    // Montgomery’s Trick and Fast Implementation of Masked AES
    // Genelle, Prouff and Quisquater
    // Section 3.2

    // First pass: compute [a, ab, abc, ...]
    let mut prod = Vec::with_capacity(v.len());
    let mut tmp = F::one();
    for f in v.iter().filter(|f| f != &&F::zero()) {
        tmp.mul_assign(f);
        prod.push(tmp);
    }

    // Invert `tmp`.
    tmp = tmp.inverse().unwrap(); // Guaranteed to be nonzero.

    // Second pass: iterate backwards to compute inverses
    for (f, s) in v
        .iter_mut()
        // Backwards
        .rev()
        // Ignore normalized elements
        .filter(|f| f != &&F::zero())
        // Backwards, skip last element, fill in one for last term.
        .zip(prod.into_iter().rev().skip(1).chain(Some(F::one())))
    {
        // tmp := tmp * f; f := tmp * s = 1/f
        let new_tmp = tmp * *f;
        *f = tmp * s;
        tmp = new_tmp;
    }
}
#[cfg(test)]
mod test {
    use super::*;
    use algebra::bls12_381::Fr;
    #[test]
    fn test_batch_inversion() {
        let one = <Fr as From<u64>>::from(1);
        let two = <Fr as From<u64>>::from(2);
        let three = <Fr as From<u64>>::from(3);
        let four = <Fr as From<u64>>::from(4);
        let five = <Fr as From<u64>>::from(5);

        let original_scalars = vec![one, two, three, four, five];
        let mut inverted_scalars = vec![one, two, three, four, five];

        batch_inversion(&mut inverted_scalars);
        for (x, x_inv) in original_scalars.iter().zip(inverted_scalars.iter()) {
            assert_eq!(x.inverse().unwrap(), *x_inv);
        }
    }
}
