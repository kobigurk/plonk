mod proverkey;
mod verifierkey;

pub use proverkey::ProverKey;
pub use verifierkey::VerifierKey;
use algebra::PrimeField;

/// Common functionality across both the ProverKey and VerifierKey are listed below
///
///
///
///

// Computes f(f-1)(f-2)(f-3)
fn delta<F: PrimeField>(f: F) -> F {
    let f_1 = f - F::one();
    let f_2 = f - F::from(2 as u64);
    let f_3 = f - F::from(3 as u64);
    f * f_1 * f_2 * f_3
}
