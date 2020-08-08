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

// The identity we want to check is q_logic * A = 0
// A = B + E
// B = q_c * [9c - 3(a+b)]
// E = 3(a+b+c) - 2F
// F = w[w(4w - 18(a+b) + 81) + 18(a^2 + b^2) - 81(a+b) + 83]
#[allow(non_snake_case)]
fn delta_xor_and<F: PrimeField>(a: &F, b: &F, w: &F, c: &F, q_c: &F) -> F {
    let nine = F::from(9 as u64);
    let two = F::from(2 as u64);
    let three = F::from(3 as u64);
    let four = F::from(4 as u64);
    let eighteen = F::from(18 as u64);
    let eighty_one = F::from(81 as u64);
    let eighty_three = F::from(83 as u64);

    let F = *w
        * &(*w * &(four * w - &(eighteen * &(*a + b)) + &eighty_one) + &(eighteen * &((*a * a) + &(*b * b)))
            - eighty_one * &(*a + b)
            + eighty_three);
    let E = three * &(*a + b + c) - &(two * &F);
    let B = *q_c * &((nine * c) - &(three * &(*a + b)));
    B + &E
}
