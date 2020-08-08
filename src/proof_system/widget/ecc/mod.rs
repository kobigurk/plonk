mod proverkey;
mod verifierkey;

pub use proverkey::ProverKey;
pub use verifierkey::VerifierKey;

// Note: The ECC gadget does not check that the initial point is on the curve for two reasons:
// - We constrain the accumulator to start from the identity point, which the verifier knows is on the curve
// - We are adding multiples of the generator to the accumulator which the verifier also knows is on the curve and is prime order
// - We do allow arbitrary Scalar multiplication, and possibly XXX: may add constraints to ensure the generator is correct (prime order)

// Bits are accumulated in base2. So we use d(Xw) - 2d(X) to extract the base2 bit

use algebra::PrimeField;

fn extract_bit<F: PrimeField>(curr_acc: F, next_acc: F) -> F {
    // Next - 2 * current
    next_acc - (curr_acc + curr_acc)
}

// Ensures that the bit is either +1, -1 or 0
fn check_bit_consistency<F: PrimeField>(bit: F) -> F {
    let one = F::one();
    bit * (bit - one) * (bit + one)
}
