use algebra::PrimeField;

/// Constants used in the permutation argument to ensure that the wire subsets are disjoint.
pub(crate) fn K1<F: PrimeField>() -> F {
    F::from(7 as u64)
}
pub(crate) fn K2<F: PrimeField>() -> F {
    F::from(13 as u64)
}
pub(crate) fn K3<F: PrimeField>() -> F {
    F::from(17 as u64)
}