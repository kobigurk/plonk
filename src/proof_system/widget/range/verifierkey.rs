use super::delta;
use crate::commitment_scheme::kzg10::Commitment;
use crate::proof_system::linearisation_poly::ProofEvaluations;
use algebra::{Field, PairingEngine};

#[derive(Debug)]
pub struct VerifierKey<E: PairingEngine> {
    pub q_range: Commitment<E>,
}

impl<E: PairingEngine> VerifierKey<E> {
    pub(crate) fn compute_linearisation_commitment(
        &self,
        range_separation_challenge: &E::Fr,
        scalars: &mut Vec<E::Fr>,
        points: &mut Vec<E::G1Affine>,
        evaluations: &ProofEvaluations<E::Fr>,
    ) {
        let four = <E::Fr as From<u64>>::from(4 as u64);

        let kappa = range_separation_challenge.square();
        let kappa_sq = kappa.square();
        let kappa_cu = kappa_sq * &kappa;

        let b_1 = delta(evaluations.c_eval - &(four * &evaluations.d_eval));
        let b_2 = delta(evaluations.b_eval - &(four * &evaluations.c_eval)) * &kappa;
        let b_3 = delta(evaluations.a_eval - &(four * &evaluations.b_eval)) * &kappa_sq;
        let b_4 = delta(evaluations.d_next_eval - &(four * &evaluations.a_eval)) * &kappa_cu;

        scalars.push((b_1 + &b_2 + &b_3 + &b_4) * range_separation_challenge);
        points.push(self.q_range.0);
    }
}
