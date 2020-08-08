use super::{delta, delta_xor_and};
use crate::commitment_scheme::kzg10::Commitment;
use crate::proof_system::linearisation_poly::ProofEvaluations;
use algebra::PairingEngine;

#[derive(Debug)]
pub struct VerifierKey<E: PairingEngine> {
    pub q_c: Commitment<E>,
    pub q_logic: Commitment<E>,
}

impl<E: PairingEngine> VerifierKey<E> {
    pub(crate) fn compute_linearisation_commitment(
        &self,
        logic_separation_challenge: &E::Fr,
        scalars: &mut Vec<E::Fr>,
        points: &mut Vec<E::G1Affine>,
        evaluations: &ProofEvaluations<E::Fr>,
    ) {
        let four = <E::Fr as From<u64>>::from(4 as u64);

        let kappa = *logic_separation_challenge * logic_separation_challenge;
        let kappa_sq = kappa * &kappa;
        let kappa_cu = kappa_sq * &kappa;
        let kappa_qu = kappa_cu * &kappa;

        let a = evaluations.a_next_eval - &(four * &evaluations.a_eval);
        let c_0 = delta(a);

        let b = evaluations.b_next_eval - &(four * &evaluations.b_eval);
        let c_1 = delta(b) * &kappa;

        let d = evaluations.d_next_eval - &(four * &evaluations.d_eval);
        let c_2 = delta(d) * &kappa_sq;

        let w = evaluations.c_eval;
        let c_3 = (w - &(a * &b)) * &kappa_cu;

        let c_4 = delta_xor_and(&a, &b, &w, &d, &evaluations.q_c_eval) * &kappa_qu;
        scalars.push((c_0 + &c_1 + &c_2 + &c_3 + &c_4) * logic_separation_challenge);
        points.push(self.q_logic.0);
    }
}
