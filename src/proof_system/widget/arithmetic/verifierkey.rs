use crate::commitment_scheme::kzg10::Commitment;
use crate::proof_system::linearisation_poly::ProofEvaluations;
use algebra::{PrimeField, PairingEngine};

#[derive(Debug)]
pub struct VerifierKey<E: PairingEngine> {
    pub q_m: Commitment<E>,
    pub q_l: Commitment<E>,
    pub q_r: Commitment<E>,
    pub q_o: Commitment<E>,
    pub q_c: Commitment<E>,
    pub q_4: Commitment<E>,
    pub q_arith: Commitment<E>,
}

impl<E: PairingEngine> VerifierKey<E> {
    pub(crate) fn compute_linearisation_commitment(
        &self,
        scalars: &mut Vec<E::Fr>,
        points: &mut Vec<E::G1Affine>,
        evaluations: &ProofEvaluations<E::Fr>,
    ) {
        let q_arith_eval = evaluations.q_arith_eval;

        scalars.push(evaluations.a_eval * &evaluations.b_eval * &q_arith_eval);
        points.push(self.q_m.0);

        scalars.push(evaluations.a_eval * &q_arith_eval);
        points.push(self.q_l.0);

        scalars.push(evaluations.b_eval * &q_arith_eval);
        points.push(self.q_r.0);

        scalars.push(evaluations.c_eval * &q_arith_eval);
        points.push(self.q_o.0);

        scalars.push(evaluations.d_eval * &q_arith_eval);
        points.push(self.q_4.0);

        scalars.push(q_arith_eval);
        points.push(self.q_c.0);
    }
}
