//! A Proof stores the commitments to all of the elements that
//! are needed to univocally identify a prove of some statement.
//!
//! This module contains the implementation of the `StandardComposer`s
//! `Proof` structure and it's methods.

use super::linearisation_poly::ProofEvaluations;
use super::proof_system_errors::ProofErrors;
use crate::commitment_scheme::kzg10::{AggregateProof, Commitment, OpeningKey};
use crate::fft::EvaluationDomain;
use crate::proof_system::widget::VerifierKey;
use crate::transcript::TranscriptProtocol;
use anyhow::{Error, Result};
use merlin::Transcript;
use algebra::{PairingEngine, VariableBaseMSM, PrimeField, AffineCurve, ProjectiveCurve, Field};
use std::ops::{Add, Mul};

/// A Proof is a composition of `Commitments` to the witness, permutation,
/// quotient, shifted and opening polynomials as well as the
/// `ProofEvaluations`.
///
/// It's main goal is to have a `verify()` method attached which contains the
/// logic of the operations that the `Verifier` will need to do in order to
/// formally verify the `Proof`.
#[derive(Debug, Eq, PartialEq)]
pub struct Proof<E: PairingEngine> {
    /// Commitment to the witness polynomial for the left wires.
    pub a_comm: Commitment<E>,
    /// Commitment to the witness polynomial for the right wires.
    pub b_comm: Commitment<E>,
    /// Commitment to the witness polynomial for the output wires.
    pub c_comm: Commitment<E>,
    /// Commitment to the witness polynomial for the fourth wires.
    pub d_comm: Commitment<E>,

    /// Commitment to the permutation polynomial.
    pub z_comm: Commitment<E>,

    /// Commitment to the quotient polynomial.
    pub t_1_comm: Commitment<E>,
    /// Commitment to the quotient polynomial.
    pub t_2_comm: Commitment<E>,
    /// Commitment to the quotient polynomial.
    pub t_3_comm: Commitment<E>,
    /// Commitment to the quotient polynomial.
    pub t_4_comm: Commitment<E>,
    /// Commitment to the opening polynomial.
    pub w_z_comm: Commitment<E>,
    /// Commitment to the shifted opening polynomial.
    pub w_zw_comm: Commitment<E>,
    /// Subset of all of the evaluations added to the proof.
    pub evaluations: ProofEvaluations<E::Fr>,
}

impl<E: PairingEngine> Proof<E> {
    /// Performs the verification of a `Proof` returning a boolean result.
    pub(crate) fn verify(
        &self,
        verifier_key: &VerifierKey<E>,
        transcript: &mut Transcript,
        opening_key: &OpeningKey<E>,
        pub_inputs: &[E::Fr],
    ) -> Result<(), Error> {
        let domain = EvaluationDomain::new(verifier_key.n)?;

        // Subgroup checks are done when the proof is deserialised.

        // In order for the Verifier and Prover to have the same view in the non-interactive setting
        // Both parties must commit the same elements into the transcript
        // Below the verifier will simulate an interaction with the prover by adding the same elements
        // that the prover added into the transcript, hence generating the same challenges
        //
        // Add commitment to witness polynomials to transcript
        transcript.append_commitment(b"w_l", &self.a_comm);
        transcript.append_commitment(b"w_r", &self.b_comm);
        transcript.append_commitment(b"w_o", &self.c_comm);
        transcript.append_commitment(b"w_4", &self.d_comm);

        // Compute beta and gamma challenges
        let beta = TranscriptProtocol::<E>::challenge_scalar(transcript, b"beta");
        TranscriptProtocol::<E>::append_scalar(transcript, b"beta", &beta);
        let gamma = TranscriptProtocol::<E>::challenge_scalar(transcript, b"gamma");
        // Add commitment to permutation polynomial to transcript
        transcript.append_commitment(b"z", &self.z_comm);

        // Compute quotient challenge
        let alpha = TranscriptProtocol::<E>::challenge_scalar(transcript, b"alpha");
        let range_sep_challenge = TranscriptProtocol::<E>::challenge_scalar(transcript, b"range separation challenge");
        let logic_sep_challenge = TranscriptProtocol::<E>::challenge_scalar(transcript, b"logic separation challenge");
        let ecc_sep_challenge = TranscriptProtocol::<E>::challenge_scalar(transcript, b"ecc separation challenge");

        // Add commitment to quotient polynomial to transcript
        transcript.append_commitment(b"t_1", &self.t_1_comm);
        transcript.append_commitment(b"t_2", &self.t_2_comm);
        transcript.append_commitment(b"t_3", &self.t_3_comm);
        transcript.append_commitment(b"t_4", &self.t_4_comm);

        // Compute evaluation challenge
        let z_challenge = TranscriptProtocol::<E>::challenge_scalar(transcript, b"z");

        // Compute zero polynomial evaluated at `z_challenge`
        let z_h_eval = domain.evaluate_vanishing_polynomial(&z_challenge);

        // Compute first lagrange polynomial evaluated at `z_challenge`
        let l1_eval = compute_first_lagrange_evaluation(&domain, &z_h_eval, &z_challenge);

        // Compute quotient polynomial evaluated at `z_challenge`
        let t_eval = self.compute_quotient_evaluation(
            &domain,
            pub_inputs,
            &alpha,
            &beta,
            &gamma,
            &z_challenge,
            &z_h_eval,
            &l1_eval,
            &self.evaluations.perm_eval,
        );

        // Compute commitment to quotient polynomial
        // This method is necessary as we pass the `un-splitted` variation to our commitment scheme
        let t_comm = self.compute_quotient_commitment(&z_challenge, domain.size());

        // Add evaluations to transcript
        TranscriptProtocol::<E>::append_scalar(transcript, b"a_eval", &self.evaluations.a_eval);
        TranscriptProtocol::<E>::append_scalar(transcript, b"b_eval", &self.evaluations.b_eval);
        TranscriptProtocol::<E>::append_scalar(transcript, b"c_eval", &self.evaluations.c_eval);
        TranscriptProtocol::<E>::append_scalar(transcript, b"d_eval", &self.evaluations.d_eval);
        TranscriptProtocol::<E>::append_scalar(transcript, b"a_next_eval", &self.evaluations.a_next_eval);
        TranscriptProtocol::<E>::append_scalar(transcript, b"b_next_eval", &self.evaluations.b_next_eval);
        TranscriptProtocol::<E>::append_scalar(transcript, b"d_next_eval", &self.evaluations.d_next_eval);
        TranscriptProtocol::<E>::append_scalar(transcript, b"left_sig_eval", &self.evaluations.left_sigma_eval);
        TranscriptProtocol::<E>::append_scalar(transcript, b"right_sig_eval", &self.evaluations.right_sigma_eval);
        TranscriptProtocol::<E>::append_scalar(transcript, b"out_sig_eval", &self.evaluations.out_sigma_eval);
        TranscriptProtocol::<E>::append_scalar(transcript, b"q_arith_eval", &self.evaluations.q_arith_eval);
        TranscriptProtocol::<E>::append_scalar(transcript, b"q_c_eval", &self.evaluations.q_c_eval);
        //TranscriptProtocol::<E>::append_scalar(transcript, b"q_l_eval", &self.evaluations.q_l_eval);
        //TranscriptProtocol::<E>::append_scalar(transcript, b"q_r_eval", &self.evaluations.q_r_eval);
        TranscriptProtocol::<E>::append_scalar(transcript, b"perm_eval", &self.evaluations.perm_eval);
        TranscriptProtocol::<E>::append_scalar(transcript, b"t_eval", &t_eval);
        TranscriptProtocol::<E>::append_scalar(transcript, b"r_eval", &self.evaluations.lin_poly_eval);

        // Compute linearisation commitment
        let r_comm = self.compute_linearisation_commitment(
            &alpha,
            &beta,
            &gamma,
            (
                &range_sep_challenge,
                &logic_sep_challenge,
                &ecc_sep_challenge,
            ),
            &z_challenge,
            l1_eval,
            &verifier_key,
        );

        // Commitment Scheme
        // Now we delegate computation to the commitment scheme by batch checking two proofs
        // The `AggregateProof`, which is a proof that all the necessary polynomials evaluated at `z_challenge` are correct
        // and a `SingleProof` which is proof that the permutation polynomial evaluated at the shifted root of unity is correct

        // Compose the Aggregated Proof
        //
        let mut aggregate_proof = AggregateProof::with_witness(self.w_z_comm);
        aggregate_proof.add_part((t_eval, t_comm));
        aggregate_proof.add_part((self.evaluations.lin_poly_eval, r_comm));
        aggregate_proof.add_part((self.evaluations.a_eval, self.a_comm));
        aggregate_proof.add_part((self.evaluations.b_eval, self.b_comm));
        aggregate_proof.add_part((self.evaluations.c_eval, self.c_comm));
        aggregate_proof.add_part((self.evaluations.d_eval, self.d_comm));
        aggregate_proof.add_part((
            self.evaluations.left_sigma_eval,
            verifier_key.permutation.left_sigma,
        ));
        aggregate_proof.add_part((
            self.evaluations.right_sigma_eval,
            verifier_key.permutation.right_sigma,
        ));
        aggregate_proof.add_part((
            self.evaluations.out_sigma_eval,
            verifier_key.permutation.out_sigma,
        ));
        // Flatten proof with opening challenge
        let flattened_proof_a = aggregate_proof.flatten(transcript);

        // Compose the shifted aggregate proof
        let mut shifted_aggregate_proof = AggregateProof::with_witness(self.w_zw_comm);
        shifted_aggregate_proof.add_part((self.evaluations.perm_eval, self.z_comm));
        shifted_aggregate_proof.add_part((self.evaluations.a_next_eval, self.a_comm));
        shifted_aggregate_proof.add_part((self.evaluations.b_next_eval, self.b_comm));
        shifted_aggregate_proof.add_part((self.evaluations.d_next_eval, self.d_comm));
        let flattened_proof_b = shifted_aggregate_proof.flatten(transcript);

        // Add commitment to openings to transcript
        transcript.append_commitment(b"w_z", &self.w_z_comm);
        transcript.append_commitment(b"w_z_w", &self.w_zw_comm);

        // Batch check
        if opening_key
            .batch_check(
                &[z_challenge, (z_challenge * &domain.group_gen)],
                &[flattened_proof_a, flattened_proof_b],
                transcript,
            )
            .is_err()
        {
            return Err(ProofErrors::ProofVerificationError.into());
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn compute_quotient_evaluation(
        &self,
        domain: &EvaluationDomain<E::Fr>,
        pub_inputs: &[E::Fr],
        alpha: &E::Fr,
        beta: &E::Fr,
        gamma: &E::Fr,
        z_challenge: &E::Fr,
        z_h_eval: &E::Fr,
        l1_eval: &E::Fr,
        z_hat_eval: &E::Fr,
    ) -> E::Fr {
        // Compute the public input polynomial evaluated at `z_challenge`
        let pi_eval = compute_barycentric_eval(pub_inputs, z_challenge, domain);

        let alpha_sq = alpha.square();
        // r + PI(z)
        let a = self.evaluations.lin_poly_eval + &pi_eval;

        // a + beta * sigma_1 + gamma
        let beta_sig1 = *beta * &self.evaluations.left_sigma_eval;
        let b_0 = self.evaluations.a_eval + &beta_sig1 + gamma;

        // b+ beta * sigma_2 + gamma
        let beta_sig2 = *beta * &self.evaluations.right_sigma_eval;
        let b_1 = self.evaluations.b_eval + &beta_sig2 + gamma;

        // c+ beta * sigma_3 + gamma
        let beta_sig3 = *beta * &self.evaluations.out_sigma_eval;
        let b_2 = self.evaluations.c_eval + &beta_sig3 + gamma;

        // ((d + gamma) * z_hat) * alpha
        let b_3 = (self.evaluations.d_eval + gamma) * z_hat_eval * alpha;

        let b = b_0 * &b_1 * &b_2 * &b_3;

        // l_1(z) * alpha^2
        let c = *l1_eval * &alpha_sq;

        // Return t_eval
        (a - &b - &c) * &z_h_eval.inverse().unwrap()
    }

    fn compute_quotient_commitment(&self, z_challenge: &E::Fr, n: usize) -> Commitment<E> {
        let z_n = z_challenge.pow(&[n as u64, 0, 0, 0]);
        let z_two_n = z_challenge.pow(&[2 * n as u64, 0, 0, 0]);
        let z_three_n = z_challenge.pow(&[3 * n as u64, 0, 0, 0]);
        let t_comm = self.t_1_comm.0.into_projective()
            + &(self.t_2_comm.0.mul(z_n.into_repr()))
            + &(self.t_3_comm.0.mul(z_two_n.into_repr()))
            + &(self.t_4_comm.0.mul(z_three_n.into_repr()));
        Commitment::from_projective(t_comm)
    }

    // Commitment to [r]_1
    #[allow(clippy::too_many_arguments)]
    fn compute_linearisation_commitment(
        &self,
        alpha: &E::Fr,
        beta: &E::Fr,
        gamma: &E::Fr,
        (range_sep_challenge, logic_sep_challenge, ecc_sep_challenge): (&E::Fr, &E::Fr, &E::Fr),
        z_challenge: &E::Fr,
        l1_eval: E::Fr,
        verifier_key: &VerifierKey<E>,
    ) -> Commitment<E> {
        let mut scalars: Vec<_> = Vec::with_capacity(6);
        let mut points: Vec<E::G1Affine> = Vec::with_capacity(6);

        verifier_key.arithmetic.compute_linearisation_commitment(
            &mut scalars,
            &mut points,
            &self.evaluations,
        );

        verifier_key.range.compute_linearisation_commitment(
            &range_sep_challenge,
            &mut scalars,
            &mut points,
            &self.evaluations,
        );

        verifier_key.logic.compute_linearisation_commitment(
            &logic_sep_challenge,
            &mut scalars,
            &mut points,
            &self.evaluations,
        );

        /*
        verifier_key.ecc.compute_linearisation_commitment(
            &ecc_sep_challenge,
            &mut scalars,
            &mut points,
            &self.evaluations,
        );
         */

        verifier_key.permutation.compute_linearisation_commitment(
            &mut scalars,
            &mut points,
            &self.evaluations,
            z_challenge,
            (alpha, beta, gamma),
            &l1_eval,
            self.z_comm.0,
        );

        Commitment::from_projective(VariableBaseMSM::multi_scalar_mul(&points, &scalars.into_iter().map(|s| s.into_repr()).collect::<Vec<_>>()))
    }
}

fn compute_first_lagrange_evaluation<F: PrimeField>(
    domain: &EvaluationDomain<F>,
    z_h_eval: &F,
    z_challenge: &F,
) -> F {
    let n_fr = F::from(domain.size() as u128);
    let denom = n_fr * &(*z_challenge - &F::one());
    *z_h_eval * &denom.inverse().unwrap()
}

#[warn(clippy::needless_range_loop)]
fn compute_barycentric_eval<F: PrimeField>(
    evaluations: &[F],
    point: &F,
    domain: &EvaluationDomain<F>,
) -> F {
    use crate::util::batch_inversion;
    use rayon::iter::IntoParallelIterator;
    use rayon::prelude::*;

    let numerator = (point.pow(&[domain.size() as u64, 0, 0, 0]) - F::one()) * domain.size_inv;

    // Indices with non-zero evaluations
    let non_zero_evaluations: Vec<usize> = (0..evaluations.len())
        .into_par_iter()
        .filter(|&i| {
            let evaluation = &evaluations[i];
            evaluation != &F::zero()
        })
        .collect();

    // Only compute the denominators with non-zero evaluations
    let mut denominators: Vec<F> = (0..non_zero_evaluations.len())
        .into_par_iter()
        .map(|i| {
            // index of non-zero evaluation
            let index = non_zero_evaluations[i];

            (domain.group_gen_inv.pow(&[index as u64, 0, 0, 0]) * point) - F::one()
        })
        .collect();
    batch_inversion(&mut denominators);

    let result: F = (0..non_zero_evaluations.len())
        .into_par_iter()
        .map(|i| {
            let eval_index = non_zero_evaluations[i];
            let eval = evaluations[eval_index];

            denominators[i] * eval
        })
        .sum();

    result * numerator
}
