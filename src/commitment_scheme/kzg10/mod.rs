//! Implementation of the KZG10 polynomial commitment scheme.
pub mod errors;
pub mod key;
pub mod srs;

pub use key::{CommitKey, OpeningKey};
pub use srs::PublicParameters;

use crate::transcript::TranscriptProtocol;
use crate::util::powers_of;
use merlin::Transcript;
use algebra::{PairingEngine, AffineCurve, Zero, One, ProjectiveCurve, Field, PrimeField};
use std::ops::Mul;

#[derive(Copy, Clone, Debug)]
/// Proof that a polynomial `p` was correctly evaluated at a point `z`
/// producing the evaluated point p(z).
pub struct Proof<E: PairingEngine> {
    /// This is a commitment to the witness polynomial.
    pub commitment_to_witness: Commitment<E>,
    /// This is the result of evaluating a polynomial at the point `z`.
    pub evaluated_point: E::Fr,
    /// This is the commitment to the polynomial that you want to prove a statement about.
    pub commitment_to_polynomial: Commitment<E>,
}

/// Proof that multiple polynomials were correctly evaluated at a point `z`,
/// each producing their respective evaluated points p_i(z).
#[derive(Debug)]
pub struct AggregateProof<E: PairingEngine> {
    /// This is a commitment to the aggregated witness polynomial.
    pub commitment_to_witness: Commitment<E>,
    /// These are the results of the evaluating each polynomial at the point `z`.
    pub evaluated_points: Vec<E::Fr>,
    /// These are the commitments to the polynomials which you want to prove a statement about.
    pub commitments_to_polynomials: Vec<Commitment<E>>,
}

impl<E: PairingEngine> AggregateProof<E> {
    /// Initialises an `AggregatedProof` with the commitment to the witness.
    pub fn with_witness(witness: Commitment<E>) -> AggregateProof<E> {
        AggregateProof {
            commitment_to_witness: witness,
            evaluated_points: Vec::new(),
            commitments_to_polynomials: Vec::new(),
        }
    }

    /// Adds an evaluated point with the commitment to the polynomial which produced it.
    pub fn add_part(&mut self, part: (E::Fr, Commitment<E>)) {
        self.evaluated_points.push(part.0);
        self.commitments_to_polynomials.push(part.1);
    }

    /// Flattens an `AggregateProof` into a `Proof`.
    /// The transcript must have the same view as the transcript that was used to aggregate the witness in the proving stage.
    pub fn flatten(&self, transcript: &mut Transcript) -> Proof<E> {
        let challenge = TranscriptProtocol::<E>::challenge_scalar(transcript, b"aggregate_witness");
        let powers = powers_of::<E::Fr>(&challenge, self.commitments_to_polynomials.len() - 1);

        // Flattened polynomial commitments using challenge
        let flattened_poly_commitments: E::G1Projective = self
            .commitments_to_polynomials
            .iter()
            .zip(powers.iter())
            .map(|(poly, challenge)| poly.0.mul(challenge.into_repr()))
            .sum();
        // Flattened evaluation points
        let flattened_poly_evaluations: E::Fr = self
            .evaluated_points
            .iter()
            .zip(powers.iter())
            .map(|(eval, challenge)| *eval * challenge)
            .fold(E::Fr::zero(), |acc, current_val| acc + &current_val);

        Proof {
            commitment_to_witness: self.commitment_to_witness,
            evaluated_point: flattened_poly_evaluations,
            commitment_to_polynomial: Commitment::from_projective(flattened_poly_commitments),
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
/// Holds a commitment to a polynomial in a form of a `G1Affine` Bls12_381 point.
pub struct Commitment<E: PairingEngine>(
    /// The commitment is a group element.
    pub E::G1Affine,
);

impl<E: PairingEngine> Commitment<E> {
    /// Builds a `Commitment` from a Bls12_381 `G1Projective` point.
    pub fn from_projective(g: E::G1Projective) -> Self {
        Self(g.into())
    }
    /// Builds a `Commitment` from a Bls12_381 `G1Affine` point.
    pub fn from_affine(g: E::G1Affine) -> Self {
        Self(g)
    }
    /// Builds an empty `Commitment` which is equivalent to the
    /// `G1Affine` identity point in Bls12_381.
    pub fn empty() -> Self {
        Commitment(E::G1Affine::prime_subgroup_generator())
    }
}

impl<E: PairingEngine> Default for Commitment<E> {
    fn default() -> Self {
        Commitment::empty()
    }
}
