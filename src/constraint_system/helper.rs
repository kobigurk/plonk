use super::StandardComposer;
use crate::commitment_scheme::kzg10::PublicParameters;
use crate::proof_system::{Prover, Verifier};
use anyhow::{Error, Result};
use algebra::{PairingEngine, Zero, One, Field, PrimeField, AffineCurve, ProjectiveCurve};

/// Adds dummy constraints using arithmetic gates
pub(crate) fn dummy_gadget<E: PairingEngine>(n: usize, composer: &mut StandardComposer<E>) {
    let one = E::Fr::one();

    let var_one = composer.add_input(one);

    for _ in 0..n {
        composer.big_add(
            var_one.into(),
            var_one.into(),
            None,
            E::Fr::zero(),
            E::Fr::zero(),
        );
    }
}

/// Takes a generic gadget function with no auxillary input and
/// tests whether it passes an end-to-end test
pub(crate) fn gadget_tester<E: PairingEngine>(
    gadget: fn(composer: &mut StandardComposer<E>),
    n: usize,
) -> Result<(), Error> {
    // Common View
    let public_parameters = PublicParameters::setup(2 * n, &mut rand::thread_rng())?;
    // Provers View
    let (proof, public_inputs) = {
        // Create a prover struct
        let mut prover = Prover::new(b"demo");

        // Additionally key the transcript
        prover.key_transcript(b"key", b"additional seed information");

        // Add gadgets
        gadget(&mut prover.mut_cs());

        // Commit Key
        let (ck, _) = public_parameters.trim(2 * prover.cs.circuit_size().next_power_of_two())?;

        // Preprocess circuit
        prover.preprocess(&ck)?;

        // Once the prove method is called, the public inputs are cleared
        // So pre-fetch these before calling Prove
        let public_inputs = prover.cs.public_inputs.clone();

        // Compute Proof
        (prover.prove(&ck)?, public_inputs)
    };
    // Verifiers view
    //
    // Create a Verifier object
    let mut verifier = Verifier::new(b"demo");

    // Additionally key the transcript
    verifier.key_transcript(b"key", b"additional seed information");

    // Add gadgets
    gadget(&mut verifier.mut_cs());

    // Compute Commit and Verifier Key
    let (ck, vk) = public_parameters.trim(verifier.cs.circuit_size().next_power_of_two())?;

    // Preprocess circuit
    verifier.preprocess(&ck)?;

    // Verify proof
    verifier.verify(&proof, &vk, &public_inputs)
}
