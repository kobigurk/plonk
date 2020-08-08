//! This is an extension over the [Merlin Transcript](https://docs.rs/merlin/1.0.3/merlin/struct.Transcript.html)
//! which adds a few extra functionalities.
use crate::commitment_scheme::kzg10::Commitment;
use algebra::{PairingEngine, PrimeField, Field, CanonicalSerialize};
use merlin::Transcript;

/// Transcript adds an abstraction over the Merlin transcript
/// For convenience
pub trait TranscriptProtocol<E: PairingEngine> {
    /// Append a `commitment` with the given `label`.
    fn append_commitment(&mut self, label: &'static [u8], comm: &Commitment<E>);

    /// Append a `Scalar` with the given `label`.
    fn append_scalar(&mut self, label: &'static [u8], s: &E::Fr);

    /// Compute a `label`ed challenge variable.
    fn challenge_scalar(&mut self, label: &'static [u8]) -> E::Fr;

    /// Append domain separator for the circuit size.
    fn circuit_domain_sep(&mut self, n: u64);
}

impl<E: PairingEngine> TranscriptProtocol<E> for Transcript {
    fn append_commitment(&mut self, label: &'static [u8], comm: &Commitment<E>) {
        let mut bytes = vec![];
        comm.0.serialize(&mut bytes).unwrap();
        self.append_message(label, &bytes);
    }

    fn append_scalar(&mut self, label: &'static [u8], s: &E::Fr) {
        let mut bytes = vec![];
        s.serialize(&mut bytes).unwrap();
        self.append_message(label, &bytes)
    }

    fn challenge_scalar(&mut self, label: &'static [u8]) -> E::Fr {
        let mut buf = [0u8; 64];
        self.challenge_bytes(label, &mut buf);

        E::Fr::from_random_bytes(&buf).unwrap()
    }

    fn circuit_domain_sep(&mut self, n: u64) {
        self.append_message(b"dom-sep", b"circuit_size");
        self.append_u64(b"n", n);
    }
}
