// This file will contain the logic required
// to perform ECC operations in a PLONK
// circuit.

// For the scalar base operations we need to
// build a look up table, where we can find
// the values of particular indexes in 
#![allow(clippy::too_many_arguments)]
use super::PreProcessedPolynomial;
use crate::commitment_scheme::kzg10::Commitment;
use crate::fft::{Evaluations, Polynomial};
use crate::proof_system::linearisation_poly::ProofEvaluations;
use jubjub::{AffinePoint, GENERATOR, Fq, Fr};

#[cfg(feature = "serde")]
use serde::{de::Visitor, ser::SerializeStruct, Deserialize, Deserializer, Serialize, Serializer};

#[derive(Debug, Eq, PartialEq)]
pub struct ECCWidget {
    pub q_ecc: PreProcessedPolynomial,
    pub q_c: PreProcessedPolynomial,
}

#[cfg(feature = "serde")]
impl Serialize for ECCWidget {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut ecc_widget = serializer.serialize_struct("struct ECCWidget", 2)?;
        ecc_widget.serialize_field("q_c", &self.q_c)?;
        ecc_widget.serialize_field("q_eee", &self.q_ecc)?;
        ecc_widget.end()
    }
}

#[cfg(feature = "serde")]
impl<'de> Deserialize<'de> for LogicWidget {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        enum Field {
            Qc,
            Qecc,
        };

        impl<'de> Deserialize<'de> for Field {
            fn deserialize<D>(deserializer: D) -> Result<Field, D::Error>
            where
                D: Deserializer<'de>,
            {
                struct FieldVisitor;

                impl<'de> Visitor<'de> for FieldVisitor {
                    type Value = Field;

                    fn expecting(
                        &self,
                        formatter: &mut ::core::fmt::Formatter,
                    ) -> ::core::fmt::Result {
                        formatter.write_str("struct ECCWidget")
                    }

                    fn visit_str<E>(self, value: &str) -> Result<Field, E>
                    where
                        E: serde::de::Error,
                    {
                        match value {
                            "q_c" => Ok(Field::Qc),
                            "q_ecc" => Ok(Field::Qecc),
                            _ => Err(serde::de::Error::unknown_field(value, FIELDS)),
                        }
                    }
                }

                deserializer.deserialize_identifier(FieldVisitor)
            }
        }

        struct ECCWidgetVisitor;

        impl<'de> Visitor<'de> for ECCWidgetVisitor {
            type Value = ECCWidget;

            fn expecting(&self, formatter: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
                formatter.write_str("struct ECCWidget")
            }

            fn visit_seq<V>(self, mut seq: V) -> Result<LogicWidget, V::Error>
            where
                V: serde::de::SeqAccess<'de>,
            {
                let q_c = seq
                    .next_element()?
                    .ok_or_else(|| serde::de::Error::invalid_length(0, &self))?;
                let q_ecc = seq
                    .next_element()?
                    .ok_or_else(|| serde::de::Error::invalid_length(0, &self))?;
                Ok(LogicWidget { q_c, q_ecc })
            }
        }

        const FIELDS: &[&str] = &["q_c", "q_ecc"];
        deserializer.deserialize_struct("ECCWidget", FIELDS, ECCWidgetVisitor)
    }
}

impl ScalarMulWidget {
    pub(crate) fn new(
        q_c: (Polynomial, Commitment, Option<Evaluations>),
        q_ecc: (Polynomial, Commitment, Option<Evaluations>),
    ) -> ScalarMulWidget {
        ScalarMulWidget {
            q_c: PreProcessedPolynomial::new(q_ecc),
            q_ecc: PreProcessedPolynomial::new(q_ecc),
        }
    }
}

fn delta(f: Scalar) -> Scalar {
    let f_1 = f - Scalar::one();
    let f_2 = f - Scalar::from(2);
    let f_3 = f - Scalar::from(3);
    f * f_1 * f_2 * f_3
}

pub(crate) fn compute_quotient_i(
    &self,
    index: usize,
    ecc_consistency_challenge: &Scalar,
    w_l_i: &Scalar,
    w_l_i_next: &Scalar,
    w_r_i: &Scalar,
    w_r_i_next: &Scalar,
    w_o_i: &Scalar,
    w_o_i_next: &Scalar,
    w_4_i: &Scalar,
    w_4_i_next: &Scalar,
    q_o_i: &Scalar,
    q_ecc_i: &Scalar,
    q_1_i: &Scalar,
    q_2_i: &Scalar
) -> Scalar {

    let q_ecc_i = &self.q_ecc.evaluations.as_ref().unwrap()[index];
    let q_o_i = &self.q_o.evaluations.as_ref().unwrap()[index];
    
    let four = Scalar::from(4);
    let nine = Scalar::from(9);
    let one = Scalar::one();

    let kappa = ecc_consistency_challenge.square();
    let kappa_2 = kappa.square();
    let kappa_3 = kappa_2 * kappa;
    let kappa_4 = kappa_3 * kappa;
    let kappa_5 = kappa_4 * kappa;
    let kappa_6 = kappa_5 * kappa;
    let kappa_7 = kappa_6 * kappa;
    let kappa_8 = kappa_7 * kappa;
    let kappa_9 = kappa_8 * kappa;

    /// Compute the accumulator which tracks the current 
    /// rounds scalar multiplier, which is depedent on the 
    /// input bit
    let acc_input = four * w_4_i;
    let accum = w_4_i_next - acc_input;

    let accum_sqr = accum.square();
    
    /// To compute the y-alpha, which is the y-coordinate that corresponds to the x which is added 
    /// in each round then we use the formula below. This y-alpha is the y-coordianate that corresponds 
    /// to the y of one of the two points in the look up table, or the y in their inverses. 
    let a = w_o_i_next * q_o_i;
    let b = a + q_ecc_i;
    let y_alpha = b * accum;
    

    /// Check that the accumulator consistency at the identity element
    /// (accum - 1)(accum - 3)(accum + 1)(accum + 3) = 0 
    let a = accum_sqr - 9; 
    let b = accum_sqr - 1;
    let scalar_accum = a * b;
    let c_1 = scalar_accum * kappa;


    /// To compute x-alpha, which is the x-coordinate that we're adding in at each round. We need to
    /// explicit formualae with selector polynomials based on the values given in the lookup table.
    let a = accum_sqr * q_1_i;
    let b = a + q_2_i;
    let x_alpha_identity = b - w_o_i_next;
    let w_1 = x_alpha_identity * kappa_2;
    
    /// Consistency check of the x_accumulator
    let a = (w_l_i_next + w_l_i + w_o_i_next);
    let b = (w_o_i_next - w_l_i);
    let c = b.square();
    let d = y_alpha - w_r_i;
    let e = d.square();
    let x_accumulator = (a + c) - e;
    let a_1 = x_accumulator * kappa_3;

    /// Consistency check of the y_accumulator;
    let a = w_r_i_next - w_r_i;
    let b = w_o_i_next - w_l_i;
    let c = y_alpha - w_r_i;
    let d = w_l_i - w_l_i_next;
    let y_accumulator = (a + b) * (c + d);
    let b_1 = y_accumulator * kappa_4;

    /// Scalar accumulator consistency check;
    let a = w_4_i - 1 - w_o_i;
    let a = 

    

  
    
    // Check for consistency in the  x_value 

    let x_alpha_minus_x_one = w_o_i_next - w_l_i; 

    let 

    // Check for consistency in the y_value 

    let a_1 = 
}



 
   



    pub(crate) fn compute_linearisation(
        &self,
        a_eval: &Scalar,
        a_next_eval: &Scalar,
        b_eval: &Scalar,
        b_next_eval: &Scalar,
        c_eval: &Scalar,
        d_eval: &Scalar,
        d_next_eval: &Scalar,
        q_c_eval: &Scalar,
    ) -> Polynomial (){}


}

/// The polynomial identity to be evaluated, will check that the
/// initialiser has been done correctly and check that the accumulating 
/// values are correct. 
/// 
/// The identity checks that q_ecc * A = 0 
/// A = B + C
/// B = q_c * (c + let a + b)
/// C = (c1 + w1 + a1 + b1)
fn delta_ecc(c: &Fr, a: &Fr, b: &Fr, c_1: &Fr, w_1: &Fr, a_1: &Fr, b1_: &Fr, q_c: &Fr) -> {
    
    let B = q_c * (c + a + b);
    let C = (c_1 + w_1 + a_1 + b_1);
    B + C

}