use crate::fft::{EvaluationDomain, Polynomial};
use crate::proof_system::widget::ProverKey;
use algebra::{PairingEngine, PrimeField};

/// Evaluations at points `z` or and `z * root of unity`
pub struct Evaluations<F: PrimeField> {
    pub proof: ProofEvaluations<F>,
    // Evaluation of the linearisation sigma polynomial at `z`
    pub quot_eval: F,
}

/// Proof Evaluations is a subset of all of the evaluations. These evaluations will be added to the proof
#[derive(Debug, Eq, PartialEq)]
pub struct ProofEvaluations<F: PrimeField> {
    // Evaluation of the witness polynomial for the left wire at `z`
    pub a_eval: F,
    // Evaluation of the witness polynomial for the right wire at `z`
    pub b_eval: F,
    // Evaluation of the witness polynomial for the output wire at `z`
    pub c_eval: F,
    // Evaluation of the witness polynomial for the fourth wire at `z`
    pub d_eval: F,
    //
    pub a_next_eval: F,
    //
    pub b_next_eval: F,
    // Evaluation of the witness polynomial for the fourth wire at `z * root of unity`
    pub d_next_eval: F,
    // Evaluation of the arithmetic selector polynomial at `z`
    pub q_arith_eval: F,
    //
    pub q_c_eval: F,
    //
    //pub q_l_eval: F,
    //
    //pub q_r_eval: F,
    // Evaluation of the left sigma polynomial at `z`
    pub left_sigma_eval: F,
    // Evaluation of the right sigma polynomial at `z`
    pub right_sigma_eval: F,
    // Evaluation of the out sigma polynomial at `z`
    pub out_sigma_eval: F,

    // Evaluation of the linearisation sigma polynomial at `z`
    pub lin_poly_eval: F,

    // (Shifted) Evaluation of the permutation polynomial at `z * root of unity`
    pub perm_eval: F,
}

#[allow(clippy::too_many_arguments)]
/// Compute the linearisation polynomial
pub fn compute<E: PairingEngine>(
    domain: &EvaluationDomain<E::Fr>,
    prover_key: &ProverKey<E::Fr>,
    (
        alpha,
        beta,
        gamma,
        range_separation_challenge,
        logic_separation_challenge,
        ecc_separation_challenge,
        z_challenge,
    ): &(E::Fr, E::Fr, E::Fr, E::Fr, E::Fr, E::Fr, E::Fr),
    w_l_poly: &Polynomial<E::Fr>,
    w_r_poly: &Polynomial<E::Fr>,
    w_o_poly: &Polynomial<E::Fr>,
    w_4_poly: &Polynomial<E::Fr>,
    t_x_poly: &Polynomial<E::Fr>,
    z_poly: &Polynomial<E::Fr>,
) -> (Polynomial<E::Fr>, Evaluations<E::Fr>) {
    // Compute evaluations
    let quot_eval = t_x_poly.evaluate(z_challenge);
    let a_eval = w_l_poly.evaluate(z_challenge);
    let b_eval = w_r_poly.evaluate(z_challenge);
    let c_eval = w_o_poly.evaluate(z_challenge);
    let d_eval = w_4_poly.evaluate(z_challenge);
    let left_sigma_eval = prover_key.permutation.left_sigma.0.evaluate(z_challenge);
    let right_sigma_eval = prover_key.permutation.right_sigma.0.evaluate(z_challenge);
    let out_sigma_eval = prover_key.permutation.out_sigma.0.evaluate(z_challenge);
    let q_arith_eval = prover_key.arithmetic.q_arith.0.evaluate(z_challenge);
    let q_c_eval = prover_key.logic.q_c.0.evaluate(z_challenge);
    //let q_l_eval = prover_key.ecc.q_l.0.evaluate(z_challenge);
    //let q_r_eval = prover_key.ecc.q_r.0.evaluate(z_challenge);

    let a_next_eval = w_l_poly.evaluate(&(*z_challenge * &domain.group_gen));
    let b_next_eval = w_r_poly.evaluate(&(*z_challenge * &domain.group_gen));
    let d_next_eval = w_4_poly.evaluate(&(*z_challenge * &domain.group_gen));
    let perm_eval = z_poly.evaluate(&(*z_challenge * &domain.group_gen));

    let f_1 = compute_circuit_satisfiability::<E>(
        (
            range_separation_challenge,
            logic_separation_challenge,
            ecc_separation_challenge,
        ),
        &a_eval,
        &b_eval,
        &c_eval,
        &d_eval,
        &a_next_eval,
        &b_next_eval,
        &d_next_eval,
        &q_arith_eval,
        &q_c_eval,
        //&q_l_eval,
        //&q_r_eval,
        prover_key,
    );

    let f_2 = prover_key.permutation.compute_linearisation(
        z_challenge,
        (alpha, beta, gamma),
        (&a_eval, &b_eval, &c_eval, &d_eval),
        (&left_sigma_eval, &right_sigma_eval, &out_sigma_eval),
        &perm_eval,
        z_poly,
    );

    let lin_poly = &f_1 + &f_2;

    // Evaluate linearisation polynomial at z_challenge
    let lin_poly_eval = lin_poly.evaluate(z_challenge);

    (
        lin_poly,
        Evaluations {
            proof: ProofEvaluations {
                a_eval,
                b_eval,
                c_eval,
                d_eval,
                a_next_eval,
                b_next_eval,
                d_next_eval,
                q_arith_eval,
                q_c_eval,
                //q_l_eval,
                //q_r_eval,
                left_sigma_eval,
                right_sigma_eval,
                out_sigma_eval,
                lin_poly_eval,
                perm_eval,
            },
            quot_eval,
        },
    )
}

#[allow(clippy::too_many_arguments)]
fn compute_circuit_satisfiability<E: PairingEngine>(
    (range_separation_challenge, logic_separation_challenge, ecc_separation_challenge): (
        &E::Fr,
        &E::Fr,
        &E::Fr,
    ),
    a_eval: &E::Fr,
    b_eval: &E::Fr,
    c_eval: &E::Fr,
    d_eval: &E::Fr,
    a_next_eval: &E::Fr,
    b_next_eval: &E::Fr,
    d_next_eval: &E::Fr,
    q_arith_eval: &E::Fr,
    q_c_eval: &E::Fr,
    //q_l_eval: &E::Fr,
    //q_r_eval: &E::Fr,
    prover_key: &ProverKey<E::Fr>,
) -> Polynomial<E::Fr> {
    let a =
        prover_key
            .arithmetic
            .compute_linearisation(a_eval, b_eval, c_eval, d_eval, q_arith_eval);

    let b = prover_key.range.compute_linearisation(
        range_separation_challenge,
        a_eval,
        b_eval,
        c_eval,
        d_eval,
        &d_next_eval,
    );

    let c = prover_key.logic.compute_linearisation(
        logic_separation_challenge,
        a_eval,
        a_next_eval,
        b_eval,
        b_next_eval,
        c_eval,
        d_eval,
        d_next_eval,
        q_c_eval,
    );

    /*
    let d = prover_key.ecc.compute_linearisation(
        ecc_separation_challenge,
        a_eval,
        a_next_eval,
        b_eval,
        b_next_eval,
        c_eval,
        d_eval,
        d_next_eval,
        q_l_eval,
        q_r_eval,
        q_c_eval,
    );
     */

    let mut linearisation_poly = &a + &b;
    linearisation_poly += &c;
    //linearisation_poly += &d;

    linearisation_poly
}
