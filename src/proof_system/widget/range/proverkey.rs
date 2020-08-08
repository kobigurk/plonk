use super::delta;
use crate::fft::{Evaluations, Polynomial};
use algebra::{PrimeField, Field};

#[derive(Debug, Eq, PartialEq)]
pub struct ProverKey<F: PrimeField> {
    pub q_range: (Polynomial<F>, Evaluations<F>),
}

impl<F: PrimeField> ProverKey<F> {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn compute_quotient_i(
        &self,
        index: usize,
        range_separation_challenge: &F,
        w_l_i: &F,
        w_r_i: &F,
        w_o_i: &F,
        w_4_i: &F,
        w_4_i_next: &F,
    ) -> F {
        let four = <F as From<u64>>::from(4 as u64);
        let q_range_i = &self.q_range.1[index];

        let kappa = range_separation_challenge.square();
        let kappa_sq = kappa.square();
        let kappa_cu = kappa_sq * &kappa;

        // Delta([c(X) - 4 * d(X)]) + Delta([b(X) - 4 * c(X)]) + Delta([a(X) - 4 * b(X)]) + Delta([d(Xg) - 4 * a(X)]) * Q_Range(X)
        //
        let b_1 = delta(*w_o_i - four * w_4_i);
        let b_2 = delta(*w_r_i - four * w_o_i) * &kappa;
        let b_3 = delta(*w_l_i - four * w_r_i) * &kappa_sq;
        let b_4 = delta(*w_4_i_next - four * w_l_i) * &kappa_cu;
        (b_1 + &b_2 + &b_3 + &b_4) * q_range_i * range_separation_challenge
    }

    pub(crate) fn compute_linearisation(
        &self,
        range_separation_challenge: &F,
        a_eval: &F,
        b_eval: &F,
        c_eval: &F,
        d_eval: &F,
        d_next_eval: &F,
    ) -> Polynomial<F> {
        let four = <F as From<u64>>::from(4u64);
        let q_range_poly = &self.q_range.0;

        let kappa = range_separation_challenge.square();
        let kappa_sq = kappa.square();
        let kappa_cu = kappa_sq * kappa;

        // Delta([c_eval - 4 * d_eval]) + Delta([b_eval - 4 * c_eval]) + Delta([a_eval - 4 * b_eval]) + Delta([d_next_eval - 4 * a_eval]) * Q_Range(X)
        let b_1 = delta(*c_eval - &(four * d_eval));
        let b_2 = delta(*b_eval - &(four * c_eval)) * &kappa;
        let b_3 = delta(*a_eval - &(four * b_eval)) * &kappa_sq;
        let b_4 = delta(*d_next_eval - &(four * a_eval)) * &kappa_cu;

        let t = (b_1 + &b_2 + &b_3 + &b_4) * range_separation_challenge;

        q_range_poly * &t
    }
}
