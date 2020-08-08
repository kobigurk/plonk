#![allow(clippy::too_many_arguments)]

use crate::constraint_system::StandardComposer;
use crate::constraint_system::Variable;
use algebra::{PrimeField, PairingEngine, Zero, One};

impl<E: PairingEngine> StandardComposer<E> {
    /// Adds a width-3 add gate to the circuit, linking the addition of the
    /// provided inputs, scaled by the selector coefficients with the output
    /// provided.
    pub fn add_gate(
        &mut self,
        a: Variable,
        b: Variable,
        c: Variable,
        q_l: E::Fr,
        q_r: E::Fr,
        q_o: E::Fr,
        q_c: E::Fr,
        pi: E::Fr,
    ) -> Variable {
        self.big_add_gate(a, b, c, None, q_l, q_r, q_o, E::Fr::zero(), q_c, pi)
    }

    /// Adds a width-4 add gate to the circuit and it's corresponding
    /// constraint.
    ///
    /// This type of gate is usually used when we need to have
    /// the largest amount of performance and the minimum circuit-size
    /// possible. Since it allows the end-user to set every selector coefficient
    /// as scaling value on the gate eq.
    pub fn big_add_gate(
        &mut self,
        a: Variable,
        b: Variable,
        c: Variable,
        d: Option<Variable>,
        q_l: E::Fr,
        q_r: E::Fr,
        q_o: E::Fr,
        q_4: E::Fr,
        q_c: E::Fr,
        pi: E::Fr,
    ) -> Variable {
        // Check if advice wire has a value
        let d = match d {
            Some(var) => var,
            None => self.zero_var,
        };

        self.w_l.push(a);
        self.w_r.push(b);
        self.w_o.push(c);
        self.w_4.push(d);

        // For an add gate, q_m is zero
        self.q_m.push(E::Fr::zero());

        // Add selector vectors
        self.q_l.push(q_l);
        self.q_r.push(q_r);
        self.q_o.push(q_o);
        self.q_c.push(q_c);
        self.q_4.push(q_4);
        self.q_arith.push(E::Fr::one());
        self.q_range.push(E::Fr::zero());
        self.q_logic.push(E::Fr::zero());
        //self.q_ecc.push(E::Fr::zero());

        self.public_inputs.push(pi);

        self.perm.add_variables_to_map(a, b, c, d, self.n);

        self.n += 1;

        c
    }
    /// Adds a width-3 add gate to the circuit linking the product of the
    /// provided inputs scaled by the selector coefficient `q_m` with the output
    /// provided scaled by `q_o`.
    ///
    /// Note that this gate requires to provide the actual result of the gate
    /// (output wire) since it will just add a `mul constraint` to the circuit.
    pub fn mul_gate(
        &mut self,
        a: Variable,
        b: Variable,
        c: Variable,
        q_m: E::Fr,
        q_o: E::Fr,
        q_c: E::Fr,
        pi: E::Fr,
    ) -> Variable {
        self.big_mul_gate(a, b, c, None, q_m, q_o, q_c, E::Fr::zero(), pi)
    }

    /// Adds a width-4 `big_mul_gate` with the left, right and fourth inputs
    /// and it's scaling factors, computing & returning the output (result)
    /// `Variable` and adding the corresponding mul constraint.
    ///
    /// This type of gate is usually used when we need to have
    /// the largest amount of performance and the minimum circuit-size
    /// possible. Since it allows the end-user to setup all of the selector
    /// coefficients.
    ///
    /// Forces `q_l * (w_l + w_r) + w_4 * q_4 + + q_c + PI = q_o * w_o(computed by the gate)`.
    /// XXX: Maybe make these tuples instead of individual field?
    pub fn big_mul_gate(
        &mut self,
        a: Variable,
        b: Variable,
        c: Variable,
        d: Option<Variable>,
        q_m: E::Fr,
        q_o: E::Fr,
        q_c: E::Fr,
        q_4: E::Fr,
        pi: E::Fr,
    ) -> Variable {
        // Check if advice wire has a value
        let d = match d {
            Some(var) => var,
            None => self.zero_var,
        };

        self.w_l.push(a);
        self.w_r.push(b);
        self.w_o.push(c);
        self.w_4.push(d);

        // For a mul gate q_L and q_R is zero
        self.q_l.push(E::Fr::zero());
        self.q_r.push(E::Fr::zero());

        // Add selector vectors
        self.q_m.push(q_m);
        self.q_o.push(q_o);
        self.q_c.push(q_c);
        self.q_4.push(q_4);
        self.q_arith.push(E::Fr::one());

        self.q_range.push(E::Fr::zero());
        self.q_logic.push(E::Fr::zero());
        //self.q_ecc.push(E::Fr::zero());

        self.public_inputs.push(pi);

        self.perm.add_variables_to_map(a, b, c, d, self.n);

        self.n += 1;

        c
    }

    /// Adds a `big_addition_gate` with the left and right inputs
    /// and it's scaling factors, computing & returning the output (result)
    /// `Variable`, and adding the corresponding addition constraint.
    ///
    /// This type of gate is usually used when we don't need to have
    /// the largest amount of performance as well as the minimum circuit-size
    /// possible. Since it defaults some of the selector coeffs = 0 in order
    /// to reduce the verbosity and complexity.
    ///
    /// Forces `q_l * w_l + q_r * w_r + q_c + PI = w_o(computed by the gate)`.
    pub fn add(
        &mut self,
        q_l_a: (E::Fr, Variable),
        q_r_b: (E::Fr, Variable),
        q_c: E::Fr,
        pi: E::Fr,
    ) -> Variable {
        self.big_add(q_l_a, q_r_b, None, q_c, pi)
    }

    /// Adds a `big_addition_gate` with the left, right and fourth inputs
    /// and it's scaling factors, computing & returning the output (result)
    /// `Variable` and adding the corresponding addition constraint.
    ///
    /// This type of gate is usually used when we don't need to have
    /// the largest amount of performance and the minimum circuit-size
    /// possible. Since it defaults some of the selector coeffs = 0 in order
    /// to reduce the verbosity and complexity.
    ///
    /// Forces `q_l * w_l + q_r * w_r + q_4 * w_4 + q_c + PI = w_o(computed by the gate)`.
    pub fn big_add(
        &mut self,
        q_l_a: (E::Fr, Variable),
        q_r_b: (E::Fr, Variable),
        q_4_d: Option<(E::Fr, Variable)>,
        q_c: E::Fr,
        pi: E::Fr,
    ) -> Variable {
        // Check if advice wire is available
        let (q_4, d) = match q_4_d {
            Some((q_4, var)) => (q_4, var),
            None => (E::Fr::zero(), self.zero_var),
        };

        let (q_l, a) = q_l_a;
        let (q_r, b) = q_r_b;

        let q_o = -E::Fr::one();

        // Compute the output wire
        let a_eval = self.variables[&a];
        let b_eval = self.variables[&b];
        let d_eval = self.variables[&d];
        let c_eval = (q_l * &a_eval) + &(q_r * &b_eval) + &(q_4 * &d_eval) + &q_c + &pi;
        let c = self.add_input(c_eval);

        self.big_add_gate(a, b, c, Some(d), q_l, q_r, q_o, q_4, q_c, pi)
    }

    /// Adds a simple and basic addition to the circuit between to `Variable`s
    /// returning the resulting `Variable`.
    pub fn mul(
        &mut self,
        q_m: E::Fr,
        a: Variable,
        b: Variable,
        q_c: E::Fr,
        pi: E::Fr,
    ) -> Variable {
        self.big_mul(q_m, a, b, None, q_c, pi)
    }

    /// Adds a width-4 `big_mul_gate` with the left, right and fourth inputs
    /// and it's scaling factors, computing & returning the output (result)
    /// `Variable` and adding the corresponding mul constraint.
    ///
    /// This type of gate is usually used when we don't need to have
    /// the largest amount of performance and the minimum circuit-size
    /// possible. Since it defaults some of the selector coeffs = 0 in order
    /// to reduce the verbosity and complexity.
    ///
    /// Forces `q_l * (w_l + w_r) + w_4 * q_4 + q_c + PI = w_o(computed by the gate)`.
    /// XXX: This API is not consistent. It should use tuples and not individual fields
    pub fn big_mul(
        &mut self,
        q_m: E::Fr,
        a: Variable,
        b: Variable,
        q_4_d: Option<(E::Fr, Variable)>,
        q_c: E::Fr,
        pi: E::Fr,
    ) -> Variable {
        let q_o = -E::Fr::one();

        // Check if advice wire is available
        let (q_4, d) = match q_4_d {
            Some((q_4, var)) => (q_4, var),
            None => (E::Fr::zero(), self.zero_var),
        };

        // Compute output wire
        let a_eval = self.variables[&a];
        let b_eval = self.variables[&b];
        let d_eval = self.variables[&d];
        let c_eval = (q_m * &a_eval * &b_eval) + &(q_4 * &d_eval) + &q_c + &pi;
        let c = self.add_input(c_eval);

        self.big_mul_gate(a, b, c, Some(d), q_m, q_o, q_c, q_4, pi)
    }
}

#[cfg(test)]
mod tests {
    use super::super::helper::*;
    use dusk_bls12_381::F;

    #[test]
    fn test_public_inputs() {
        let res = gadget_tester(
            |composer| {
                let var_one = composer.add_input(E::Fr::one());

                let should_be_three = composer.big_add(
                    var_one.into(),
                    var_one.into(),
                    None,
                    E::Fr::zero(),
                    E::Fr::one(),
                );
                composer.constrain_to_constant(should_be_three, E::Fr::from(3), E::Fr::zero());
                let should_be_four = composer.big_add(
                    var_one.into(),
                    var_one.into(),
                    None,
                    E::Fr::zero(),
                    E::Fr::from(2),
                );
                composer.constrain_to_constant(should_be_four, E::Fr::from(4), E::Fr::zero());
            },
            200,
        );
        assert!(res.is_ok());
    }

    #[test]
    fn test_correct_add_mul_gate() {
        let res = gadget_tester(
            |composer| {
                // Verify that (4+5+5) * (6+7+7) = 280
                let four = composer.add_input(E::Fr::from(4));
                let five = composer.add_input(E::Fr::from(5));
                let six = composer.add_input(E::Fr::from(6));
                let seven = composer.add_input(E::Fr::from(7));

                let fourteen = composer.big_add(
                    four.into(),
                    five.into(),
                    Some(five.into()),
                    E::Fr::zero(),
                    E::Fr::zero(),
                );

                let twenty = composer.big_add(
                    six.into(),
                    seven.into(),
                    Some(seven.into()),
                    E::Fr::zero(),
                    E::Fr::zero(),
                );

                // There are quite a few ways to check the equation is correct, depending on your circumstance
                // If we already have the output wire, we can constrain the output of the mul_gate to be equal to it
                // If we do not, we can compute it using the `mul`
                // If the output is public, we can also constrain the output wire of the mul gate to it. This is what this test does
                let output = composer.mul(
                    E::Fr::one(),
                    fourteen,
                    twenty,
                    E::Fr::zero(),
                    E::Fr::zero(),
                );
                composer.constrain_to_constant(output, E::Fr::from(280), E::Fr::zero());
            },
            200,
        );
        assert!(res.is_ok());
    }

    #[test]
    fn test_correct_add_gate() {
        let res = gadget_tester(
            |composer| {
                let zero = composer.add_input(E::Fr::zero());
                let one = composer.add_input(E::Fr::one());

                let c = composer.add(
                    (E::Fr::one(), one),
                    (E::Fr::zero(), zero),
                    E::Fr::from(2u64),
                    E::Fr::zero(),
                );
                composer.constrain_to_constant(c, E::Fr::from(3), E::Fr::zero());
            },
            32,
        );
        assert!(res.is_ok())
    }

    #[test]
    fn test_correct_big_add_mul_gate() {
        let res = gadget_tester(
            |composer| {
                // Verify that (4+5+5) * (6+7+7) + (8*9) = 352
                let four = composer.add_input(E::Fr::from(4));
                let five = composer.add_input(E::Fr::from(5));
                let six = composer.add_input(E::Fr::from(6));
                let seven = composer.add_input(E::Fr::from(7));
                let nine = composer.add_input(E::Fr::from(9));

                let fourteen = composer.big_add(
                    four.into(),
                    five.into(),
                    Some(five.into()),
                    E::Fr::zero(),
                    E::Fr::zero(),
                );

                let twenty = composer.big_add(
                    six.into(),
                    seven.into(),
                    Some(seven.into()),
                    E::Fr::zero(),
                    E::Fr::zero(),
                );

                let output = composer.big_mul(
                    E::Fr::one(),
                    fourteen,
                    twenty,
                    Some((E::Fr::from(8), nine)),
                    E::Fr::zero(),
                    E::Fr::zero(),
                );
                composer.constrain_to_constant(output, E::Fr::from(352), E::Fr::zero());
            },
            200,
        );
        assert!(res.is_ok());
    }

    #[test]
    fn test_incorrect_add_mul_gate() {
        let res = gadget_tester(
            |composer| {
                // Verify that (5+5) * (6+7) != 117
                let five = composer.add_input(E::Fr::from(5));
                let six = composer.add_input(E::Fr::from(6));
                let seven = composer.add_input(E::Fr::from(7));

                let five_plus_five = composer.big_add(
                    five.into(),
                    five.into(),
                    None,
                    E::Fr::zero(),
                    E::Fr::zero(),
                );

                let six_plus_seven = composer.big_add(
                    six.into(),
                    seven.into(),
                    None,
                    E::Fr::zero(),
                    E::Fr::zero(),
                );

                let output = composer.mul(
                    E::Fr::one(),
                    five_plus_five,
                    six_plus_seven,
                    E::Fr::zero(),
                    E::Fr::zero(),
                );
                composer.constrain_to_constant(output, E::Fr::from(117), E::Fr::zero());
            },
            200,
        );
        assert!(res.is_err());
    }
}
