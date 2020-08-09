use crate::constraint_system::StandardComposer;
use crate::constraint_system::Variable;
use algebra::{PairingEngine, Zero, One, Field, PrimeField, BigInteger};

impl<E: PairingEngine> StandardComposer<E> {
    /// Adds a boolean constraint (also known as binary constraint) where
    /// the gate eq. will enforce that the `Variable` received is either `0`
    /// or `1` by adding a constraint in the circuit.
    ///
    /// Note that using this constraint with whatever `Variable` that is not
    /// representing a value equalling 0 or 1, will always force the equation to fail.
    pub fn boolean_gate(&mut self, a: Variable) -> Variable {
        self.w_l.push(a);
        self.w_r.push(a);
        self.w_o.push(a);
        self.w_4.push(self.zero_var);

        self.q_m.push(E::Fr::one());
        self.q_l.push(E::Fr::zero());
        self.q_r.push(E::Fr::zero());
        self.q_o.push(-E::Fr::one());
        self.q_c.push(E::Fr::zero());
        self.q_4.push(E::Fr::zero());
        self.q_arith.push(E::Fr::one());

        self.q_range.push(E::Fr::zero());
        self.q_logic.push(E::Fr::zero());
        //self.q_ecc.push(E::Fr::zero());

        self.public_inputs.push(E::Fr::zero());

        self.perm
            .add_variables_to_map(a, a, a, self.zero_var, self.n);

        self.n += 1;

        a
    }
}
#[cfg(test)]
mod tests {
    use super::super::helper::*;
    use algebra::{
        bls12_381::{Fr, Bls12_381},
        Zero, One,
    };
    #[test]
    fn test_correct_bool_gate() {
        let res = gadget_tester::<Bls12_381>(
            |composer| {
                let zero = composer.add_input(Fr::zero());
                let one = composer.add_input(Fr::one());

                composer.boolean_gate(zero);
                composer.boolean_gate(one);
            },
            32,
        );
        assert!(res.is_ok())
    }

    #[test]
    fn test_incorrect_bool_gate() {
        let res = gadget_tester::<Bls12_381>(
            |composer| {
                let zero = composer.add_input(<Fr as From<u64>>::from(5));
                let one = composer.add_input(Fr::one());

                composer.boolean_gate(zero);
                composer.boolean_gate(one);
            },
            32,
        );
        assert!(res.is_err())
    }
}
