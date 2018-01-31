extern crate ndarray;
extern crate vdw_stack;

use ndarray::Array1;
use vdw_stack::enumerate::enumerate;

pub fn main() {
    let a_graphene = 2.46; // Angstrom [PRB 87, 205404]
    let a_wse2 = 3.32; // Angstrom [J Phys Chem C 119, 13169]

    let a1_hat = Array1::from_vec(vec![0.5, -0.5 * 3.0_f64.sqrt()]);
    let a2_hat = Array1::from_vec(vec![0.5, 0.5 * 3.0_f64.sqrt()]);

    let a1_graphene = a_graphene * &a1_hat;
    let a2_graphene = a_graphene * &a2_hat;
    let a1_wse2 = a_wse2 * &a1_hat;
    let a2_wse2 = a_wse2 * &a2_hat;

    let lat_vecs = [[a1_graphene, a2_graphene], [a1_wse2, a2_wse2]];

    let max_repeat = 10;
    let max_error = 0.0123;

    let result = enumerate(&lat_vecs, max_repeat, max_error);

    for v in result {
        println!("{:?}", v);
    }
}
