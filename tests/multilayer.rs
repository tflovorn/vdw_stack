extern crate vdw_stack;

use vdw_stack::atoms::{LengthUnits, Shift};
use vdw_stack::layers::Layer;
use vdw_stack::layers::tmd::{TMDStacking, TMD};
use vdw_stack::layers::flat::FlatHoneycomb;
use vdw_stack::layers::multilayer::Multilayer;

#[test]
pub fn construct_graphene_wse2_supercell() {
    let a_graphene = 2.46;
    let a_wse2 = (4.0 / 3.0) * a_graphene;
    let h_wse2 = 3.36;

    let c_bulk_wse2 = 12.976;
    let vertical_shift = c_bulk_wse2 / 2.0 - h_wse2 / 2.0;

    let graphene = FlatHoneycomb::new("C", "C", LengthUnits::Angstrom, a_graphene);
    let wse2 = TMD::new(
        "W",
        "Se",
        LengthUnits::Angstrom,
        a_wse2,
        h_wse2,
        TMDStacking::ABA,
    );

    let shifts = vec![
        Shift {
            cartesian: [0.0, 0.0, vertical_shift],
            units: LengthUnits::Angstrom,
        },
    ];
    let repeat_count = vec![[4, 4], [3, 3]];
    let expected_num_c = repeat_count[0][0] * repeat_count[0][1] * 2;
    let expected_num_w = repeat_count[1][0] * repeat_count[1][1];
    let expected_num_se = repeat_count[1][0] * repeat_count[1][1] * 2;

    let bilayer = Multilayer {
        layers: vec![Box::new(graphene), Box::new(wse2)],
        shifts,
        repeat_count,
    };

    let atoms = bilayer.atoms();

    assert_eq!(
        atoms.len(),
        expected_num_c + expected_num_w + expected_num_se
    );

    //for atom in bilayer.atoms() {
    //    println!("{:?}", atom);
    //}
}
