use ndarray::{Array1, Array2, arr2};
use tightbinding::float::is_near_float;
use atoms::{AtomPosition, LengthUnits, PositionList, Shift};
use layers::Layer;

/// vdW material multilayer.
///
/// `shift` must have one fewer entry than `layers`. The bottom layer is unshifted.
///
/// # TODO
///
/// Const generics would allow for the expected size of the `Vec` fields to be expressed at
/// compile-time using arrays.
pub struct Multilayer {
    /// Each layer in the multilayer.
    pub layers: Vec<Box<Layer>>,
    /// The displacement by which each layer is shifted relative to the bottom layer.
    /// (x, y) components give the lateral shift (allowing the various local configurations
    /// seen in the moire pattern) and the z component gives the total vertical displacement.
    pub shifts: Vec<Shift>,
    /// The unit cell of each layer is repeated [r1, r2] times in the [a1, a2] direction
    /// to form a supercell.
    ///
    /// r1 * a1 must be identical for each layer, and similarly for r2 * a2.
    pub repeat_count: Vec<[usize; 2]>,
}

impl Layer for Multilayer {
    /// # Panics
    ///
    /// Panics if self.shifts.len() != self.layers.len() - 1.
    ///
    /// Panics if all layer lattice vectors do not have the same units.
    ///
    /// Panics if supercell lattice vectors are not the same for each layer.
    fn atoms(&self) -> Vec<AtomPosition> {
        if self.shifts.len() != self.layers.len() - 1 {
            panic!("must have one fewer shift than number of layers (bottom layer is unshifted)");
        }
        self.verify_lat_vecs_constraint();

        let mut atoms = Vec::new();

        for layer_index in 0..self.layers.len() {
            let layer = &self.layers[layer_index];

            let layer_atoms_base = if layer_index == 0 {
                layer.atoms()
            } else {
                layer.atoms().shift(&self.shifts[layer_index - 1])
            };

            let (latvecs, units) = layer.lat_vecs();
            let layer_a1 = latvecs.slice(s![.., 0]);
            let layer_a2 = latvecs.slice(s![.., 1]);

            for r1 in 0..self.repeat_count[layer_index][0] {
                for r2 in 0..self.repeat_count[layer_index][1] {
                    let shift_xy = (r1 as f64) * &layer_a1 + (r2 as f64) * &layer_a2;
                    let shift = Shift {
                        cartesian: [shift_xy[0_usize], shift_xy[1_usize], 0.0],
                        units: units,
                    };

                    atoms.extend(layer_atoms_base.shift(&shift).iter().cloned());
                }
            }
        }

        atoms
    }

    /// # Panics
    ///
    /// Panics if all layer lattice vectors do not have the same units.
    ///
    /// Panics if supercell lattice vectors are not the same for each layer.
    fn lat_vecs(&self) -> (Array2<f64>, LengthUnits) {
        self.verify_lat_vecs_constraint();

        let (supercell_a1, supercell_a2, units) = self.supercell_lat_vecs(0);

        let d_supercell = arr2(&[
            [supercell_a1[0], supercell_a2[0]],
            [supercell_a1[1], supercell_a2[1]],
        ]);

        (d_supercell, units)
    }
}

impl Multilayer {
    fn supercell_lat_vecs(&self, layer: usize) -> (Array1<f64>, Array1<f64>, LengthUnits) {
        let (lat_vecs, units) = self.layers[layer].lat_vecs();
        let layer_a1 = lat_vecs.slice(s![.., 0]);
        let layer_a2 = lat_vecs.slice(s![.., 1]);

        let supercell_a1 = (self.repeat_count[layer][0] as f64) * &layer_a1;
        let supercell_a2 = (self.repeat_count[layer][1] as f64) * &layer_a2;

        (supercell_a1, supercell_a2, units)
    }

    fn verify_lat_vecs_constraint(&self) {
        let (supercell_a1_base, supercell_a2_base, units_base) = self.supercell_lat_vecs(0);

        let norm = |v: &Array1<f64>| (v[0].powi(2) + v[1].powi(2)).sqrt();

        let eps_rel = 2.0 * ::std::f64::EPSILON;
        let eps_abs_a1 = norm(&supercell_a1_base) * eps_rel;
        let eps_abs_a2 = norm(&supercell_a2_base) * eps_rel;

        for layer_index in 1..self.layers.len() {
            let (supercell_a1, supercell_a2, units) = self.supercell_lat_vecs(layer_index);

            if units != units_base {
                panic!("multilayer units not same for all layers");
            }

            for &(ref this_vec, ref base_vec, eps_abs) in [
                (&supercell_a1, &supercell_a1_base, eps_abs_a1),
                (&supercell_a2, &supercell_a2_base, eps_abs_a2),
            ].iter()
            {
                if !(is_near_float(this_vec[0], base_vec[0], eps_abs, eps_rel)
                    && is_near_float(this_vec[1], base_vec[1], eps_abs, eps_rel))
                {
                    panic!("Multilayer lattice vectors are inconsistent");
                }
            }
        }
    }
}
