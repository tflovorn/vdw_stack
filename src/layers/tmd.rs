use ndarray::{Array2, arr1, arr2};
use atoms::{AtomPosition, LengthUnits, PositionList, Shift};
use layers::Layer;

/// Transition metal dichalcogenide in triangular lattice with minimal (3-atom) unit cell.
#[derive(Debug, Clone)]
pub struct TMD {
    /// Metal atom (Mo, W, ...).
    pub m_species: String,
    /// Chalcogen atom (S, Se, Te).
    pub x_species: String,
    /// Length units for `a` and `h`.
    pub units: LengthUnits,
    /// Lattice constant.
    pub a: f64,
    /// Vertical distance between chalcogen atoms.
    pub h: f64,
    /// Stacking configuration within the TMD layer. See `TMDStacking` for details.
    pub stacking_type: TMDStacking,
}

/// Stacking configuration within the TMD layer.
/// Value names give the triangular lattice unit cell positions of the [X, M, X] atoms,
/// in order from bottom to top.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum TMDStacking {
    ABA,
    BAB,
}

impl TMD {
    pub fn new(
        m_species: &str,
        x_species: &str,
        units: LengthUnits,
        a: f64,
        h: f64,
        stacking_type: TMDStacking,
    ) -> TMD {
        let m_species = String::from(m_species);
        let x_species = String::from(x_species);

        TMD {
            m_species,
            x_species,
            units,
            a,
            h,
            stacking_type,
        }
    }

    /// The matrix of lattice vectors,
    /// ```
    /// D = [[a1x, a2x,
    ///       a1y, a2y]]
    /// ```
    fn d(&self) -> Array2<f64> {
        arr2(&[
            [1.0 / 2.0, 1.0 / 2.0],
            [-3.0_f64.sqrt() / 2.0, 3.0_f64.sqrt() / 2.0],
        ]) * self.a
    }
}

impl Layer for TMD {
    fn atoms(&self) -> Vec<AtomPosition> {
        let a_lat = arr1(&[0.0, 0.0]);
        let b_lat = arr1(&[1.0 / 3.0, 2.0 / 3.0]);
        let d = self.d();

        let a_cart = d.dot(&a_lat);
        let b_cart = d.dot(&b_lat);

        let in_plane = match self.stacking_type {
            TMDStacking::ABA => [a_cart.clone(), b_cart.clone(), a_cart.clone()],
            TMDStacking::BAB => [b_cart.clone(), a_cart.clone(), b_cart.clone()],
        };

        let x_bot = AtomPosition::new(
            &self.x_species,
            self.units,
            [in_plane[0][0], in_plane[0][1], -self.h / 2.0],
        );
        let m = AtomPosition::new(
            &self.m_species,
            self.units,
            [in_plane[1][0], in_plane[1][1], 0.0],
        );
        let x_top = AtomPosition::new(
            &self.x_species,
            self.units,
            [in_plane[2][0], in_plane[2][1], self.h / 2.0],
        );

        vec![x_bot, m, x_top]
    }
}

#[cfg(test)]
mod tests {
    use tightbinding::float::is_near_float;
    use atoms::LengthUnits;
    use layers::Layer;
    use super::{TMDStacking, TMD};

    #[test]
    fn tmd_atom_conventions() {
        let a = 3.32;
        let h = 3.36;

        let atoms_aba = TMD::new("W", "Se", LengthUnits::Angstrom, a, h, TMDStacking::ABA).atoms();
        let atoms_bab = TMD::new("W", "Se", LengthUnits::Angstrom, a, h, TMDStacking::BAB).atoms();

        let eps_rel = 2.0 * ::std::f64::EPSILON;
        let eps_abs = a * eps_rel;

        // TMD uses the conventional lattice vectors and places A site at (0, 0).
        for pos in vec![&atoms_aba[0], &atoms_aba[2], &atoms_bab[1]] {
            assert!(pos.units == LengthUnits::Angstrom);
            assert!(is_near_float(pos.cartesian[0], 0.0, eps_abs, eps_rel));
            assert!(is_near_float(pos.cartesian[1], 0.0, eps_abs, eps_rel));
        }

        let b_x = a / 2.0;
        let b_y = a / (2.0 * 3.0_f64.sqrt());

        for pos in vec![&atoms_aba[1], &atoms_bab[0], &atoms_bab[2]] {
            assert!(pos.units == LengthUnits::Angstrom);
            assert!(is_near_float(pos.cartesian[0], b_x, eps_abs, eps_rel));
            assert!(is_near_float(pos.cartesian[1], b_y, eps_abs, eps_rel));
        }

        // TMD centers atoms around z = 0 in the out-of-plane direction.
        for atoms in vec![&atoms_aba, &atoms_bab] {
            let avg = atoms.iter().map(|ref pos| pos.cartesian[2]).sum::<f64>() / 3.0;
            assert!(is_near_float(avg, 0.0, eps_abs, eps_rel));
        }
    }
}

/// Transition metal dichalcogenide multilayer.
///
/// `shift` must have one fewer entry than `layers`. The bottom layer is unshifted.
///
/// # TODO
///
/// Const generics would allow for the expected size of the `Vec` fields to be expressed at
/// compile-time using arrays.
pub struct Multilayer {
    /// Each TMD layer in the multilayer.
    pub layers: Vec<TMD>,
    /// The displacement by which each layer is shifted relative to the bottom layer.
    /// (x, y) components give the lateral shift (allowing the various local configurations
    /// seen in the moire pattern) and the z component gives the total vertical displacement.
    pub shifts: Vec<Shift>,
}

impl Multilayer {
    /// Create a TMD bilayer where both layers are the same material.
    ///
    /// # Arguments
    ///
    /// * `bottom_layer` - The bottom layer of the bilayer.
    /// The top layer is identical, except that it is chosen to have the opposite stacking
    /// configuration (ABA <=> BAB) and is shifted upward.
    /// * `metal_distance` - vertical distance between the metal atoms.
    pub fn new_2h_bilayer(bottom_layer: &TMD, metal_distance: f64) -> Multilayer {
        let top_stacking = match bottom_layer.stacking_type {
            TMDStacking::ABA => TMDStacking::BAB,
            TMDStacking::BAB => TMDStacking::ABA,
        };
        let top_layer = TMD::new(
            &bottom_layer.m_species,
            &bottom_layer.x_species,
            bottom_layer.units,
            bottom_layer.a,
            bottom_layer.h,
            top_stacking,
        );

        let layers = vec![bottom_layer.clone(), top_layer];
        let shifts = vec![
            Shift {
                cartesian: [0.0, 0.0, metal_distance],
                units: bottom_layer.units,
            },
        ];

        Multilayer { layers, shifts }
    }
}

impl Layer for Multilayer {
    /// # Panics
    ///
    /// Panics if self.shifts.len() != self.layers.len() - 1.
    fn atoms(&self) -> Vec<AtomPosition> {
        if self.shifts.len() != self.layers.len() - 1 {
            panic!("must have one fewer shift than number of layers (bottom layer is unshifted)");
        }

        let mut atoms = Vec::new();

        for layer_index in 0..self.layers.len() {
            let layer = &self.layers[layer_index];

            let layer_atoms = if layer_index == 0 {
                layer.atoms()
            } else {
                layer.atoms().shift(&self.shifts[layer_index - 1])
            };

            atoms.extend(layer_atoms.iter().cloned());
        }

        atoms
    }
}

/// Given the vertical distance `metal_distance` between metal atoms of two adjacent TMD layers and
/// the vertical distance `h` between chalcogens within each layer, return the vertical distance
/// between the top chalcogen of the bottom layer and the bottom chalcogen of the top layer.
pub fn interlayer_distance(metal_distance: f64, h: [f64; 2]) -> f64 {
    metal_distance - h[0] / 2.0 - h[1] / 2.0
}

/// For two adjacent TMD layers, given the vertical distance `interlayer_distance` between the top
/// chalcogen of the bottom layer and the bottom chalcogen of the top layer and the vertical
/// distance `h` between chalcogens within each layer, return the vertical distance between metal atoms
/// of the two layers.
pub fn metal_distance(interlayer_distance: f64, h: [f64; 2]) -> f64 {
    h[0] / 2.0 + interlayer_distance + h[1] / 2.0
}
