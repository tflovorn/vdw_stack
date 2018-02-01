use ndarray::{Array2, arr1};
use atoms::{AtomPosition, LengthUnits};
use layers::{d_triangular, Layer};

#[derive(Debug, Clone)]
pub struct FlatHoneycomb {
    /// Species of atom at the A site.
    pub a_species: String,
    /// Species of atom at the B site.
    pub b_species: String,
    /// Length units for `a`.
    pub units: LengthUnits,
    /// Lattice constant.
    pub a: f64,
}

impl FlatHoneycomb {
    pub fn new(a_species: &str, b_species: &str, units: LengthUnits, a: f64) -> FlatHoneycomb {
        let a_species = String::from(a_species);
        let b_species = String::from(b_species);

        FlatHoneycomb {
            a_species,
            b_species,
            units,
            a,
        }
    }
}

impl Layer for FlatHoneycomb {
    fn atoms(&self) -> Vec<AtomPosition> {
        let a_lat = arr1(&[0.0, 0.0]);
        let b_lat = arr1(&[1.0 / 3.0, 2.0 / 3.0]);
        let d = d_triangular(self.a);

        let a_cart = d.dot(&a_lat);
        let b_cart = d.dot(&b_lat);

        let a_atom = AtomPosition::new(&self.a_species, self.units, [a_cart[0], a_cart[1], 0.0]);
        let b_atom = AtomPosition::new(&self.b_species, self.units, [b_cart[0], b_cart[1], 0.0]);

        vec![a_atom, b_atom]
    }

    fn lat_vecs(&self) -> (Array2<f64>, LengthUnits) {
        (d_triangular(self.a), self.units)
    }
}

#[cfg(test)]
mod tests {
    use tightbinding::float::is_near_float;
    use atoms::LengthUnits;
    use layers::Layer;
    use super::FlatHoneycomb;

    #[test]
    fn flat_honeycomb_atom_conventions() {
        let a = 2.46;

        let atoms = FlatHoneycomb::new("C", "C", LengthUnits::Angstrom, a).atoms();

        let eps_rel = 2.0 * ::std::f64::EPSILON;
        let eps_abs = a * eps_rel;

        // All atoms are placed at z = 0.
        for pos in vec![&atoms[0], &atoms[1]] {
            assert!(pos.units == LengthUnits::Angstrom);
            assert!(is_near_float(pos.cartesian[2], 0.0, eps_abs, eps_rel));
        }

        // A and B atoms are in the expected right lateral positions.
        let a_pos = atoms[0].cartesian;
        assert!(is_near_float(a_pos[0], 0.0, eps_abs, eps_rel));
        assert!(is_near_float(a_pos[1], 0.0, eps_abs, eps_rel));

        let b_x = a / 2.0;
        let b_y = a / (2.0 * 3.0_f64.sqrt());
        let b_pos = atoms[1].cartesian;
        assert!(is_near_float(b_pos[0], b_x, eps_abs, eps_rel));
        assert!(is_near_float(b_pos[1], b_y, eps_abs, eps_rel));
    }
}
