use ndarray::{Array2, arr2};
use atoms::{AtomPosition, LengthUnits};

pub mod tmd;
pub mod flat;
pub mod multilayer;

/// A vdW material layer, which knows how to emit a list of the atoms forming its unit cell.
/// Some conventions are assumed, to allow for interoperability between different materials:
///
/// * Materials have a triangular lattice with lattice vectors `a_1 = a (1/2, -\sqrt{3}/2)` and
/// `a_2 = a (1/2, \sqrt{3}/2)`.
/// * The triangular lattice A site is at in-plane position (x, y) = (0, 0).
/// * The atoms are centered around out-of-plane position z = 0.
pub trait Layer {
    /// A list of the atoms forming the unit cell.
    fn atoms(&self) -> Vec<AtomPosition>;

    /// In-plane lattice vectors of the layer unit cell.
    fn lat_vecs(&self) -> (Array2<f64>, LengthUnits);
}

/// The conventional matrix of lattice vectors for the triangular lattice,
/// ```
/// D = [[a1x, a2x,
///       a1y, a2y]]
/// ```
pub fn d_triangular(a: f64) -> Array2<f64> {
    arr2(&[
        [1.0 / 2.0, 1.0 / 2.0],
        [-3.0_f64.sqrt() / 2.0, 3.0_f64.sqrt() / 2.0],
    ]) * a
}
