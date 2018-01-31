use atoms::AtomPosition;

pub mod tmd;

/// A vdW material layer, which knows how to emit a list of the atoms forming its unit cell.
/// Some conventions are assumed, to allow for interoperability between different materials:
///
/// * Materials have a triangular lattice with lattice vectors `a_1 = a (1/2, -\sqrt{3}/2)` and
/// `a_2 = a (1/2, \sqrt{3}/2)`.
/// * The triangular lattice A site is at in-plane position (x, y) = (0, 0).
/// * The atoms are centered around out-of-plane position z = 0.
pub trait Layer {
    fn atoms(&self) -> Vec<AtomPosition>;
}
