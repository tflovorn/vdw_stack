#[derive(Debug, Clone)]
pub struct AtomPosition {
    pub species: String,
    pub units: LengthUnits,
    pub cartesian: [f64; 3],
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum LengthUnits {
    Angstrom,
    Bohr,
}

impl AtomPosition {
    pub fn new(species: &str, units: LengthUnits, cartesian: [f64; 3]) -> AtomPosition {
        let species = String::from(species);
        AtomPosition {
            species,
            units,
            cartesian,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Shift {
    pub cartesian: [f64; 3],
    pub units: LengthUnits,
}

pub trait PositionList {
    /// Return a copy of `self` with each member shifted by `s`.
    fn shift(&self, s: &Shift) -> Self;
}

impl PositionList for Vec<AtomPosition> {
    /// # Panics
    ///
    /// Will panic if the `AtomPosition` units are not the same as the `Shift` units.
    ///
    /// # TODO
    ///
    /// Implement unit conversion to eliminate panic.
    fn shift(&self, s: &Shift) -> Vec<AtomPosition> {
        self.iter()
            .map(|ref pos| {
                if s.units != pos.units {
                    panic!(
                        "Incompatible units in shift -- shift {:?}, position {:?}",
                        s.units, pos.units
                    )
                }

                let shifted = [
                    pos.cartesian[0] + s.cartesian[0],
                    pos.cartesian[1] + s.cartesian[1],
                    pos.cartesian[2] + s.cartesian[2],
                ];

                AtomPosition {
                    species: pos.species.clone(),
                    units: pos.units,
                    cartesian: shifted,
                }
            })
            .collect()
    }
}
