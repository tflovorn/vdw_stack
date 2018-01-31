use ndarray::Array1;

/// Enumerate all rational approximant supercells for a pair of layers.
///
/// # Arguments
///
/// * lat_vecs - The lattice vectors of the two layers: `lat_vecs[i]` for the i'th layer.
/// * max_repeat - The largest multiple of lattice vectors to consider for constructing the
/// supercell.
/// * max_error - Include all supercells in the result where the relative difference between the
/// supercell lattice vectors for the two layers is greater than `max_error`.
///
/// # Panics
///
/// Will panic if max_error <= 0.
pub fn enumerate(
    lat_vecs: &[[Array1<f64>; 2]; 2],
    max_repeat: u16,
    max_error: f64,
) -> Vec<([[i64; 2]; 2], f64, f64, f64)> {
    if max_error <= 0.0 {
        panic!("must have `max_error` > 0");
    }

    let mut result = Vec::new();

    let min_bound = -(max_repeat as i64);
    let max_bound = (max_repeat as i64) + 1;

    for m0 in min_bound..max_bound {
        for n0 in min_bound..max_bound {
            for m1 in min_bound..max_bound {
                for n1 in min_bound..max_bound {
                    if (m0 == 0 && n0 == 0) || (m1 == 0 && n1 == 0) {
                        continue;
                    }

                    let v0 = (m0 as f64) * &lat_vecs[0][0] + (n0 as f64) * &lat_vecs[0][1];
                    let v1 = (m1 as f64) * &lat_vecs[1][0] + (n1 as f64) * &lat_vecs[1][1];
                    let v_diff = &v0 - &v1;

                    let rel_err = norm(&v_diff) / norm(&v0);

                    if rel_err < max_error {
                        result.push(([[m0, n0], [m1, n1]], norm(&v0), norm(&v1), rel_err));
                    }
                }
            }
        }
    }

    result
}

fn norm(v: &Array1<f64>) -> f64 {
    v.map(|x| *x * *x).scalar_sum().sqrt()
}
