//! BioJava element handling and CDK van der Waals radius lookup (as used by the surface).

use crate::elements_data::{CHAR_CDK_VDW, ELEMENT_CDK_VDW, ELEMENT_NAMES};

pub const ELEMENT_H: u8 = 0;
pub const ELEMENT_C: u8 = 1;
pub const ELEMENT_N: u8 = 2;
pub const ELEMENT_O: u8 = 3;
pub const ELEMENT_D: u8 = 4;
pub const ELEMENT_R: u8 = 105;

/// BioJava Element.valueOfIgnoreCase
pub fn element_value_of_ignore_case(sym: &str) -> Option<u8> {
    let l = sym.to_lowercase();
    ELEMENT_NAMES.iter().position(|n| n.to_lowercase() == l).map(|i| i as u8)
}

/// Van der Waals radius used by the surface calculation for an atom, following
/// CdkUtils.bioJavaToCDKAtom (element symbol -> fallback to first char of atom name -> "C")
/// and FasterNumericalSurface.getVdwRadius.
pub fn surface_vdw_radius(element: u8, atom_name: &str) -> Result<f64, String> {
    if let Some(bits) = ELEMENT_CDK_VDW[element as usize] {
        return Ok(f64::from_bits(bits));
    }
    // element unknown to CDK: use first char of the atom name
    let c = atom_name.chars().next().ok_or_else(|| "empty atom name".to_string())?;
    let code = c as u32;
    if (32..=126).contains(&code) {
        match CHAR_CDK_VDW[(code - 32) as usize] {
            0 => {}
            1 => return Err(format!("CDK cannot create atom with symbol '{}'", c)),
            bits => return Ok(f64::from_bits(bits)),
        }
    }
    // fallback "C"
    Ok(f64::from_bits(ELEMENT_CDK_VDW[ELEMENT_C as usize].unwrap()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_indices() {
        assert_eq!(ELEMENT_NAMES[ELEMENT_R as usize], "R");
        assert_eq!(ELEMENT_NAMES[ELEMENT_H as usize], "H");
        assert_eq!(element_value_of_ignore_case("SE"), Some(35));
        assert_eq!(ELEMENT_NAMES[35], "Se");
    }
}
