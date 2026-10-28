//! Feature extraction for the default P2Rank model:
//! features = [chem, volsite, protrusion, bfactor] + atom_table[apRawValids, apRawInvalids, atomicHydrophobicity]

use crate::atoms::{PointSet, P3};
use crate::elements::{ELEMENT_C, ELEMENT_N, ELEMENT_O};
use crate::protein::Protein;
use std::collections::HashMap;
use std::sync::OnceLock;

pub const NUM_FEATURES: usize = 35;
pub const IDX_CHEM_ATOMS: usize = 17;
pub const IDX_PROTRUSION: usize = 30;
pub const IDX_BFACTOR: usize = 31;
pub const IDX_ATOM_TABLE: usize = 32;

pub const FEATURE_HEADER: [&str; NUM_FEATURES] = [
    "chem.hydrophobic",
    "chem.hydrophilic",
    "chem.hydrophatyIndex",
    "chem.aliphatic",
    "chem.aromatic",
    "chem.sulfur",
    "chem.hydroxyl",
    "chem.basic",
    "chem.acidic",
    "chem.amide",
    "chem.posCharge",
    "chem.negCharge",
    "chem.hBondDonor",
    "chem.hBondAcceptor",
    "chem.hBondDonorAcceptor",
    "chem.polar",
    "chem.ionizable",
    "chem.atoms",
    "chem.atomDensity",
    "chem.atomC",
    "chem.atomO",
    "chem.atomN",
    "chem.hDonorAtoms",
    "chem.hAcceptorAtoms",
    "volsite.vsAromatic",
    "volsite.vsCation",
    "volsite.vsAnion",
    "volsite.vsHydrophobic",
    "volsite.vsAcceptor",
    "volsite.vsDonor",
    "protrusion.protrusion",
    "bfactor.bfactor",
    "atom_table.apRawValids",
    "atom_table.apRawInvalids",
    "atom_table.atomicHydrophobicity",
];

pub type FVec = [f64; NUM_FEATURES];

/// AminoAcidMapper.map in "minimal" mode
pub fn correct_residue_code(code: &str) -> String {
    if code.is_empty() {
        return code.to_string();
    }
    let upper = code.to_uppercase();
    match upper.as_str() {
        "MSE" => "MET".to_string(),
        "MEN" => "ASN".to_string(),
        _ => upper,
    }
}

// ---------------------------------------------------------------- chem

const AACODES: &[&str] =
    &["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "STP"];

fn mval(aa: &str, map: &[(&str, f64)]) -> f64 {
    map.iter().find(|e| e.0 == aa).map(|e| e.1).unwrap_or(0.0)
}

const HYDROPHOBIC: &[(&str, f64)] = &[
    ("PHE", 1.0),
    ("ILE", 1.0),
    ("TRP", 1.0),
    ("GLY", 1.0),
    ("LEU", 1.0),
    ("VAL", 1.0),
    ("MET", 1.0),
    ("ALA", 1.0),
    ("CYS", 1.0),
    ("TYR", 1.0),
    ("ARG", -1.0),
    ("ASN", -1.0),
    ("ASP", -1.0),
    ("GLN", -1.0),
    ("GLU", -1.0),
    ("LYS", -1.0),
    ("PRO", -1.0),
];
const HYDROPHATY_INDEX: &[(&str, f64)] = &[
    ("ALA", 1.8),
    ("ARG", -4.5),
    ("ASN", -3.5),
    ("ASP", -3.5),
    ("CYS", 2.5),
    ("GLU", -3.5),
    ("GLN", -3.5),
    ("GLY", -0.4),
    ("HIS", -3.2),
    ("ILE", 4.5),
    ("LEU", 3.8),
    ("LYS", -3.9),
    ("MET", 1.9),
    ("PHE", 2.8),
    ("PRO", -1.6),
    ("SER", -0.8),
    ("THR", -0.7),
    ("TRP", -0.9),
    ("TYR", -1.3),
    ("VAL", 4.2),
];
const ALIPHATIC: &[(&str, f64)] = &[("ALA", 1.0), ("LEU", 1.0), ("ILE", 1.0), ("VAL", 1.0), ("GLY", 1.0), ("PRO", 1.0)];
const AROMATIC: &[(&str, f64)] = &[("PHE", 1.0), ("TRP", 1.0), ("TYR", 1.0)];
const SULFUR: &[(&str, f64)] = &[("CYS", 1.0), ("MET", 1.0)];
const HYDROXYL: &[(&str, f64)] = &[("SER", 1.0), ("THR", 1.0)];
const BASIC: &[(&str, f64)] = &[("ARG", 3.0), ("LYS", 2.0), ("HIS", 1.0)];
const ACIDIC: &[(&str, f64)] = &[("ASP", 1.0), ("GLU", 1.0)];
const AMIDE: &[(&str, f64)] = &[("ASN", 1.0), ("GLN", 1.0)];
const CHARGE: &[(&str, f64)] = &[("ASP", -1.0), ("GLU", -1.0), ("ARG", 1.0), ("HIS", 1.0), ("LYS", 1.0)];
const POLAR: &[(&str, f64)] = &[
    ("ARG", 1.0),
    ("ASN", 1.0),
    ("ASP", 1.0),
    ("GLN", 1.0),
    ("GLU", 1.0),
    ("HIS", 1.0),
    ("LYS", 1.0),
    ("SER", 1.0),
    ("THR", 1.0),
    ("TYR", 1.0),
    ("CYS", 1.0),
];
const IONIZABLE: &[(&str, f64)] = &[("ASP", 1.0), ("GLU", 1.0), ("HIS", 1.0), ("LYS", 1.0), ("ARG", 1.0), ("CYS", 1.0), ("TYR", 1.0)];
const HB_DONOR: &[(&str, f64)] = &[("ARG", 1.0), ("LYS", 1.0), ("TRY", 1.0)];
const HB_ACCEPTOR: &[(&str, f64)] = &[("ASP", 1.0), ("GLU", 1.0)];
const HB_DONOR_ACCEPTOR: &[(&str, f64)] = &[("ASN", 1.0), ("GLN", 1.0), ("HIS", 1.0), ("SER", 1.0), ("THR", 1.0), ("TYR", 1.0)];

/// ChemDefaults.setAAProperties -> 24 chem values (index layout of ChemVector)
fn chem_defaults(aa: &str) -> [f64; 24] {
    let mut p = [0.0f64; 24];
    p[0] = mval(aa, HYDROPHOBIC);
    if p[0] < 0.0 {
        p[0] = 0.0;
    }
    p[1] = -mval(aa, HYDROPHOBIC);
    if p[1] < 0.0 {
        p[1] = 0.0;
    }
    p[2] = mval(aa, HYDROPHATY_INDEX);
    p[3] = mval(aa, ALIPHATIC);
    p[4] = mval(aa, AROMATIC);
    p[5] = mval(aa, SULFUR);
    p[6] = mval(aa, HYDROXYL);
    p[7] = mval(aa, BASIC);
    p[8] = mval(aa, ACIDIC);
    p[9] = mval(aa, AMIDE);
    let charge = mval(aa, CHARGE);
    if charge >= 0.0 {
        p[10] = charge;
    } else {
        p[11] = -charge;
    }
    p[12] = mval(aa, HB_DONOR);
    p[13] = mval(aa, HB_ACCEPTOR);
    p[14] = mval(aa, HB_DONOR_ACCEPTOR);
    p[15] = mval(aa, POLAR);
    p[16] = mval(aa, IONIZABLE);
    p
}

/// ChemVector.forAtom
fn chem_for_atom(element: u8, name: &str, code: &str) -> [f64; 24] {
    let mut p = if AACODES.contains(&code) { chem_defaults(code) } else { [0.0f64; 24] };
    // setFromResidueAtom
    p[18] = 1.0; // atomDensity
    let an = name;
    let (mut hd, mut ha) = (0.0f64, 0.0f64);
    match code {
        "ARG" => {
            if an == "NE" || an == "NH1" || an == "NH2" {
                hd += 1.0;
            }
        }
        "ASN" => {
            if an == "ND2" {
                hd += 1.0;
            } else if an == "OD1" {
                ha += 1.0;
            }
        }
        "ASP" => {
            if an == "OD1" || an == "OD2" {
                ha += 1.0;
            }
        }
        "GLN" => {
            if an == "NE2" {
                hd += 1.0;
            } else if an == "OE1" {
                ha += 1.0;
            }
        }
        "GLU" => {
            if an == "OE1" || an == "OE2" {
                ha += 1.0;
            }
        }
        "HIS" => {
            if an == "ND1" || an == "NE2" {
                hd += 1.0;
                ha += 1.0;
            }
        }
        "LYS" => {
            if an == "NZ" {
                hd += 1.0;
            }
        }
        "SER" => {
            if an == "OG" {
                hd += 1.0;
                ha += 1.0;
            }
        }
        "THR" => {
            if an == "OG1" {
                hd += 1.0;
                ha += 1.0;
            }
        }
        "TRP" => {
            if an == "NE1" {
                hd += 1.0;
            }
        }
        "TYR" => {
            if an == "OH" {
                hd += 1.0;
                ha += 1.0;
            }
        }
        _ => {}
    }
    p[22] += hd;
    p[23] += ha;
    if element == ELEMENT_C {
        p[19] = 1.0;
    } else if element == ELEMENT_O {
        p[20] = 1.0;
    } else if element == ELEMENT_N {
        p[21] = 1.0;
    } else if name.starts_with('C') {
        p[19] = 1.0;
    } else if name.starts_with('O') {
        p[20] = 1.0;
    } else if name.starts_with('N') {
        p[21] = 1.0;
    }
    p
}

// ---------------------------------------------------------------- volsite

/// VolSitePharmacophore.getAtomProperties -> [aromatic, cation, anion, hydrophobic, acceptor, donor]
fn volsite(atom_code: &str, residue_code: &str) -> [f64; 6] {
    let name = atom_code.to_uppercase();
    let resname = residue_code.to_uppercase();
    let name = name.as_str();
    let resname = resname.as_str();
    let app = |a: &str, list: &str| list.contains(a);
    const AROM: usize = 0;
    const CATION: usize = 1;
    const ANION: usize = 2;
    const HYD: usize = 3;
    const ACC: usize = 4;
    const DON: usize = 5;
    let mut r = [0.0f64; 6];
    let set = |r: &mut [f64; 6], idx: &[usize]| {
        for &i in idx {
            r[i] = 1.0;
        }
    };

    if name == "C" {
        set(&mut r, &[HYD]);
        return r;
    }
    if app(name, "CA, CB") && resname != "CA" {
        set(&mut r, &[HYD]);
        return r;
    }
    if name == "CD" && app(resname, "ARG, GLN, GLU, LYS, PRO") {
        set(&mut r, &[HYD]);
        return r;
    }
    if name == "CD1" {
        if app(resname, "ILE, LEU") {
            set(&mut r, &[HYD]);
            return r;
        } else if app(resname, "PHE, TRP, TYR") {
            set(&mut r, &[AROM]);
            return r;
        }
    }
    if name == "CD2" {
        if resname == "LEU" {
            set(&mut r, &[HYD]);
            return r;
        } else if app(resname, "PHE, TRP, TYR") {
            set(&mut r, &[AROM]);
            return r;
        } else if app(resname, "HIS,HID,HIE") {
            set(&mut r, &[AROM]);
            return r;
        }
    }
    if name == "CE" && app(resname, "LYS, MET") {
        set(&mut r, &[HYD]);
        return r;
    }
    if name == "CE1" {
        if app(resname, "HIS,HID,HIE") {
            set(&mut r, &[AROM]);
            return r;
        } else if app(resname, "PHE, TYR") {
            set(&mut r, &[AROM]);
            return r;
        }
    }
    if name == "CE2" && app(resname, "PHE, TRP, TYR") {
        set(&mut r, &[AROM]);
        return r;
    }
    if name == "CE3" && resname == "TRP" {
        set(&mut r, &[AROM]);
        return r;
    }
    if name == "CG" {
        if app(resname, "ARG, ASN, ASP, CYS, CYX, GLN, GLU, LEU, LYS, MET, PRO") {
            set(&mut r, &[HYD]);
            return r;
        } else if app(resname, "HIS, HID, HIE, PHE, TRP, TYR") {
            set(&mut r, &[AROM]);
            return r;
        }
    }
    if (name == "CG1" && app(resname, "ILE, VAL")) || (name == "CG2" && app(resname, "ILE, THR, VAL")) {
        set(&mut r, &[HYD]);
        return r;
    }
    if name == "CH2" && resname == "TRP" {
        set(&mut r, &[AROM]);
        return r;
    }
    if name == "CZ" {
        if resname == "ARG" {
            set(&mut r, &[HYD]);
            return r;
        } else if app(resname, "PHE, TYR") {
            set(&mut r, &[AROM]);
            return r;
        }
    }
    if app(name, "CZ2, CZ3") && resname == "TRP" {
        set(&mut r, &[AROM]);
        return r;
    }
    if name == "N" {
        set(&mut r, &[DON]);
        return r;
    }
    if name == "ND1" && app(resname, "HIS, HID, HIE") {
        set(&mut r, &[DON, ACC]);
        return r;
    }
    if name == "ND2" && resname == "ASN" {
        set(&mut r, &[DON]);
        return r;
    }
    if name == "NE" && app(resname, "ARG, LYS") {
        set(&mut r, &[CATION]);
        return r;
    }
    if name == "NE1" && resname == "TRP" {
        set(&mut r, &[DON]);
        return r;
    }
    if name == "NE2" {
        if resname == "GLN" {
            set(&mut r, &[DON]);
            return r;
        } else if app(resname, "HIS, HID, HIE") {
            set(&mut r, &[DON, ACC]);
            return r;
        }
    }
    if (app(name, "NH1, NH2") && resname == "ARG") || (name == "NZ" && resname == "LYS") {
        set(&mut r, &[CATION]);
        return r;
    }
    if name == "O" {
        set(&mut r, &[ACC]);
        return r;
    }
    if app(name, "OD1, OD2") && resname == "ASP" {
        set(&mut r, &[ANION]);
        return r;
    }
    if name == "OD1" && resname == "ASN" {
        set(&mut r, &[ACC]);
        return r;
    }
    if name == "OE1" && resname == "GLN" {
        set(&mut r, &[ACC]);
        return r;
    }
    if app(name, "OE1, OE2") && resname == "GLU" {
        set(&mut r, &[ANION]);
        return r;
    }
    if (name == "OG" && resname == "SER") || (name == "OG1" && resname == "THR") || (name == "OH" && resname == "TYR") {
        set(&mut r, &[DON, ACC]);
        return r;
    }
    if name == "OXT" {
        set(&mut r, &[ANION]);
        return r;
    }
    if app(name, "SD, SG") {
        if app(resname, "CYS, CYX") {
            set(&mut r, &[ACC]);
            return r;
        } else if resname == "MET" {
            set(&mut r, &[HYD]);
            return r;
        }
    }
    if app(name, "FE, MG, MG, MN, ZN, CO") {
        set(&mut r, &[CATION]);
        return r;
    }
    r
}

// ---------------------------------------------------------------- atom table

struct AtomTable {
    /// atomName -> [apRawValids, apRawInvalids, atomicHydrophobicity] (None = missing)
    values: HashMap<String, [Option<f64>; 3]>,
}

fn atom_table() -> &'static AtomTable {
    static T: OnceLock<AtomTable> = OnceLock::new();
    T.get_or_init(|| {
        let text = include_str!("data/atomic-properties.csv");
        let mut values = HashMap::new();
        let mut header: Option<Vec<String>> = None;
        for line in text.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let fields: Vec<String> = line.split(',').map(|s| s.trim().to_string()).collect();
            match &header {
                None => header = Some(fields[1..].to_vec()),
                Some(h) => {
                    let mut v = [None; 3];
                    for (k, prop) in ["apRawValids", "apRawInvalids", "atomicHydrophobicity"].iter().enumerate() {
                        if let Some(pi) = h.iter().position(|x| x == prop) {
                            if let Some(s) = fields.get(pi + 1) {
                                v[k] = crate::jutil::parse_java_double(s).filter(|x| !x.is_nan());
                            }
                        }
                    }
                    values.insert(fields[0].clone(), v);
                }
            }
        }
        AtomTable { values }
    })
}

/// Math.pow(x, 2.0) as computed by HotSpot (special-cased to x*x)
#[inline(always)]
pub fn pow2(x: f64) -> f64 {
    x * x
}

// ---------------------------------------------------------------- extractor

pub struct FeatureContext {
    /// surface layer (exposed) atoms with their per-atom feature vectors
    pub exposed: PointSet,
    /// per position in `exposed`: index into `atom_vectors` (resolving duplicate PDB serials)
    vec_index: Vec<u32>,
    atom_vectors: Vec<FVec>,
    pub deep_layer: PointSet,
}

fn atom_vector(p: &Protein, a: usize) -> FVec {
    let atom = &p.s.atoms[a];
    let code = correct_residue_code(p.s.atom_resname(a));
    let mut v = [0.0f64; NUM_FEATURES];
    let chem = chem_for_atom(atom.element, &atom.name, &code);
    v[..24].copy_from_slice(&chem);
    let vs = volsite(&atom.name, &code);
    v[24..30].copy_from_slice(&vs);
    v[IDX_BFACTOR] = atom.temp_factor as f64;
    let key = format!("{}.{}", code, atom.name);
    let row = atom_table().values.get(&key);
    for k in 0..3 {
        let val = row.and_then(|r| r[k]).unwrap_or(0.0);
        v[IDX_ATOM_TABLE + k] = pow2(val);
    }
    v
}

impl FeatureContext {
    pub fn new(p: &Protein, exposed: PointSet, deep_layer: PointSet) -> FeatureContext {
        let mut by_serial: HashMap<i32, u32> = HashMap::new();
        let mut atom_vectors = Vec::with_capacity(exposed.len());
        for (i, &a) in exposed.ids.iter().enumerate() {
            atom_vectors.push(atom_vector(p, a as usize));
            by_serial.insert(p.s.atoms[a as usize].serial, i as u32);
        }
        let vec_index = exposed.ids.iter().map(|&a| by_serial[&p.s.atoms[a as usize].serial]).collect();
        FeatureContext { exposed, vec_index, atom_vectors, deep_layer }
    }

    /// PrankFeatureExtractor.calcFeatureVector for a SAS point
    pub fn sas_vector(&self, point: &P3, scratch: &mut Vec<u32>) -> FVec {
        self.exposed.cutout_sphere_idx(point, 6.0, scratch);
        let mut res = [0.0f64; NUM_FEATURES];
        for &i in scratch.iter() {
            let i = i as usize;
            let props = &self.atom_vectors[self.vec_index[i] as usize];
            let a = &self.exposed.pts[i];
            let x = point[0] - a[0];
            let y = point[1] - a[1];
            let z = point[2] - a[2];
            let d = (x * x + y * y + z * z).sqrt();
            let weight = 1.0 - d / 6.0;
            for k in 0..NUM_FEATURES {
                res[k] += props[k] * weight;
            }
        }
        res[IDX_CHEM_ATOMS] = scratch.len() as f64;
        res[IDX_PROTRUSION] = self.deep_layer.count_within_sphere(point, 10.0) as f64;
        res
    }
}
