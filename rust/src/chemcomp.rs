//! Chemical component dictionary, mirroring the parts of BioJava's ChemComp / ResidueType /
//! PolymerType logic that influence PDB parsing.
//!
//! By default the dictionary matches BioJava's `ReducedChemCompProvider` (the bundled set of
//! chem comps that BioJava falls back to when the RCSB download is not available). A full CCD
//! table produced by `tools/gen_chemcomp.py` can be loaded instead to mirror online behaviour.

use std::collections::HashMap;
use std::sync::OnceLock;

#[allow(dead_code)] // mirrors the Java data model
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PolymerType {
    Peptide,
    DPeptide,
    Dna,
    Rna,
    DnaRna,
    Polysaccharide,
    LPolysaccharide,
    OtherPolymer,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResidueType {
    Atomn,
    DPeptideLinking,
    LPeptideLinking,
    Glycine,
    PeptideLike,
    DPeptideAminoTerminus,
    LPeptideAminoTerminus,
    DPeptideCarboxyTerminus,
    LPeptideCarboxyTerminus,
    DnaLinking,
    RnaLinking,
    Dna3PrimeTerminus,
    Rna3PrimeTerminus,
    Dna5PrimeTerminus,
    Rna5PrimeTerminus,
    DSaccharide,
    DSaccharide14and14linking,
    DSaccharide14and16linking,
    LSaccharide,
    LSaccharide14and14linking,
    LSaccharide14and16linking,
    Saccharide,
    DBetaPeptideCGammaLinking,
    DGammaPeptideCDeltaLinking,
    LBetaPeptideCGammaLinking,
    LGammaPeptideCDeltaLinking,
    LDnaLinking,
    LRnaLinking,
    NonPolymer,
    OtherChemComp,
}

use PolymerType as P;
use ResidueType as R;

/// (residue type, polymer type, chem_comp_type string) in BioJava enum declaration order
const RESIDUE_TYPES: &[(ResidueType, Option<PolymerType>, &str)] = &[
    (R::Atomn, None, "null"),
    (R::DPeptideLinking, Some(P::DPeptide), "D-peptide linking"),
    (R::LPeptideLinking, Some(P::Peptide), "L-peptide linking"),
    (R::Glycine, Some(P::Peptide), "PEPTIDE LINKING"),
    (R::PeptideLike, Some(P::OtherPolymer), "peptide-like"),
    (R::DPeptideAminoTerminus, Some(P::DPeptide), "D-peptide NH3 amino terminus"),
    (R::LPeptideAminoTerminus, Some(P::Peptide), "L-peptide NH3 amino terminus"),
    (R::DPeptideCarboxyTerminus, Some(P::DPeptide), "D-peptide COOH carboxy terminus"),
    (R::LPeptideCarboxyTerminus, Some(P::Peptide), "L-peptide COOH carboxy terminus"),
    (R::DnaLinking, Some(P::Dna), "DNA linking"),
    (R::RnaLinking, Some(P::Rna), "RNA linking"),
    (R::Dna3PrimeTerminus, Some(P::Dna), "DNA OH 3 prime terminus"),
    (R::Rna3PrimeTerminus, Some(P::Rna), "RNA OH 3 prime terminus"),
    (R::Dna5PrimeTerminus, Some(P::Dna), "DNA OH 5 prime terminus"),
    (R::Rna5PrimeTerminus, Some(P::Rna), "RNA OH 5 prime terminus"),
    (R::DSaccharide, Some(P::Polysaccharide), "D-saccharide"),
    (R::DSaccharide14and14linking, Some(P::Polysaccharide), "D-saccharide 1,4 and 1,4 linking"),
    (R::DSaccharide14and16linking, Some(P::Polysaccharide), "D-saccharide 1,4 and 1,6 linking"),
    (R::LSaccharide, Some(P::LPolysaccharide), "L-saccharide"),
    (R::LSaccharide14and14linking, Some(P::LPolysaccharide), "L-saccharide 1,4 and 1,4 linking"),
    (R::LSaccharide14and16linking, Some(P::LPolysaccharide), "L-saccharide 1,4 and 1,6 linking"),
    (R::Saccharide, Some(P::Polysaccharide), "saccharide"),
    (R::DBetaPeptideCGammaLinking, Some(P::DPeptide), "D-beta-peptide, C-gamma linking"),
    (R::DGammaPeptideCDeltaLinking, Some(P::DPeptide), "D-gamma-peptide, C-delta linking"),
    (R::LBetaPeptideCGammaLinking, Some(P::Peptide), "L-beta-peptide, C-gamma linking"),
    (R::LGammaPeptideCDeltaLinking, Some(P::Peptide), "L-gamma-peptide, C-delta linking"),
    (R::LDnaLinking, Some(P::Dna), "L-DNA linking"),
    (R::LRnaLinking, Some(P::Dna), "L-RNA linking"),
    (R::NonPolymer, None, "non-polymer"),
    (R::OtherChemComp, None, "other"),
];

impl ResidueType {
    pub fn polymer_type(self) -> Option<PolymerType> {
        RESIDUE_TYPES.iter().find(|r| r.0 == self).unwrap().1
    }

    /// Port of BioJava ResidueType.getResidueTypeFromString
    pub fn from_chem_comp_type(t: &str) -> Option<ResidueType> {
        if t.eq_ignore_ascii_case("L-peptide linking") {
            return Some(R::LPeptideLinking);
        }
        // lookupTable contains both the original and the lower-cased chem_comp_type strings
        let lookup = |s: &str| -> Option<ResidueType> {
            for r in RESIDUE_TYPES {
                if r.2 == s || r.2.to_lowercase() == s {
                    return Some(r.0);
                }
            }
            None
        };
        if let Some(r) = lookup(t) {
            return Some(r);
        }
        let tl = t.to_lowercase();
        if let Some(r) = lookup(&tl) {
            return Some(r);
        }
        for r in RESIDUE_TYPES {
            let rl = r.2.to_lowercase();
            if r.2.eq_ignore_ascii_case(t) {
                return Some(r.0);
            }
            if rl.starts_with(&tl) {
                return Some(r.0);
            }
            if tl.starts_with(&rl) {
                return Some(r.0);
            }
        }
        None
    }
}

pub fn is_protein_only(p: Option<PolymerType>) -> bool {
    p == Some(P::Peptide)
}

pub fn is_polynucleotide_only(p: Option<PolymerType>) -> bool {
    matches!(p, Some(P::Dna) | Some(P::Rna) | Some(P::DnaRna))
}

#[derive(Debug, Clone)]
pub struct ChemComp {
    pub id: Option<String>,
    pub three_letter_code: Option<String>,
    pub residue_type: Option<ResidueType>,
    pub polymer_type: Option<PolymerType>,
    /// (atom_id, type_symbol)
    pub atoms: Vec<(String, String)>,
}

impl ChemComp {
    /// BioJava ChemComp.getEmptyChemComp() with id set
    pub fn empty(id: &str) -> ChemComp {
        ChemComp {
            id: Some(id.to_string()),
            three_letter_code: Some("???".to_string()),
            residue_type: Some(R::Atomn),
            polymer_type: Some(P::Unknown),
            atoms: Vec::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.id.is_none() || self.three_letter_code.is_none() || self.three_letter_code.as_deref() == Some("???")
    }
}

fn cif_value(v: &str) -> Option<String> {
    if v == "?" || v == "." {
        None
    } else {
        Some(v.to_string())
    }
}

pub struct ChemCompDictionary {
    comps: HashMap<String, ChemComp>,
}

impl ChemCompDictionary {
    pub fn parse_table(text: &str) -> ChemCompDictionary {
        let mut comps = HashMap::new();
        for line in text.lines() {
            if line.is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            let id = f[0].to_string();
            let ctype = cif_value(f[1]);
            let residue_type = ctype.as_deref().and_then(ResidueType::from_chem_comp_type);
            // ChemComp.setType: polymerType = residueType.polymerType when residueType != null
            let polymer_type = residue_type.and_then(|r| r.polymer_type());
            let three = cif_value(f.get(4).copied().unwrap_or("?"));
            let mut atoms = Vec::new();
            if let Some(a) = f.get(5) {
                for kv in a.split(';') {
                    if kv.is_empty() {
                        continue;
                    }
                    let (k, v) = kv.split_once('=').unwrap();
                    atoms.push((k.to_string(), v.to_string()));
                }
            }
            let cc = ChemComp { id: Some(id.clone()), three_letter_code: three, residue_type, polymer_type, atoms };
            comps.insert(id, cc);
        }
        ChemCompDictionary { comps }
    }

    pub fn reduced() -> &'static ChemCompDictionary {
        static D: OnceLock<ChemCompDictionary> = OnceLock::new();
        D.get_or_init(|| ChemCompDictionary::parse_table(include_str!("data/reduced_chemcomp.tsv")))
    }

    /// The full wwPDB Chemical Component Dictionary (types, codes, atoms), mirroring BioJava's
    /// `DownloadChemCompProvider` when RCSB is reachable.
    pub fn full() -> &'static ChemCompDictionary {
        static D: OnceLock<ChemCompDictionary> = OnceLock::new();
        D.get_or_init(|| {
            use std::io::Read;
            let gz: &[u8] = include_bytes!("data/ccd_full.tsv.gz");
            let mut text = String::new();
            flate2::read::GzDecoder::new(gz).read_to_string(&mut text).expect("embedded CCD table");
            ChemCompDictionary::parse_table(&text)
        })
    }

    /// ChemCompGroupFactory.getChemComp(recordName): never returns None for the reduced provider
    /// (unknown components get an empty chem comp).
    pub fn get(&self, record_name: &str) -> ChemComp {
        let name = crate::jutil::jtrim(&record_name.to_uppercase()).to_string();
        match self.comps.get(&name) {
            Some(cc) => cc.clone(),
            None => ChemComp::empty(&name),
        }
    }
}
