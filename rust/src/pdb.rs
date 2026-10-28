//! PDB file parser reproducing the behaviour of BioJava's `PDBFileParser` (as configured by P2Rank)
//! for everything that influences P2Rank predictions: chains, groups (residues), alt locs, atoms,
//! group types, element assignment and entity (polymer / non-polymer / water) classification.

use crate::chemcomp::{is_polynucleotide_only, is_protein_only, ChemComp, ChemCompDictionary, ResidueType};
use crate::elements::{element_value_of_ignore_case, ELEMENT_R};
use crate::jutil::{is_java_space, jtrim, parse_java_double, parse_java_float, parse_java_int, JLine};
use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ResNum {
    pub chain: String,
    pub seq: i32,
    pub icode: Option<char>,
}

impl ResNum {
    /// ResidueNumber.toString()
    pub fn to_jstring(&self) -> String {
        match self.icode {
            Some(c) if c != ' ' => format!("{}{}", self.seq, c),
            _ => format!("{}", self.seq),
        }
    }

    /// ResidueNumber.printFull()
    pub fn print_full(&self) -> String {
        format!("{}_{}", self.chain, self.to_jstring())
    }
}

#[allow(dead_code)] // mirrors the Java data model
#[derive(Debug, Clone)]
pub struct Atom {
    pub serial: i32,
    pub name: String,
    pub alt_loc: char,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub occupancy: f32,
    pub temp_factor: f32,
    pub element: u8,
    /// index of the main group containing this atom
    pub group: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GroupKind {
    AminoAcid,
    Nucleotide,
    Hetatm,
}

#[derive(Debug, Clone)]
pub struct Group {
    pub pdb_name: String,
    pub resnum: ResNum,
    pub het_in_file: bool,
    pub kind: GroupKind,
    pub atoms: Vec<usize>,
    pub alt_locs: Vec<usize>,
    pub chain: Option<usize>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EntityType {
    Polymer,
    NonPolymer,
    Water,
}

#[derive(Debug, Clone)]
pub struct Chain {
    pub name: String,
    pub id: String,
    pub groups: Vec<usize>,
    pub entity: Option<EntityType>,
}

#[derive(Debug, Default)]
pub struct Structure {
    pub atoms: Vec<Atom>,
    pub groups: Vec<Group>,
    pub chains: Vec<Chain>,
    /// chains of each model, in BioJava order (polymers, non-polymers, waters)
    pub models: Vec<Vec<usize>>,
}

#[derive(Debug)]
pub struct ParseError(pub String);

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Outcome of a record handler: Java exceptions caught by the main loop (StringIndexOutOfBounds)
/// vs. fatal ones that abort parsing (e.g. NumberFormatException).
enum HErr {
    Caught,
    Fatal(String),
}

impl From<crate::jutil::IndexOutOfBounds> for HErr {
    fn from(_: crate::jutil::IndexOutOfBounds) -> Self {
        HErr::Caught
    }
}

const ATOM_CA_THRESHOLD: usize = 500000;

const COMPND_FIELD_VALUES: &[&str] = &[
    "MOL_ID:", "MOLECULE:", "CHAIN:", "SYNONYM:", "EC:", "FRAGMENT:", "ENGINEERED:", "MUTATION:", "BIOLOGICAL_UNIT:",
    "OTHER_DETAILS:",
];

const SOURCE_FIELD_VALUES: &[&str] = &[
    "ENGINEERED:",
    "MOL_ID:",
    "SYNTHETIC:",
    "FRAGMENT:",
    "ORGANISM_SCIENTIFIC:",
    "ORGANISM_COMMON:",
    "ORGANISM_TAXID:",
    "STRAIN:",
    "VARIANT:",
    "CELL_LINE:",
    "ATCC:",
    "ORGAN:",
    "TISSUE:",
    "CELL:",
    "ORGANELLE:",
    "SECRETION:",
    "GENE:",
    "CELLULAR_LOCATION:",
    "EXPRESSION_SYSTEM:",
    "EXPRESSION_SYSTEM_TAXID:",
    "EXPRESSION_SYSTEM_STRAIN:",
    "EXPRESSION_SYSTEM_VARIANT:",
    "EXPRESSION_SYSTEM_CELL_LINE:",
    "EXPRESSION_SYSTEM_ATCC_NUMBER:",
    "EXPRESSION_SYSTEM_ORGAN:",
    "EXPRESSION_SYSTEM_TISSUE:",
    "EXPRESSION_SYSTEM_CELL:",
    "EXPRESSION_SYSTEM_ORGANELLE:",
    "EXPRESSION_SYSTEM_CELLULAR_LOCATION:",
    "EXPRESSION_SYSTEM_VECTOR_TYPE:",
    "EXPRESSION_SYSTEM_VECTOR:",
    "EXPRESSION_SYSTEM_PLASMID:",
    "EXPRESSION_SYSTEM_GENE:",
    "OTHER_DETAILS:",
];

const WATER_NAMES: &[&str] = &["HOH", "DOD", "WAT"];

pub fn is_water_name(name: &str) -> bool {
    WATER_NAMES.contains(&name)
}

fn amino_acid_1letter(code3: &str) -> Option<char> {
    Some(match code3 {
        "GLY" => 'G',
        "ALA" => 'A',
        "VAL" => 'V',
        "LEU" => 'L',
        "ILE" => 'I',
        "PHE" => 'F',
        "TYR" => 'Y',
        "TRP" => 'W',
        "PRO" => 'P',
        "HIS" => 'H',
        "LYS" => 'K',
        "ARG" => 'R',
        "SER" => 'S',
        "THR" => 'T',
        "GLU" => 'E',
        "GLN" => 'Q',
        "ASP" => 'D',
        "ASN" => 'N',
        "CYS" => 'C',
        "MET" => 'M',
        "MSE" => 'M',
        "CSE" => 'U',
        "SEC" => 'U',
        "PYH" => 'O',
        "PYL" => 'O',
        _ => return None,
    })
}

const NUCLEOTIDES30: &[(&str, char)] = &[
    ("DA", 'A'),
    ("DC", 'C'),
    ("DG", 'G'),
    ("DT", 'T'),
    ("DI", 'I'),
    ("A", 'A'),
    ("G", 'G'),
    ("C", 'C'),
    ("U", 'U'),
    ("I", 'I'),
    ("TAF", 'X'),
    ("TC1", 'X'),
    ("TFE", 'X'),
    ("TFO", 'X'),
    ("TGP", 'X'),
    ("THX", 'X'),
    ("TLC", 'X'),
    ("TLN", 'X'),
    ("LCG", 'X'),
    ("TP1", 'X'),
    ("CP1", 'X'),
    ("TPN", 'X'),
    ("CPN", 'X'),
    ("GPN", 'X'),
    ("APN", 'X'),
    ("TPC", 'X'),
];
const NUCLEOTIDES23: &[&str] = &["C", "G", "A", "T", "U", "I", "+C", "+G", "+A", "+T", "+U", "+I"];

pub fn amino_acid_1letter_pub(code3: &str) -> Option<char> {
    amino_acid_1letter(code3)
}

pub fn is_nucleotide_pub(code3: &str) -> bool {
    is_nucleotide(code3)
}

/// StructureTools.isNucleotide
fn is_nucleotide(code3: &str) -> bool {
    let c = jtrim(code3);
    NUCLEOTIDES30.iter().any(|n| n.0 == c) || NUCLEOTIDES23.contains(&c)
}

/// StructureTools.get1LetterCode
fn get_1letter_code(code3: &str) -> char {
    if let Some(c) = amino_acid_1letter(code3) {
        return c;
    }
    let c3 = jtrim(code3);
    if is_nucleotide(c3) {
        if let Some(n) = NUCLEOTIDES30.iter().find(|n| n.0 == c3) {
            return n.1;
        }
        if NUCLEOTIDES23.contains(&c3) {
            return c3.chars().last().unwrap();
        }
        return 'X';
    }
    'X'
}

/// Java String.split("\\s+") semantics
fn java_split_ws(s: &str) -> Vec<String> {
    let mut parts: Vec<String> = Vec::new();
    let mut cur = String::new();
    let mut in_ws = false;
    let mut started = false;
    for c in s.chars() {
        if is_java_space(c) {
            if !in_ws {
                // a match begins: emit current token (including a leading empty one)
                parts.push(std::mem::take(&mut cur));
                in_ws = true;
            }
        } else {
            in_ws = false;
            cur.push(c);
        }
        started = true;
    }
    if !started {
        return vec![String::new()];
    }
    parts.push(cur);
    while let Some(last) = parts.last() {
        if last.is_empty() {
            parts.pop();
        } else {
            break;
        }
    }
    if parts.is_empty() && !s.chars().any(is_java_space) {
        return vec![s.to_string()];
    }
    parts
}

/// java.util.StringTokenizer with given delimiters
fn tokenize(s: &str, delims: &[char]) -> Vec<String> {
    s.split(|c| delims.contains(&c)).filter(|t| !t.is_empty()).map(|t| t.to_string()).collect()
}

const TOKENIZER_DEFAULT_DELIMS: &[char] = &[' ', '\t', '\n', '\r', '\u{0C}'];

struct Compound {
    mol_id: i32,
}

struct Parser<'d> {
    dict: &'d ChemCompDictionary,
    cc_cache: HashMap<String, ChemComp>,
    s: Structure,
    all_models: Vec<Vec<usize>>,
    current_model: Option<Vec<usize>>,
    current_chain: Option<usize>,
    current_group: Option<usize>,
    start_of_molecule: bool,
    start_of_model: bool,
    atom_count: usize,
    // header
    is_legacy_format: bool,
    pdb_id: Option<String>,
    // compound parsing
    compnd_lines: Vec<String>,
    source_lines: Vec<String>,
    entities: Vec<Compound>,
    compound_chain_ids: HashMap<i32, Vec<String>>,
    current_compound: Option<usize>,
    /// whether current_compound refers to a compound not (yet) added to `entities`
    pending_compound: Option<Compound>,
    prev_mol_id: i32,
    previous_continuation_field: String,
    continuation_field: String,
    continuation_string: String,
    is_last_compnd_line: bool,
    is_last_source_line: bool,
}

impl<'d> Parser<'d> {
    fn chem_comp(&mut self, name: &str) -> &ChemComp {
        let key = jtrim(&name.to_uppercase()).to_string();
        if !self.cc_cache.contains_key(&key) {
            let cc = self.dict.get(&key);
            self.cc_cache.insert(key.clone(), cc);
        }
        self.cc_cache.get(&key).unwrap()
    }

    /// PDBFileParser.getNewGroup
    fn new_group(&mut self, amino_code1: Option<char>, code3: &str, resnum: ResNum, het: bool) -> usize {
        let cc = self.chem_comp(code3);
        let kind = if !cc.is_empty() {
            if is_protein_only(cc.polymer_type) {
                GroupKind::AminoAcid
            } else if is_polynucleotide_only(cc.polymer_type) {
                GroupKind::Nucleotide
            } else {
                GroupKind::Hetatm
            }
        } else if amino_code1.is_none() || amino_code1 == Some('X') {
            GroupKind::Hetatm
        } else if is_nucleotide(code3) {
            GroupKind::Nucleotide
        } else {
            GroupKind::AminoAcid
        };
        self.s.groups.push(Group {
            pdb_name: code3.to_string(),
            resnum,
            het_in_file: het,
            kind,
            atoms: Vec::new(),
            alt_locs: Vec::new(),
            chain: None,
        });
        self.s.groups.len() - 1
    }

    fn add_group_to_chain(&mut self, chain: usize, group: usize) {
        self.s.groups[group].chain = Some(chain);
        let alts = self.s.groups[group].alt_locs.clone();
        for a in alts {
            self.s.groups[a].chain = Some(chain);
        }
        self.s.chains[chain].groups.push(group);
    }

    fn new_chain(&mut self, name: &str) -> usize {
        self.s.chains.push(Chain { name: name.to_string(), id: name.to_string(), groups: Vec::new(), entity: None });
        self.s.chains.len() - 1
    }

    fn group_has_atom(&self, g: usize, name: &str) -> bool {
        let name = jtrim(name);
        self.s.groups[g].atoms.iter().any(|&a| self.s.atoms[a].name == name)
    }

    fn add_atom(&mut self, g: usize, atom: usize) {
        self.s.atoms[atom].group = g;
        self.s.groups[g].atoms.push(atom);
    }

    /// PDBFileParser.getCorrectAltLocGroup
    fn correct_alt_loc_group(&mut self, alt_loc: char, amino_code1: Option<char>, code3: &str) -> usize {
        let cg = self.current_group.unwrap();
        if let Some(&a1) = self.s.groups[cg].atoms.first() {
            if self.s.atoms[a1].alt_loc == alt_loc {
                return cg;
            }
        }
        for &alt in &self.s.groups[cg].alt_locs {
            for &a in &self.s.groups[alt].atoms {
                if self.s.atoms[a].alt_loc == alt_loc {
                    return alt;
                }
            }
        }
        if code3 == self.s.groups[cg].pdb_name {
            if self.s.groups[cg].atoms.is_empty() {
                return cg;
            }
            let mut clone = self.s.groups[cg].clone();
            clone.atoms = Vec::new();
            clone.alt_locs = Vec::new();
            self.s.groups.push(clone);
            let id = self.s.groups.len() - 1;
            self.s.groups[cg].alt_locs.push(id);
            return id;
        }
        let resnum = self.s.groups[cg].resnum.clone();
        let het = false; // altLoc groups created here never get the het flag set in BioJava
        let id = self.new_group(amino_code1, code3, resnum, het);
        self.s.groups[cg].alt_locs.push(id);
        id
    }

    fn atom_handler(&mut self, line: &JLine) -> Result<(), HErr> {
        let chain_name = line.sub(21, 22)?;

        if let Some(cc) = self.current_chain {
            if self.s.chains[cc].name != chain_name {
                self.start_of_molecule = true;
            }
        }

        if self.start_of_molecule {
            if let Some(cc) = self.current_chain {
                self.current_model.as_mut().unwrap().push(cc);
                if let Some(cg) = self.current_group {
                    self.add_group_to_chain(cc, cg);
                }
            }
            let c = self.new_chain(&chain_name);
            self.current_chain = Some(c);
        }

        if self.start_of_model {
            if let Some(m) = self.current_model.take() {
                self.all_models.push(m);
            }
            self.current_model = Some(Vec::new());
        }

        let group_code3 = jtrim(&line.sub(17, 20)?).to_string();
        let res_num_str = jtrim(&line.sub(22, 26)?).to_string();
        let icode_c = line.sub(26, 27)?.chars().next().unwrap_or(' ');
        let icode = if icode_c == ' ' { None } else { Some(icode_c) };
        let seq = match parse_java_int(&res_num_str) {
            Some(v) => v,
            None => return Err(HErr::Fatal(format!("NumberFormatException: For input string: \"{}\"", res_num_str))),
        };
        let residue_number = ResNum { chain: chain_name.clone(), seq, icode };

        let mut amino_code1 = Some(get_1letter_code(&group_code3));
        let record_name = jtrim(&line.sub(0, 6)?).to_string();
        let mut is_het = false;
        if record_name == "HETATM" {
            if amino_code1 == Some('X') {
                amino_code1 = None;
            }
            is_het = true;
        }

        if self.start_of_molecule {
            let g = self.new_group(amino_code1, &group_code3, residue_number.clone(), is_het);
            self.current_group = Some(g);
        }

        self.start_of_model = false;
        self.start_of_molecule = false;

        let alt_loc = line.sub(16, 17)?.chars().next().unwrap_or(' ');
        let mut alt_group: Option<usize> = None;

        let cg = self.current_group.unwrap();
        if residue_number != self.s.groups[cg].resnum {
            let chain = self.current_chain.unwrap();
            self.add_group_to_chain(chain, cg);
            let g = self.new_group(amino_code1, &group_code3, residue_number.clone(), is_het);
            self.current_group = Some(g);
        } else if alt_loc != ' ' {
            let ag = self.correct_alt_loc_group(alt_loc, amino_code1, &group_code3);
            if self.s.groups[ag].chain.is_none() {
                self.s.groups[ag].chain = self.current_chain;
            }
            alt_group = Some(ag);
        }

        self.atom_count += 1;
        if self.atom_count == ATOM_CA_THRESHOLD {
            return Err(HErr::Fatal(format!("structures with more than {} atoms are not supported", ATOM_CA_THRESHOLD)));
        }

        let fullname = line.sub(12, 16)?;

        let serial_str = jtrim(&line.sub(6, 11)?).to_string();
        let serial = match parse_java_int(&serial_str) {
            Some(v) => v,
            None => return Err(HErr::Fatal(format!("NumberFormatException: For input string: \"{}\"", serial_str))),
        };

        let parse_coord = |s: String| -> Result<f64, HErr> {
            parse_java_double(&s).ok_or_else(|| {
                if s.is_empty() {
                    HErr::Fatal("NumberFormatException: empty String".to_string())
                } else {
                    HErr::Fatal(format!("NumberFormatException: For input string: \"{}\"", s))
                }
            })
        };
        let x = parse_coord(jtrim(&line.sub(30, 38)?).to_string())?;
        let y = parse_coord(jtrim(&line.sub(38, 46)?).to_string())?;
        let z = parse_coord(jtrim(&line.sub(46, 54)?).to_string())?;

        let mut occu = 1.0f32;
        if line.len() > 59 {
            if let Some(v) = parse_java_float(jtrim(&line.sub(54, 60)?)) {
                occu = v;
            }
        }
        let mut tempf = 0.0f32;
        if line.len() > 65 {
            if let Some(v) = parse_java_float(jtrim(&line.sub(60, 66)?)) {
                tempf = v;
            }
        }

        let name = jtrim(&fullname).to_string();
        let mut element = ELEMENT_R;
        let mut guess = true;
        if line.len() > 77 {
            let sym = jtrim(&line.sub(76, 78)?).to_string();
            if !sym.is_empty() {
                if let Some(e) = element_value_of_ignore_case(&sym) {
                    element = e;
                    guess = false;
                }
            }
        }
        if guess {
            let cg = self.current_group.unwrap();
            let gname = self.s.groups[cg].pdb_name.clone();
            let cc = self.chem_comp(&gname);
            if let Some((_, sym)) = cc.atoms.iter().find(|(id, _)| *id == name) {
                if let Some(e) = element_value_of_ignore_case(sym) {
                    element = e;
                }
            }
        }

        let cg = self.current_group.unwrap();
        self.s.atoms.push(Atom {
            serial,
            name: name.clone(),
            alt_loc,
            x,
            y,
            z,
            occupancy: occu,
            temp_factor: tempf,
            element,
            group: cg,
        });
        let atom = self.s.atoms.len() - 1;

        if let Some(ag) = alt_group {
            self.add_atom(ag, atom);
        } else {
            self.add_atom(cg, atom);
        }
        if !self.group_has_atom(cg, &name) {
            self.add_atom(cg, atom);
        }
        Ok(())
    }

    fn header_handler(&mut self, line: &JLine) -> Result<(), HErr> {
        let trimmed_len = jtrim(line.as_str()).encode_utf16().count();
        if trimmed_len > 62 {
            let code = jtrim(&line.sub(62, trimmed_len.min(66))?).to_string();
            self.pdb_id = Some(code);
        }
        if trimmed_len > 66 {
            let id = self.pdb_id.clone().unwrap_or_default();
            if id == line.sub(72, 76)? {
                self.is_legacy_format = true;
            }
        }
        Ok(())
    }

    // ---------------- COMPND / SOURCE ----------------

    fn compnd_handler(&mut self, line0: &str) -> Result<(), String> {
        let mut line = JLine::new(line0.to_string());
        if self.is_legacy_format {
            line = JLine::new(line.sub(0, 72).map_err(|_| "COMPND line too short (legacy format)".to_string())?);
        }
        let rest = line.sub_from(10).map_err(|_| "COMPND line too short".to_string())?;
        let field_list = java_split_ws(jtrim(&rest));
        if !field_list.is_empty() {
            let field0 = &field_list[0];
            if COMPND_FIELD_VALUES.contains(&field0.as_str()) {
                self.continuation_field = field0.clone();
                if self.previous_continuation_field.is_empty() {
                    self.previous_continuation_field = self.continuation_field.clone();
                }
            } else if field0.ends_with(';') && COMPND_FIELD_VALUES.contains(&&field0[..field0.len() - 1]) {
                return Ok(());
            }
        }
        let rest = if self.continuation_field.is_empty() {
            jtrim(&rest).to_string()
        } else {
            jtrim(&rest.replace(&self.continuation_field, "")).to_string()
        };
        for token in tokenize(&rest, TOKENIZER_DEFAULT_DELIMS) {
            if self.previous_continuation_field.is_empty() {
                self.previous_continuation_field = self.continuation_field.clone();
            }
            if self.previous_continuation_field == self.continuation_field
                && COMPND_FIELD_VALUES.contains(&self.continuation_field.as_str())
            {
                self.continuation_string.push_str(&token);
                self.continuation_string.push(' ');
            }
            if self.continuation_field != self.previous_continuation_field {
                if self.continuation_string.is_empty() {
                    self.continuation_string = token.clone();
                } else {
                    let f = self.previous_continuation_field.clone();
                    let v = self.continuation_string.clone();
                    self.compnd_value_setter(&f, &v);
                    self.previous_continuation_field = self.continuation_field.clone();
                    self.continuation_string = format!("{} ", token);
                }
            }
        }
        if self.is_last_compnd_line {
            let f = self.continuation_field.clone();
            let v = self.continuation_string.clone();
            self.compnd_value_setter(&f, &v);
            self.continuation_string = String::new();
            self.flush_pending_compound();
        }
        Ok(())
    }

    fn flush_pending_compound(&mut self) {
        if let Some(c) = self.pending_compound.take() {
            self.entities.push(c);
            self.current_compound = Some(self.entities.len() - 1);
        } else if let Some(ci) = self.current_compound {
            // re-adding an already added compound (BioJava adds the same object again)
            let mol_id = self.entities[ci].mol_id;
            self.entities.push(Compound { mol_id });
        }
    }

    fn current_mol_id(&self) -> Option<i32> {
        if let Some(c) = &self.pending_compound {
            return Some(c.mol_id);
        }
        self.current_compound.map(|i| self.entities[i].mol_id)
    }

    fn compnd_value_setter(&mut self, field: &str, value: &str) {
        let value = jtrim(value).replace(';', "");
        if field == "MOL_ID:" {
            let i = parse_java_int(jtrim(&value)).unwrap_or(-1);
            if i > 0 && self.prev_mol_id != i {
                if self.pending_compound.is_some() || self.current_compound.is_some() {
                    self.flush_pending_compound();
                }
                self.pending_compound = Some(Compound { mol_id: i });
                self.prev_mol_id = i;
            }
        }
        let mol_id = match self.current_mol_id() {
            Some(m) => m,
            None => return,
        };
        if field == "CHAIN:" {
            let mut chains = Vec::new();
            for t in tokenize(&value, &[',']) {
                let mut id = jtrim(&t).to_string();
                if id == "NULL" {
                    id = " ".to_string();
                }
                chains.push(id);
            }
            self.compound_chain_ids.insert(mol_id, chains);
        }
    }

    fn source_handler(&mut self, line0: &str) -> Result<(), String> {
        let mut line = JLine::new(line0.to_string());
        let err = |_| "SOURCE line too short".to_string();
        let continuation_nr = jtrim(&line.sub(9, 10).map_err(err)?).to_string();
        if line.len() > 79 {
            line = JLine::new(line.sub(0, 79).map_err(err)?);
        }
        let rest = line.sub_from(10).map_err(err)?;
        let field_list = java_split_ws(&rest);
        if field_list.is_empty() {
            return Err("ArrayIndexOutOfBoundsException in SOURCE handler".to_string());
        }
        if !field_list[0].is_empty() && SOURCE_FIELD_VALUES.contains(&field_list[0].as_str()) {
            self.continuation_field = field_list[0].clone();
            if self.previous_continuation_field.is_empty() {
                self.previous_continuation_field = self.continuation_field.clone();
            }
        } else if field_list.len() > 1 && SOURCE_FIELD_VALUES.contains(&field_list[1].as_str()) {
            self.continuation_field = field_list[1].clone();
            if self.previous_continuation_field.is_empty() {
                self.previous_continuation_field = self.continuation_field.clone();
            }
        } else if continuation_nr.is_empty() {
            self.continuation_field = "MOLECULE:".to_string();
            if self.previous_continuation_field.is_empty() {
                self.previous_continuation_field = self.continuation_field.clone();
            }
        }
        let rest = if self.continuation_field.is_empty() {
            jtrim(&rest).to_string()
        } else {
            jtrim(&rest.replace(&self.continuation_field, "")).to_string()
        };
        for token in tokenize(&rest, TOKENIZER_DEFAULT_DELIMS) {
            if self.previous_continuation_field.is_empty() {
                self.previous_continuation_field = self.continuation_field.clone();
            }
            if self.previous_continuation_field == self.continuation_field
                && SOURCE_FIELD_VALUES.contains(&self.continuation_field.as_str())
            {
                self.continuation_string.push_str(&token);
                self.continuation_string.push(' ');
            }
            if self.continuation_field != self.previous_continuation_field {
                if self.continuation_string.is_empty() {
                    self.continuation_string = token.clone();
                } else {
                    let f = self.previous_continuation_field.clone();
                    let v = self.continuation_string.clone();
                    self.source_value_setter(&f, &v)?;
                    self.previous_continuation_field = self.continuation_field.clone();
                    self.continuation_string = format!("{} ", token);
                }
            }
        }
        if self.is_last_source_line {
            let f = self.continuation_field.clone();
            let v = self.continuation_string.clone();
            self.source_value_setter(&f, &v)?;
            self.continuation_string = String::new();
        }
        Ok(())
    }

    fn source_value_setter(&mut self, field: &str, value: &str) -> Result<(), String> {
        let value = jtrim(value).replace(';', "");
        if field == "MOL_ID:" {
            match parse_java_int(jtrim(&value)) {
                None => return Ok(()),
                Some(i) => {
                    let idx = i as i64 - 1;
                    if idx < 0 || idx as usize >= self.entities.len() {
                        return Err(format!("IndexOutOfBoundsException in SOURCE MOL_ID {}", i));
                    }
                    self.current_compound = Some(idx as usize);
                }
            }
        }
        Ok(())
    }

    fn make_compounds(&mut self) -> Result<(), String> {
        let compnd = std::mem::take(&mut self.compnd_lines);
        let n = compnd.len();
        for line in compnd.iter() {
            let first_idx = compnd.iter().position(|l| l == line).unwrap();
            if first_idx + 1 == n {
                self.is_last_compnd_line = true;
            }
            self.compnd_handler(line)?;
        }
        if self.entities.is_empty() {
            // current_compound = new EntityInfo() (not added to entities)
            self.pending_compound = None;
            self.current_compound = None;
        } else {
            self.pending_compound = None;
            self.current_compound = Some(0);
        }
        let source = std::mem::take(&mut self.source_lines);
        let n = source.len();
        for line in source.iter() {
            let first_idx = source.iter().position(|l| l == line).unwrap();
            if first_idx + 1 == n {
                self.is_last_source_line = true;
            }
            self.source_handler(line)?;
        }
        Ok(())
    }

    // ---------------- end of file ----------------

    fn is_water_only(&self, chain: usize) -> bool {
        self.s.chains[chain].groups.iter().all(|&g| is_water_name(&self.s.groups[g].pdb_name))
    }

    fn is_polymeric(&mut self, g: usize) -> Result<bool, String> {
        let name = self.s.groups[g].pdb_name.clone();
        let cc = self.chem_comp(&name);
        let rt = match cc.residue_type {
            Some(r) => r,
            None => return Err(format!("unknown chem comp type for {}", name)),
        };
        if rt == ResidueType::NonPolymer {
            return Ok(false);
        }
        let pt = rt.polymer_type();
        Ok(is_protein_only(pt) || is_polynucleotide_only(pt) || rt == ResidueType::LPeptideLinking)
    }

    fn is_pure_non_polymer(&mut self, chain: usize) -> Result<bool, String> {
        let groups = self.s.chains[chain].groups.clone();
        for g in groups {
            if self.is_polymeric(g)? && !self.s.groups[g].het_in_file {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn split_non_poly_chain(&mut self, chain: usize) -> (Vec<usize>, Vec<usize>) {
        let mut split_non_polys = Vec::new();
        let mut water_chains = Vec::new();
        let mut split: Option<usize> = None;
        let mut previous_is_water = false;
        let name = self.s.chains[chain].name.clone();
        let groups = self.s.chains[chain].groups.clone();
        for g in groups {
            let is_water = is_water_name(&self.s.groups[g].pdb_name);
            if !previous_is_water {
                if let Some(sp) = split {
                    split_non_polys.push(sp);
                }
                split = Some(self.new_chain(&name));
            } else if !is_water {
                water_chains.push(split.unwrap());
                split = Some(self.new_chain(&name));
            }
            previous_is_water = is_water;
            self.add_group_to_chain(split.unwrap(), g);
        }
        if let Some(sp) = split {
            if previous_is_water {
                water_chains.push(sp);
            } else {
                split_non_polys.push(sp);
            }
        }
        (split_non_polys, water_chains)
    }

    fn assign_chains_and_entities(&mut self) -> Result<(), String> {
        let mut poly_models: Vec<Vec<usize>> = Vec::new();
        let mut non_poly_models: Vec<Vec<usize>> = Vec::new();
        let mut water_models: Vec<Vec<usize>> = Vec::new();
        let all_models = self.all_models.clone();
        for model in &all_models {
            let mut poly = Vec::new();
            let mut non_poly = Vec::new();
            let mut water = Vec::new();
            for &c in model {
                if self.is_water_only(c) {
                    water.push(c);
                } else if self.is_pure_non_polymer(c)? {
                    non_poly.push(c);
                } else {
                    poly.push(c);
                }
            }
            poly_models.push(poly);
            non_poly_models.push(non_poly);
            water_models.push(water);
        }

        let mut split_non_poly_models: Vec<Vec<usize>> = Vec::new();
        for i in 0..non_poly_models.len() {
            let mut split_non_polys = Vec::new();
            for &np in &non_poly_models[i].clone() {
                let (sn, sw) = self.split_non_poly_chain(np);
                split_non_polys.extend(sn);
                water_models[i].extend(sw);
            }
            split_non_poly_models.push(split_non_polys);
        }

        // asym ids
        for i in 0..poly_models.len() {
            let mut asym = "A".to_string();
            for &c in poly_models[i].iter().chain(split_non_poly_models[i].iter()).chain(water_models[i].iter()) {
                self.s.chains[c].id = asym.clone();
                asym = next_asym_id(&asym);
            }
        }

        if !self.entities.is_empty() {
            let entities: Vec<i32> = self.entities.iter().map(|e| e.mol_id).collect();
            for mol_id in entities {
                let chain_ids = match self.compound_chain_ids.get(&mol_id) {
                    Some(c) => c.clone(),
                    None => continue,
                };
                for cid in chain_ids {
                    for model in &poly_models {
                        for &c in model {
                            if self.s.chains[c].name == cid {
                                self.s.chains[c].entity = Some(EntityType::Polymer);
                            }
                        }
                    }
                }
            }
        } else if !poly_models.is_empty() {
            // EntityFinder.findPolyEntities: every chain of the first model gets a polymer entity,
            // chains of further models are matched by id (NPE in BioJava when there is no match)
            let first_ids: Vec<String> = poly_models[0].iter().map(|&c| self.s.chains[c].id.clone()).collect();
            for &c in &poly_models[0] {
                self.s.chains[c].entity = Some(EntityType::Polymer);
            }
            for model in poly_models.iter().skip(1) {
                for &c in model {
                    if !first_ids.contains(&self.s.chains[c].id) {
                        return Err(format!("EntityFinder: chain {} not found in first model", self.s.chains[c].id));
                    }
                    self.s.chains[c].entity = Some(EntityType::Polymer);
                }
            }
        }

        if !split_non_poly_models.is_empty() {
            for model in &split_non_poly_models {
                for &c in model {
                    self.s.chains[c].entity = Some(EntityType::NonPolymer);
                }
            }
            for model in &water_models {
                for &c in model {
                    self.s.chains[c].entity = Some(EntityType::Water);
                }
            }
        }

        for i in 0..all_models.len() {
            let mut model = Vec::new();
            model.extend(poly_models[i].iter());
            model.extend(split_non_poly_models[i].iter());
            model.extend(water_models[i].iter());
            self.s.models.push(model);
        }
        Ok(())
    }

    fn trigger_end_file_checks(&mut self) -> Result<(), String> {
        if let (Some(c), Some(g)) = (self.current_chain, self.current_group) {
            self.add_group_to_chain(c, g);
        }
        if let (Some(_), Some(c)) = (self.current_model.as_ref(), self.current_chain) {
            self.current_model.as_mut().unwrap().push(c);
        }
        if let Some(m) = self.current_model.take() {
            self.all_models.push(m);
        }
        self.assign_chains_and_entities()
    }
}

fn next_char(c: char) -> char {
    if c != 'Z' {
        char::from_u32(c as u32 + 1).unwrap_or('?')
    } else {
        'A'
    }
}

/// PDBFileParser.getNextAsymId
fn next_asym_id(id: &str) -> String {
    let c: Vec<char> = id.chars().collect();
    match c.len() {
        1 => {
            if id != "Z" {
                next_char(c[0]).to_string()
            } else {
                "AA".to_string()
            }
        }
        2 => {
            if id == "ZZ" {
                return "AAA".to_string();
            }
            let mut c = c;
            c[0] = next_char(c[0]);
            if c[0] == 'A' {
                c[1] = next_char(c[1]);
            }
            c.into_iter().collect()
        }
        3 => {
            let mut c = c;
            c[0] = next_char(c[0]);
            if c[0] == 'A' {
                c[1] = next_char(c[1]);
                if c[1] == 'A' {
                    c[2] = next_char(c[2]);
                }
            }
            c.into_iter().collect()
        }
        _ => "null".to_string(),
    }
}

/// Parse PDB file content (already decompressed).
pub fn parse_pdb(bytes: &[u8], dict: &ChemCompDictionary) -> Result<Structure, ParseError> {
    let mut p = Parser {
        dict,
        cc_cache: HashMap::new(),
        s: Structure::default(),
        all_models: Vec::new(),
        current_model: None,
        current_chain: None,
        current_group: None,
        start_of_molecule: true,
        start_of_model: true,
        atom_count: 0,
        is_legacy_format: false,
        pdb_id: None,
        compnd_lines: Vec::new(),
        source_lines: Vec::new(),
        entities: Vec::new(),
        compound_chain_ids: HashMap::new(),
        current_compound: None,
        pending_compound: None,
        prev_mol_id: -1,
        previous_continuation_field: String::new(),
        continuation_field: String::new(),
        continuation_string: String::new(),
        is_last_compnd_line: false,
        is_last_source_line: false,
    };

    for line in crate::jutil::read_lines(bytes) {
        if line.is_empty() {
            continue;
        }
        if line.starts_with("END") {
            continue;
        }
        let jl = JLine::new(line);
        if jl.len() < 6 && !jl.starts_with("TER") {
            continue;
        }
        let record_name = if jl.len() < 6 {
            jtrim(jl.as_str()).to_string()
        } else {
            jtrim(&jl.sub(0, 6).unwrap()).to_string()
        };
        let res = match record_name.as_str() {
            "ATOM" | "HETATM" => p.atom_handler(&jl),
            "MODEL" => {
                p.start_of_molecule = true;
                p.start_of_model = true;
                Ok(())
            }
            "TER" => {
                p.start_of_molecule = true;
                Ok(())
            }
            "HEADER" => p.header_handler(&jl),
            "SOURCE" => {
                p.source_lines.push(jl.as_str().to_string());
                Ok(())
            }
            "COMPND" => {
                p.compnd_lines.push(jl.as_str().to_string());
                Ok(())
            }
            _ => Ok(()),
        };
        match res {
            Ok(()) | Err(HErr::Caught) => {}
            Err(HErr::Fatal(msg)) => return Err(ParseError(msg)),
        }
    }

    p.make_compounds().map_err(ParseError)?;
    p.trigger_end_file_checks().map_err(ParseError)?;
    Ok(p.s)
}

impl Structure {
    pub fn atom_resname<'a>(&'a self, atom: usize) -> &'a str {
        &self.groups[self.atoms[atom].group].pdb_name
    }

    pub fn group_has_atom(&self, g: usize, name: &str) -> bool {
        self.groups[g].atoms.iter().any(|&a| self.atoms[a].name == name)
    }

    /// Group.hasAminoAtoms
    pub fn has_amino_atoms(&self, g: usize) -> bool {
        self.group_has_atom(g, "CA") && self.group_has_atom(g, "C") && self.group_has_atom(g, "N") && self.group_has_atom(g, "O")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_split_ws() {
        assert_eq!(java_split_ws("a b"), vec!["a", "b"]);
        assert_eq!(java_split_ws(" a b"), vec!["", "a", "b"]);
        assert_eq!(java_split_ws("a b  "), vec!["a", "b"]);
        assert_eq!(java_split_ws(""), vec![""]);
        assert!(java_split_ws("   ").is_empty());
    }

    #[test]
    fn test_asym() {
        assert_eq!(next_asym_id("A"), "B");
        assert_eq!(next_asym_id("Z"), "AA");
        assert_eq!(next_asym_id("AA"), "BA");
        assert_eq!(next_asym_id("ZA"), "AB");
    }
}
