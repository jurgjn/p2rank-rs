//! P2Rank `Protein` loading: residue chains, residues, protein heavy atoms.

use crate::atoms::{PointSet, KD_TREE_THRESHOLD};
use crate::elements::ELEMENT_H;
use crate::pdb::{EntityType, GroupKind, ResNum, Structure};

pub struct Residue {
    pub group: usize,
    pub chain_author_id: String,
}

#[allow(dead_code)] // mirrors the Java data model
pub struct Protein {
    pub name: String,
    pub s: Structure,
    /// chains of the (model 0, alt-conformation reduced) structure
    pub chains: Vec<usize>,
    pub residues: Vec<Residue>,
    /// protein heavy atoms (indices into s.atoms)
    pub protein_atoms: Vec<usize>,
}

/// Struct.isHydrogenAtom
pub fn is_hydrogen_atom(s: &Structure, a: usize) -> bool {
    let atom = &s.atoms[a];
    if atom.element == ELEMENT_H {
        return true;
    }
    if atom.name.starts_with('H') {
        return true;
    }
    let mut it = atom.name.encode_utf16();
    it.next();
    if let Some(c) = it.next() {
        if c == 'H' as u16 {
            return true;
        }
    }
    false
}

fn is_amino_acid_group(s: &Structure, g: Option<usize>) -> bool {
    match g {
        None => false,
        Some(g) => s.groups[g].kind == GroupKind::AminoAcid,
    }
}

/// Struct.isAminoAcidResidue
fn is_amino_acid_residue(s: &Structure, g: usize) -> bool {
    if s.groups[g].kind == GroupKind::AminoAcid {
        return true;
    }
    s.groups[g].pdb_name.starts_with("UNK")
}

/// Struct.getResidueGroupsFromChain
fn residue_groups_from_chain(s: &Structure, chain: usize) -> Vec<usize> {
    let groups = &s.chains[chain].groups;
    let n = groups.len();
    let mut res = Vec::new();
    for i in 0..n {
        let g = groups[i];
        let prev = if i > 0 { Some(groups[i - 1]) } else { None };
        let next = if i + 1 < n { Some(groups[i + 1]) } else { None };
        let is_res = is_amino_acid_residue(s, g)
            || (s.has_amino_atoms(g) && (prev.is_none() || is_amino_acid_group(s, prev)) && (next.is_some() && is_amino_acid_group(s, next)));
        if is_res {
            res.push(g);
        }
    }
    res
}

fn chain_atoms(s: &Structure, chain: usize) -> Vec<usize> {
    let mut res = Vec::new();
    for &g in &s.chains[chain].groups {
        res.extend(s.groups[g].atoms.iter().copied());
    }
    res
}

fn point_set(s: &Structure, atoms: &[usize]) -> PointSet {
    PointSet::new(atoms.iter().map(|&a| [s.atoms[a].x, s.atoms[a].y, s.atoms[a].z]).collect(), atoms.iter().map(|&a| a as u32).collect())
}

/// AlternateChainReducer.uniformAltLoc
fn uniform_alt_loc(s: &Structure, atoms: &[usize]) -> Option<char> {
    let total = atoms.len();
    if total == 0 {
        return None;
    }
    let mut letter: Option<char> = None;
    let mut nonblank = 0usize;
    for &a in atoms {
        let al = s.atoms[a].alt_loc;
        if al == ' ' {
            continue;
        }
        match letter {
            None => letter = Some(al),
            Some(l) if l != al => return None,
            _ => {}
        }
        nonblank += 1;
    }
    let letter = letter?;
    if (nonblank as f64) >= 0.5 * (total as f64) {
        Some(letter)
    } else {
        None
    }
}

fn overlaps(s: &Structure, cand: &[usize], refs: &[&Vec<usize>]) -> bool {
    if refs.is_empty() || cand.is_empty() {
        return false;
    }
    let mut all = Vec::new();
    for r in refs {
        all.extend(r.iter().copied());
    }
    if all.is_empty() {
        return false;
    }
    let ps = point_set(s, &all);
    if ps.len() > KD_TREE_THRESHOLD {
        ps.kd();
    }
    let mut within = 0usize;
    let need = (0.7f64 * cand.len() as f64).ceil() as usize;
    for &a in cand {
        let p = [s.atoms[a].x, s.atoms[a].y, s.atoms[a].z];
        if ps.are_within_distance(&p, 2.0) {
            within += 1;
            if within >= need {
                return true;
            }
        }
    }
    false
}

/// AlternateChainReducer.reduceAlternateConformationChains
fn reduce_alternate_conformation_chains(s: &Structure, chains: Vec<usize>) -> Vec<usize> {
    if chains.len() < 2 {
        return chains;
    }
    struct Info {
        chain: usize,
        atoms: Vec<usize>,
        alt_loc: Option<char>,
    }
    let infos: Vec<Info> = chains
        .iter()
        .map(|&c| {
            let atoms = chain_atoms(s, c);
            let alt_loc = uniform_alt_loc(s, &atoms);
            Info { chain: c, atoms, alt_loc }
        })
        .collect();
    let mut references: Vec<&Info> = Vec::new();
    let mut candidates: Vec<&Info> = Vec::new();
    for info in &infos {
        if info.alt_loc.is_some() {
            candidates.push(info);
        } else {
            references.push(info);
        }
    }
    let chain_id = |c: usize| -> String {
        let id = &s.chains[c].id;
        if id.is_empty() {
            s.chains[c].name.clone()
        } else {
            id.clone()
        }
    };
    candidates.sort_by(|a, b| {
        let c = (a.alt_loc.unwrap() as u32).cmp(&(b.alt_loc.unwrap() as u32));
        if c != std::cmp::Ordering::Equal {
            c
        } else {
            crate::jutil::java_str_cmp(&chain_id(a.chain), &chain_id(b.chain))
        }
    });
    let mut kept: Vec<&Info> = Vec::new();
    let mut removed: Vec<usize> = Vec::new();
    for cand in candidates {
        let mut refs: Vec<&Vec<usize>> = references.iter().map(|r| &r.atoms).collect();
        for r in &kept {
            if (r.alt_loc.unwrap() as u32) < (cand.alt_loc.unwrap() as u32) {
                refs.push(&r.atoms);
            }
        }
        if overlaps(s, &cand.atoms, &refs) {
            removed.push(cand.chain);
        } else {
            kept.push(cand);
        }
    }
    if removed.is_empty() {
        return chains;
    }
    chains.into_iter().filter(|c| !removed.contains(c)).collect()
}

impl Protein {
    pub fn from_structure(name: String, s: Structure, only_chains: Option<&[String]>) -> Result<Protein, String> {
        // model 0 (Struct.reduceStructureToModel0 / Structure.getChains)
        let mut chains: Vec<usize> = s.models.first().cloned().unwrap_or_default();
        if let Some(only) = only_chains {
            // PdbUtils.reduceStructureToChains: keep chains whose author id is listed
            chains.retain(|&c| only.contains(&s.chains[c].name));
        }
        let chains = reduce_alternate_conformation_chains(&s, chains);

        let mut residues = Vec::new();
        for &c in &chains {
            if s.chains[c].entity != Some(EntityType::Polymer) {
                continue;
            }
            let groups = residue_groups_from_chain(&s, c);
            if groups.is_empty() {
                continue;
            }
            for g in groups {
                let chain_author_id = match s.groups[g].chain {
                    Some(gc) => s.chains[gc].name.clone(),
                    None => "null".to_string(),
                };
                residues.push(Residue { group: g, chain_author_id });
            }
        }

        let mut protein_atoms = Vec::new();
        for r in &residues {
            for &a in &s.groups[r.group].atoms {
                if !is_hydrogen_atom(&s, a) {
                    protein_atoms.push(a);
                }
            }
        }
        Ok(Protein { name, s, chains, residues, protein_atoms })
    }

    pub fn atom_pos(&self, a: usize) -> [f64; 3] {
        let at = &self.s.atoms[a];
        [at.x, at.y, at.z]
    }

    pub fn point_set(&self, atoms: &[usize]) -> PointSet {
        point_set(&self.s, atoms)
    }

    pub fn residue_key(&self, r: &Residue) -> &ResNum {
        &self.s.groups[r.group].resnum
    }

    /// Residue.getCode(): PDB name of the group of the first atom of the residue's group
    pub fn residue_code(&self, r: &Residue) -> Option<&str> {
        let g = &self.s.groups[r.group];
        g.atoms.first().map(|&a| self.s.groups[self.s.atoms[a].group].pdb_name.as_str())
    }

    /// heavy atoms of a residue (Residue.getAtoms())
    pub fn residue_heavy_atoms(&self, r: &Residue) -> Vec<usize> {
        self.s.groups[r.group].atoms.iter().copied().filter(|&a| !is_hydrogen_atom(&self.s, a)).collect()
    }
}
