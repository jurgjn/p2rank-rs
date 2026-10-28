//! mmCIF (text) reader reproducing ciftools-java's tokenizer and BioJava's
//! `CifStructureConsumerImpl`, limited to the data P2Rank uses (atoms, groups, chains,
//! models, entity types).

use crate::chemcomp::{is_polynucleotide_only, is_protein_only, ChemComp, ChemCompDictionary};
use crate::elements::{element_value_of_ignore_case, ELEMENT_R};
use crate::jutil::{jtrim, parse_java_double, parse_java_int};
use crate::pdb::{is_water_name, Atom, Chain, EntityType, Group, GroupKind, ParseError, ResNum, Structure};
use std::collections::HashMap;

// ------------------------------------------------------------------ tokenizer (ciftools TokenizerState)

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Tok {
    Data,
    Save,
    Loop,
    Value,
    ColumnName,
    Comment,
    End,
}

struct Tokenizer<'a> {
    data: &'a [u8],
    position: usize,
    is_escaped: bool,
    in_save_frame: bool,
    token: Tok,
    start: usize,
    end: usize,
}

/// A text column: token ranges into the data
#[derive(Clone, Default)]
pub struct Column {
    ranges: Vec<(usize, usize)>,
}

pub struct Category {
    columns: Vec<(String, Column)>,
    row_count: usize,
}

type Categories = Vec<(String, Category)>;

impl<'a> Tokenizer<'a> {
    fn ch(&self, i: usize) -> u8 {
        self.data[i]
    }

    fn eat_value(&mut self) {
        while self.position < self.data.len() {
            match self.ch(self.position) {
                b'\t' | b'\n' | b'\r' | b' ' => {
                    self.end = self.position;
                    return;
                }
                _ => self.position += 1,
            }
        }
        self.end = self.position;
    }

    fn eat_escaped(&mut self, esc: u8) {
        self.position += 1;
        while self.position < self.data.len() {
            let c = self.ch(self.position);
            if c == esc {
                if self.position + 1 >= self.data.len() {
                    self.start += 1;
                    self.end = self.position;
                    self.is_escaped = true;
                    self.position += 1;
                    return;
                }
                let next = self.ch(self.position + 1);
                if next == b'\t' || next == b'\n' || next == b'\r' || next == b' ' {
                    self.start += 1;
                    self.end = self.position;
                    self.is_escaped = true;
                    self.position += 1;
                    return;
                }
            } else if c == b'\r' || c == b'\n' {
                self.end = self.position;
                return;
            }
            self.position += 1;
        }
        self.end = self.position;
    }

    fn is_triple_quote_at_position(&self) -> bool {
        if self.data.len() - self.position < 2 {
            return false;
        }
        // note: ciftools reads position+2 without a bounds check beyond this test
        if self.ch(self.position + 1) != b'\'' {
            return false;
        }
        self.position + 2 < self.data.len() && self.ch(self.position + 2) == b'\''
    }

    fn eat_triple_quote(&mut self) {
        self.position += 3;
        while self.position < self.data.len() {
            if self.ch(self.position) == b'\'' && self.is_triple_quote_at_position() {
                self.start += 3;
                self.end = self.position;
                self.is_escaped = true;
                self.position += 3;
                return;
            }
            self.position += 1;
        }
        self.end = self.position;
    }

    fn eat_multiline(&mut self) {
        let mut prev = b';';
        let mut pos = self.position + 1;
        while pos < self.data.len() {
            let c = self.ch(pos);
            if c == b';' && (prev == b'\n' || prev == b'\r') {
                self.position = pos + 1;
                self.start += 1;
                pos -= 1;
                let mut c2 = self.ch(pos);
                while c2 == b'\n' || c2 == b'\r' {
                    pos -= 1;
                    c2 = self.ch(pos);
                }
                self.end = pos + 1;
                self.is_escaped = true;
                return;
            }
            prev = c;
            pos += 1;
        }
        self.position = pos;
    }

    fn skip_comment_line(&mut self) {
        while self.position < self.data.len() {
            let c = self.ch(self.position);
            if c == b'\r' || c == b'\n' {
                return;
            }
            self.position += 1;
        }
    }

    fn skip_whitespace(&mut self) -> u8 {
        let mut prev = b'\n';
        while self.position < self.data.len() {
            let c = self.ch(self.position);
            match c {
                b'\t' | b' ' | b'\n' | b'\r' => {
                    prev = c;
                    self.position += 1;
                }
                _ => return prev,
            }
        }
        prev
    }

    fn prefix_ci(&self, p: &[u8]) -> bool {
        self.start + p.len() <= self.data.len() && self.data[self.start..self.start + p.len()].eq_ignore_ascii_case(p)
    }

    fn move_next_internal(&mut self) {
        let prev = self.skip_whitespace();
        if self.position >= self.data.len() {
            self.token = Tok::End;
            return;
        }
        self.start = self.position;
        self.end = self.position;
        self.is_escaped = false;
        let c = self.ch(self.position);
        match c {
            b'"' | b'\'' => {
                if c == b'\'' && self.is_triple_quote_at_position() {
                    self.eat_triple_quote();
                } else {
                    self.eat_escaped(c);
                }
                self.token = Tok::Value;
            }
            b'#' => {
                self.skip_comment_line();
                self.token = Tok::Comment;
            }
            b';' => {
                if prev != b'\n' && prev != b'\r' {
                    self.eat_value();
                } else {
                    self.eat_multiline();
                }
                self.token = Tok::Value;
            }
            _ => {
                self.eat_value();
                if self.is_escaped {
                    self.token = Tok::Value;
                } else if self.ch(self.start) == b'_' {
                    self.token = Tok::ColumnName;
                } else if self.end - self.start < 5 || self.ch(self.start + 4) != b'_' {
                    self.token = Tok::Value;
                } else if self.prefix_ci(b"data") {
                    self.token = Tok::Data;
                } else if self.prefix_ci(b"save") {
                    self.token = Tok::Save;
                } else if self.prefix_ci(b"loop") {
                    self.token = Tok::Loop;
                } else {
                    self.token = Tok::Value;
                }
            }
        }
    }

    fn move_next(&mut self) {
        self.move_next_internal();
        while self.token == Tok::Comment {
            self.move_next_internal();
        }
    }

    fn token_str(&self) -> String {
        String::from_utf8_lossy(&self.data[self.start..self.end]).to_string()
    }

    fn namespace_end(&self) -> usize {
        match self.data[self.start..self.end].iter().position(|&c| c == b'.') {
            Some(i) => self.start + i,
            None => self.end,
        }
    }

    fn is_namespace(&self, ns_start: usize, ns_end: usize) -> bool {
        let ns_len = ns_end - ns_start;
        let token_len = self.end - self.start;
        if token_len < ns_len {
            return false;
        }
        if self.data[ns_start..ns_end] != self.data[self.start..self.start + ns_len] {
            return false;
        }
        if ns_len == token_len {
            return true;
        }
        self.ch(self.start + ns_len) == b'.'
    }

    fn is_flat_namespace(&self) -> bool {
        !self.data[self.start..self.end].contains(&b'.')
    }

    fn handle_single(&mut self, ctx: &mut Categories) -> Result<(), String> {
        let ns_start = self.start;
        let ns_end = self.namespace_end();
        let name = String::from_utf8_lossy(&self.data[self.start..ns_end]).to_string();
        let is_flat = self.is_flat_namespace();
        let mut fields: Vec<(String, Column)> = Vec::new();
        let category_name = name[1..].to_string();
        while self.token == Tok::ColumnName && self.is_namespace(ns_start, ns_end) {
            let column_name = if is_flat { String::new() } else { self.token_str()[name.len() + 1..].to_string() };
            self.move_next();
            if self.token != Tok::Value {
                return Err("Expected value.".into());
            }
            let col = Column { ranges: vec![(self.start, self.end)] };
            put_ci(&mut fields, column_name, col);
            self.move_next();
        }
        let rc = fields.first().map(|c| c.1.ranges.len()).unwrap_or(0);
        put_ci(&mut *ctx, category_name, Category { columns: fields, row_count: rc });
        Ok(())
    }

    fn handle_loop(&mut self, ctx: &mut Categories) -> Result<(), String> {
        self.move_next();
        let ns_end = self.namespace_end();
        let name = String::from_utf8_lossy(&self.data[self.start..ns_end]).to_string();
        let is_flat = self.is_flat_namespace();
        let mut column_names: Vec<String> = Vec::new();
        let mut cols: Vec<Column> = Vec::new();
        while self.token == Tok::ColumnName {
            let cn = if is_flat { self.token_str() } else { self.token_str()[name.len() + 1..].to_string() };
            column_names.push(cn);
            cols.push(Column::default());
            self.move_next();
        }
        let mut token_count = 0usize;
        while self.token == Tok::Value {
            if column_names.is_empty() {
                return Err("ArithmeticException: / by zero".into());
            }
            let i = token_count % column_names.len();
            cols[i].ranges.push((self.start, self.end));
            self.move_next();
            token_count += 1;
        }
        if is_flat {
            for (i, c) in cols.into_iter().enumerate() {
                let flat = column_names[i][1..].to_string();
                let rc = c.ranges.len();
                put_ci(&mut *ctx, flat, Category { columns: vec![(String::new(), c)], row_count: rc });
            }
        } else {
            let category_name = name[1..].to_string();
            let mut columns: Vec<(String, Column)> = Vec::new();
            for (i, c) in cols.into_iter().enumerate() {
                put_ci(&mut columns, column_names[i].clone(), c);
            }
            let rc = columns.first().map(|c| c.1.ranges.len()).unwrap_or(0);
            put_ci(&mut *ctx, category_name, Category { columns, row_count: rc });
        }
        Ok(())
    }
}

/// LinkedCaseInsensitiveMap.put: replaces the value of an existing (case-insensitive) key in place
fn put_ci<T>(v: &mut Vec<(String, T)>, k: String, val: T) {
    if let Some(e) = v.iter_mut().find(|e| e.0.eq_ignore_ascii_case(&k)) {
        e.1 = val;
    } else {
        v.push((k, val));
    }
}

/// TextCifReader.readText: returns the categories of the first data block
fn read_first_block(data: &[u8]) -> Result<Categories, String> {
    if data.is_empty() {
        return Err("ParsingException: Cannot parse empty file.".into());
    }
    let mut t = Tokenizer { data, position: 0, is_escaped: false, in_save_frame: false, token: Tok::End, start: 0, end: 0 };
    let mut blocks: Vec<Categories> = Vec::new();
    let mut block_ctx: Categories = Vec::new();
    let mut save_ctx: Categories = Vec::new();
    t.move_next();
    while t.token != Tok::End {
        match t.token {
            Tok::Data => {
                if t.in_save_frame {
                    return Err("ParsingException: Unexpected data block inside a save frame.".into());
                }
                if !block_ctx.is_empty() {
                    blocks.push(std::mem::take(&mut block_ctx));
                }
                t.move_next();
            }
            Tok::Save => {
                let header_empty = t.end - t.start <= 5;
                if header_empty {
                    t.in_save_frame = false;
                } else {
                    if t.in_save_frame {
                        return Err("ParsingException: Save frames cannot be nested.".into());
                    }
                    t.in_save_frame = true;
                    save_ctx = Vec::new();
                }
                t.move_next();
            }
            Tok::Loop => {
                if t.in_save_frame {
                    t.handle_loop(&mut save_ctx)?;
                } else {
                    t.handle_loop(&mut block_ctx)?;
                }
            }
            Tok::ColumnName => {
                if t.in_save_frame {
                    t.handle_single(&mut save_ctx)?;
                } else {
                    t.handle_single(&mut block_ctx)?;
                }
            }
            tok => return Err(format!("ParsingException: Unexpected token ({:?}). Expected data_, loop_, or data name.", tok)),
        }
    }
    if t.in_save_frame {
        return Err("ParsingException: Unfinished save frame.".into());
    }
    if !block_ctx.is_empty() {
        blocks.push(block_ctx);
    }
    if blocks.is_empty() {
        return Err("IndexOutOfBoundsException: Index 0 out of bounds for length 0".into());
    }
    Ok(blocks.swap_remove(0))
}

// ------------------------------------------------------------------ typed column access

struct Col<'a> {
    data: &'a [u8],
    col: Option<&'a Column>,
    /// decoded BinaryCIF column data
    bin: Option<crate::bcif::Data>,
    name: &'a str,
}

/// Java Double.toString is only needed for float columns read as strings; not supported
fn bin_unsupported(name: &str) -> String {
    format!("UnsupportedOperationException: float BinaryCIF column {} read as string", name)
}

impl<'a> Col<'a> {
    fn defined(&self) -> bool {
        if let Some(b) = &self.bin {
            return match b {
                crate::bcif::Data::Ints(v, _) => !v.is_empty(),
                crate::bcif::Data::Floats(v, _) => !v.is_empty(),
                crate::bcif::Data::Strs(v) => !v.is_empty(),
                crate::bcif::Data::Bytes(v) => !v.is_empty(),
            };
        }
        self.col.map(|c| !c.ranges.is_empty()).unwrap_or(false)
    }
    fn raw(&self, row: usize) -> Result<&'a str, String> {
        let c = self.col.ok_or_else(|| format!("EmptyColumnException: column {} is undefined", self.name))?;
        let &(s, e) = c.ranges.get(row).ok_or_else(|| format!("ArrayIndexOutOfBoundsException: Index {} out of bounds", row))?;
        std::str::from_utf8(&self.data[s..e]).map_err(|_| "invalid utf-8 in cif".to_string())
    }
    fn oob(row: usize) -> String {
        format!("ArrayIndexOutOfBoundsException: Index {} out of bounds", row)
    }
    /// StrColumn.get (text: '.' and '?' read as empty; binary: raw value)
    fn str(&self, row: usize) -> Result<String, String> {
        if let Some(b) = &self.bin {
            return match b {
                crate::bcif::Data::Strs(v) => v.get(row).cloned().ok_or_else(|| Self::oob(row)),
                crate::bcif::Data::Ints(v, _) => v.get(row).map(|i| i.to_string()).ok_or_else(|| Self::oob(row)),
                _ => Err(bin_unsupported(self.name)),
            };
        }
        if self.col.is_none() && self.name.is_empty() {
            return Err("EmptyColumnException".into());
        }
        let v = self.raw(row)?;
        Ok(if v == "." || v == "?" { String::new() } else { v.to_string() })
    }
    /// Column.getStringData ('.' and '?' read as empty)
    fn string_data(&self, row: usize) -> Result<String, String> {
        if let Some(crate::bcif::Data::Strs(v)) = &self.bin {
            let x = v.get(row).ok_or_else(|| Self::oob(row))?;
            return Ok(if x == "." || x == "?" { String::new() } else { x.clone() });
        }
        self.str(row)
    }
    fn parse_int_text(t: &str) -> Result<i32, String> {
        if t.is_empty() || t == "." || t == "?" {
            return Ok(0);
        }
        let t = match t.find('.') {
            Some(i) => &t[..i],
            None => t,
        };
        parse_java_int(t).ok_or_else(|| format!("NumberFormatException: For input string: \"{}\"", t))
    }
    fn parse_float_text(t: &str) -> Result<f64, String> {
        if t.is_empty() || t == "." || t == "?" {
            return Ok(0.0);
        }
        let t = match t.find('(') {
            Some(i) => &t[..i],
            None => t,
        };
        parse_java_double(t).ok_or_else(|| {
            if jtrim(t).is_empty() {
                "NumberFormatException: empty String".to_string()
            } else {
                format!("NumberFormatException: For input string: \"{}\"", t)
            }
        })
    }
    /// IntColumn.get
    fn int(&self, row: usize) -> Result<i32, String> {
        if let Some(b) = &self.bin {
            return match b {
                crate::bcif::Data::Ints(v, _) => v.get(row).copied().ok_or_else(|| Self::oob(row)),
                crate::bcif::Data::Floats(v, _) => v.get(row).map(|&d| d as i32).ok_or_else(|| Self::oob(row)),
                crate::bcif::Data::Strs(v) => Self::parse_int_text(v.get(row).ok_or_else(|| Self::oob(row))?),
                _ => Err("ClassCastException".into()),
            };
        }
        let t = self.str(row)?;
        Self::parse_int_text(&t)
    }
    /// FloatColumn.get
    fn float(&self, row: usize) -> Result<f64, String> {
        if let Some(b) = &self.bin {
            return match b {
                crate::bcif::Data::Floats(v, _) => v.get(row).copied().ok_or_else(|| Self::oob(row)),
                crate::bcif::Data::Ints(v, _) => v.get(row).map(|&i| i as f64).ok_or_else(|| Self::oob(row)),
                crate::bcif::Data::Strs(v) => Self::parse_float_text(v.get(row).ok_or_else(|| Self::oob(row))?),
                _ => Err("ClassCastException".into()),
            };
        }
        let t = self.str(row)?;
        Self::parse_float_text(&t)
    }
}

enum Cats {
    Text(Categories),
    Bin(Vec<(String, crate::bcif::BCategory)>),
}

struct Block<'a> {
    data: &'a [u8],
    cats: Cats,
}

impl<'a> Block<'a> {
    fn row_count(&self, cat: &str) -> usize {
        match &self.cats {
            Cats::Text(c) => c.iter().find(|c| c.0.eq_ignore_ascii_case(cat)).map(|c| c.1.row_count).unwrap_or(0),
            Cats::Bin(c) => c.iter().find(|c| c.0.eq_ignore_ascii_case(cat)).map(|c| c.1.row_count).unwrap_or(0),
        }
    }
    fn col<'b>(&'b self, cat: &str, name: &'b str) -> Result<Col<'b>, String> {
        match &self.cats {
            Cats::Text(c) => {
                let col = c
                    .iter()
                    .find(|c| c.0.eq_ignore_ascii_case(cat))
                    .and_then(|c| c.1.columns.iter().find(|x| x.0.eq_ignore_ascii_case(name)).map(|x| &x.1));
                Ok(Col { data: self.data, col, bin: None, name })
            }
            Cats::Bin(c) => {
                let enc = c
                    .iter()
                    .find(|c| c.0.eq_ignore_ascii_case(cat))
                    .and_then(|c| c.1.columns.iter().find(|x| x.0.eq_ignore_ascii_case(name)).map(|x| &x.1));
                let bin = match enc {
                    Some(e) => Some(crate::bcif::decode_column_data(e.get("data").ok_or("missing column data")?)?),
                    None => None,
                };
                Ok(Col { data: self.data, col: None, bin, name })
            }
        }
    }
}

// ------------------------------------------------------------------ structure building

fn amino_1letter(code: &str) -> Option<char> {
    crate::pdb::amino_acid_1letter_pub(code)
}

struct Builder<'d> {
    dict: &'d ChemCompDictionary,
    cc_cache: HashMap<String, ChemComp>,
    s: Structure,
}

impl<'d> Builder<'d> {
    fn chem_comp(&mut self, name: &str) -> &ChemComp {
        let key = jtrim(&name.to_uppercase()).to_string();
        if !self.cc_cache.contains_key(&key) {
            let cc = self.dict.get(&key);
            self.cc_cache.insert(key.clone(), cc);
        }
        self.cc_cache.get(&key).unwrap()
    }

    /// CifStructureConsumerImpl.createGroup
    fn create_group(&mut self, record: &str, one_letter: Option<char>, code3: &str) -> GroupKind {
        let cc = self.chem_comp(code3);
        if !cc.is_empty() {
            return if is_protein_only(cc.polymer_type) {
                GroupKind::AminoAcid
            } else if is_polynucleotide_only(cc.polymer_type) {
                GroupKind::Nucleotide
            } else {
                GroupKind::Hetatm
            };
        }
        let nuc = crate::pdb::is_nucleotide_pub(code3);
        if record == "ATOM" {
            if nuc {
                GroupKind::Nucleotide
            } else if one_letter.is_none() || one_letter == Some('X') {
                GroupKind::Hetatm
            } else {
                GroupKind::AminoAcid
            }
        } else if nuc {
            GroupKind::Nucleotide
        } else if one_letter.is_some() {
            GroupKind::AminoAcid
        } else {
            GroupKind::Hetatm
        }
    }

    fn new_group(&mut self, record: &str, one_letter: Option<char>, code3: &str, resnum: ResNum, het: bool) -> usize {
        let kind = self.create_group(record, one_letter, code3);
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

    fn new_chain(&mut self, name: &str, id: &str) -> usize {
        self.s.chains.push(Chain { name: name.to_string(), id: id.to_string(), groups: Vec::new(), entity: None });
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

    fn alt_loc_group(&mut self, cg: usize, record: &str, alt_loc: char, one_letter: Option<char>, code3: &str) -> usize {
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
        let id = self.new_group(record, one_letter, code3, resnum, false);
        self.s.groups[cg].alt_locs.push(id);
        id
    }

    fn is_polymeric(&mut self, g: usize) -> Result<bool, String> {
        let name = self.s.groups[g].pdb_name.clone();
        let cc = self.chem_comp(&name);
        let rt = cc.residue_type.ok_or_else(|| format!("NullPointerException: unknown chem comp type for {}", name))?;
        if rt == crate::chemcomp::ResidueType::NonPolymer {
            return Ok(false);
        }
        let pt = rt.polymer_type();
        Ok(is_protein_only(pt) || is_polynucleotide_only(pt) || rt == crate::chemcomp::ResidueType::LPeptideLinking)
    }

    fn is_water_only(&self, c: usize) -> bool {
        self.s.chains[c].groups.iter().all(|&g| is_water_name(&self.s.groups[g].pdb_name))
    }

    fn is_pure_non_polymer(&mut self, c: usize) -> Result<bool, String> {
        for g in self.s.chains[c].groups.clone() {
            if self.is_polymeric(g)? && !self.s.groups[g].het_in_file {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

/// Parses an mmCIF text file into a Structure following BioJava semantics.
pub fn parse_cif(bytes: &[u8], dict: &ChemCompDictionary) -> Result<Structure, ParseError> {
    parse_cif_inner(bytes, dict).map_err(ParseError)
}

fn parse_cif_inner(bytes: &[u8], dict: &ChemCompDictionary) -> Result<Structure, String> {
    let text;
    let blk = if bytes.first() == Some(&0x83) {
        Block { data: bytes, cats: Cats::Bin(crate::bcif::read_first_block(bytes)?) }
    } else {
        // BufferedReader.lines() joined with "\n": normalizes line terminators, drops the last one
        text = crate::jutil::read_lines(bytes).join("\n");
        let data = text.as_bytes();
        Block { data, cats: Cats::Text(read_first_block(data)?) }
    };
    let mut b = Builder { dict, cc_cache: HashMap::new(), s: Structure::default() };

    // ---- consumeAtomSite
    const AS: &str = "atom_site";
    let label_asym_id = blk.col(AS, "label_asym_id")?;
    let auth_asym_id = blk.col(AS, "auth_asym_id")?;
    let group_pdb = blk.col(AS, "group_PDB")?;
    let auth_seq_id = blk.col(AS, "auth_seq_id")?;
    let label_comp_id = blk.col(AS, "label_comp_id")?;
    let id = blk.col(AS, "id")?;
    let label_atom_id = blk.col(AS, "label_atom_id")?;
    let cartn_x = blk.col(AS, "Cartn_x")?;
    let cartn_y = blk.col(AS, "Cartn_y")?;
    let cartn_z = blk.col(AS, "Cartn_z")?;
    let occupancy = blk.col(AS, "occupancy")?;
    let b_iso = blk.col(AS, "B_iso_or_equiv")?;
    let label_alt_id = blk.col(AS, "label_alt_id")?;
    let type_symbol = blk.col(AS, "type_symbol")?;
    let ins_code_col = blk.col(AS, "pdbx_PDB_ins_code")?;
    let label_seq_id = blk.col(AS, "label_seq_id")?;
    let model_num = blk.col(AS, "pdbx_PDB_model_num")?;

    let mut all_models: Vec<Vec<usize>> = Vec::new();
    let mut current_model: Vec<usize> = Vec::new();
    let mut current_chain: Option<usize> = None;
    let mut current_group: Option<usize> = None;
    let mut current_model_num: Option<String> = None;

    let n_rows = blk.row_count(AS);
    for i in 0..n_rows {
        let mut start_of_new_chain = false;
        let comp_id = label_comp_id.str(i)?;
        let mut one_letter = amino_1letter(&comp_id);
        let record = group_pdb.str(i)?;
        let mut het = false;
        if record != "ATOM" {
            if one_letter == Some('X') {
                one_letter = None;
            }
            het = true;
        }
        let ins_str = if ins_code_col.defined() { Some(ins_code_col.str(i)?) } else { None };
        let ins_code = match ins_str.as_deref() {
            Some(s) if !s.is_empty() && s != "?" => s.chars().next(),
            _ => None,
        };
        let seq_id = label_seq_id.int(i)?;
        // getStringData: '.' and '?' are empty
        let nmr_model = model_num.string_data(i)?;
        if current_model_num.is_none() {
            current_model_num = Some(nmr_model.clone());
        }
        if current_model_num.as_deref() != Some(nmr_model.as_str()) {
            current_model_num = Some(nmr_model.clone());
            if let Some(cc) = current_chain {
                if let Some(cg) = current_group {
                    b.add_group_to_chain(cc, cg);
                }
            }
            all_models.push(std::mem::take(&mut current_model));
            current_chain = None;
            current_group = None;
        }
        let asym_id = label_asym_id.str(i)?;
        let auth_id = if auth_asym_id.defined() { auth_asym_id.str(i)? } else { asym_id.clone() };

        if current_chain.is_none() {
            let c = b.new_chain(&auth_id, &asym_id);
            current_model.push(c);
            current_chain = Some(c);
            start_of_new_chain = true;
        }
        let cc = current_chain.unwrap();
        if asym_id != b.s.chains[cc].id {
            start_of_new_chain = true;
            if let Some(cg) = current_group {
                b.add_group_to_chain(cc, cg);
            }
            let found = current_model.iter().copied().find(|&c| b.s.chains[c].id == asym_id);
            let nc = match found {
                Some(c) => c,
                None => b.new_chain(&auth_id, &asym_id),
            };
            current_chain = Some(nc);
            if !current_model.contains(&nc) {
                current_model.push(nc);
            }
        }
        let auth_seq = if auth_seq_id.defined() { auth_seq_id.int(i)? } else { seq_id };
        let resnum = ResNum { chain: auth_id.clone(), seq: auth_seq, icode: ins_code };

        if current_group.is_none() {
            let g = b.new_group(&record, one_letter, &comp_id, resnum.clone(), het);
            current_group = Some(g);
        }
        let mut alt_group: Option<usize> = None;
        let alt_location = if label_alt_id.defined() { Some(label_alt_id.str(i)?) } else { None };
        if start_of_new_chain {
            let g = b.new_group(&record, one_letter, &comp_id, resnum.clone(), het);
            current_group = Some(g);
        } else {
            let cg = current_group.unwrap();
            if resnum != b.s.groups[cg].resnum {
                b.add_group_to_chain(current_chain.unwrap(), cg);
                let g = b.new_group(&record, one_letter, &comp_id, resnum.clone(), het);
                current_group = Some(g);
            } else if let Some(al) = alt_location.as_deref() {
                if !al.is_empty() && al != "." {
                    let ag = b.alt_loc_group(cg, &record, al.chars().next().unwrap(), one_letter, &comp_id);
                    if b.s.groups[ag].chain.is_none() {
                        b.s.groups[ag].chain = current_chain;
                    }
                    alt_group = Some(ag);
                }
            }
        }

        let serial = id.int(i)?;
        let name = label_atom_id.str(i)?;
        let x = cartn_x.float(i)?;
        let y = cartn_y.float(i)?;
        let z = cartn_z.float(i)?;
        let occ = if occupancy.defined() { occupancy.float(i)? } else { 1.0 };
        let tf = b_iso.float(i)? as f32;
        let alt_loc = match alt_location.as_deref() {
            None | Some("") | Some(".") => ' ',
            Some(s) => s.chars().next().unwrap(),
        };
        let ts = type_symbol.str(i)?;
        let element = element_value_of_ignore_case(&ts).unwrap_or(ELEMENT_R);
        let cg = current_group.unwrap();
        b.s.atoms.push(Atom { serial, name: name.clone(), alt_loc, x, y, z, occupancy: occ as f32, temp_factor: tf, element, group: cg });
        let a = b.s.atoms.len() - 1;
        if let Some(ag) = alt_group {
            b.add_atom(ag, a);
        } else {
            b.add_atom(cg, a);
        }
        if !b.group_has_atom(cg, &name) {
            let atom_group = b.s.atoms[a].group;
            if b.s.groups[cg].pdb_name == b.s.groups[atom_group].pdb_name {
                // StructureTools.hasNonDeuteratedEquiv
                let deut_equiv = b.s.atoms[a].element == crate::elements::ELEMENT_D
                    && name.starts_with('D')
                    && b.group_has_atom(cg, &name.replacen('D', "H", 1));
                if !deut_equiv {
                    b.add_atom(cg, a);
                    // the atom's group pointer is now the main group
                }
            }
        }
    }

    // ---- finish()
    if let Some(cc) = current_chain {
        if let Some(cg) = current_group {
            b.add_group_to_chain(cc, cg);
        }
        let cid = b.s.chains[cc].id.clone();
        if !current_model.iter().any(|&c| b.s.chains[c].id == cid) {
            current_model.push(cc);
        }
    }
    all_models.push(current_model);

    // asym id -> entity id (_struct_asym), entity id -> type (_entity)
    const SA: &str = "struct_asym";
    let sa_id = blk.col(SA, "id")?;
    let sa_entity = blk.col(SA, "entity_id")?;
    let ent_id = blk.col("entity", "id")?;
    let ent_type = blk.col("entity", "type")?;
    let n_ent = blk.row_count("entity");
    let mut asym2entity: HashMap<String, String> = HashMap::new();
    let mut entity_types: HashMap<i32, Option<EntityTypeX>> = HashMap::new();
    let n_sa = blk.row_count(SA);
    for r in 0..n_sa {
        asym2entity.insert(sa_id.str(r)?, sa_entity.str(r)?);
    }
    for r in 0..n_sa {
        let entity_id = sa_entity.str(r)?;
        // getEntityType: NoSuchElementException if the entity is missing in _entity
        let mut etype: Option<String> = None;
        for k in 0..n_ent {
            if ent_id.str(k)? == entity_id {
                etype = Some(ent_type.str(k)?);
                break;
            }
        }
        let etype = etype.ok_or_else(|| format!("NoSuchElementException: could not find entity with id {}", entity_id))?;
        let eid = parse_java_int(&entity_id).unwrap_or(0);
        entity_types.entry(eid).or_insert_with(|| entity_type_from_string(&etype));
    }

    let mut any_entity = !entity_types.is_empty();
    for model in &all_models {
        for &c in model {
            let cid = b.s.chains[c].id.clone();
            let entity_id = match asym2entity.get(&cid) {
                Some(e) => e.clone(),
                None => continue,
            };
            let eid = parse_java_int(&entity_id).ok_or_else(|| format!("NumberFormatException: For input string: \"{}\"", entity_id))?;
            match entity_types.get(&eid) {
                Some(t) => b.s.chains[c].entity = t.and_then(to_chain_entity),

                None => {
                    let t = if b.is_water_only(c) { EntityType::Water } else { EntityType::NonPolymer };
                    entity_types.insert(eid, Some(match t {
                        EntityType::Water => EntityTypeX::Water,
                        _ => EntityTypeX::NonPolymer,
                    }));
                    b.s.chains[c].entity = Some(t);
                    any_entity = true;
                }
            }
        }
    }

    if !any_entity {
        // no entity information: EntityFinder heuristics as in the PDB parser
        let mut poly_models = Vec::new();
        let mut non_poly_models = Vec::new();
        let mut water_models = Vec::new();
        for model in &all_models {
            let (mut p, mut np, mut w) = (Vec::new(), Vec::new(), Vec::new());
            for &c in model {
                if b.is_water_only(c) {
                    w.push(c);
                } else if b.is_pure_non_polymer(c)? {
                    np.push(c);
                } else {
                    p.push(c);
                }
            }
            poly_models.push(p);
            non_poly_models.push(np);
            water_models.push(w);
        }
        if let Some(first) = poly_models.first() {
            let first_ids: Vec<String> = first.iter().map(|&c| b.s.chains[c].id.clone()).collect();
            for (mi, model) in poly_models.iter().enumerate() {
                for &c in model {
                    if mi > 0 && !first_ids.contains(&b.s.chains[c].id) {
                        return Err(format!("NullPointerException: EntityFinder: chain {} not found in first model", b.s.chains[c].id));
                    }
                    b.s.chains[c].entity = Some(EntityType::Polymer);
                }
            }
        }
        if non_poly_models.first().map(|m| !m.is_empty()).unwrap_or(false) {
            for m in &non_poly_models {
                for &c in m {
                    b.s.chains[c].entity = Some(EntityType::NonPolymer);
                }
            }
        }
        if water_models.first().map(|m| !m.is_empty()).unwrap_or(false) {
            for m in &water_models {
                for &c in m {
                    b.s.chains[c].entity = Some(EntityType::Water);
                }
            }
        }
    }

    // structure.addModel -> Model.setChains: polymers (and chains without entity), non-polymers, water
    for model in all_models {
        let mut poly = Vec::new();
        let mut nonpoly = Vec::new();
        let mut water = Vec::new();
        for c in model {
            match b.s.chains[c].entity {
                None | Some(EntityType::Polymer) => poly.push(c),
                Some(EntityType::NonPolymer) => nonpoly.push(c),
                Some(EntityType::Water) => water.push(c),
            }
        }
        let mut m = poly;
        m.extend(nonpoly);
        m.extend(water);
        b.s.models.push(m);
    }
    Ok(b.s)
}

/// All BioJava entity types (branched/macrolide are stored as non-polymer chains)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum EntityTypeX {
    Polymer,
    Branched,
    NonPolymer,
    Water,
    Macrolide,
}

fn entity_type_from_string(s: &str) -> Option<EntityTypeX> {
    for (n, t) in [
        ("polymer", EntityTypeX::Polymer),
        ("branched", EntityTypeX::Branched),
        ("non-polymer", EntityTypeX::NonPolymer),
        ("water", EntityTypeX::Water),
        ("macrolide", EntityTypeX::Macrolide),
    ] {
        if s.eq_ignore_ascii_case(n) {
            return Some(t);
        }
    }
    None
}

fn to_chain_entity(t: EntityTypeX) -> Option<EntityType> {
    Some(match t {
        EntityTypeX::Polymer => EntityType::Polymer,
        EntityTypeX::Water => EntityType::Water,
        EntityTypeX::NonPolymer | EntityTypeX::Branched | EntityTypeX::Macrolide => EntityType::NonPolymer,
    })
}
