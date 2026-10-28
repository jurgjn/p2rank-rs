//! BinaryCIF decoding (port of ciftools-java MessagePackCodec, BinaryCifCodec and encodings).

use std::collections::BTreeMap;

#[allow(dead_code)] // mirrors the Java data model
#[derive(Debug, Clone)]
pub enum Mp {
    Nil,
    Bool(bool),
    Int(i64),
    Float(f64),
    Str(String),
    Bin(Vec<u8>),
    Arr(Vec<Mp>),
    Map(Vec<(String, Mp)>),
}

impl Mp {
    pub fn get(&self, k: &str) -> Option<&Mp> {
        match self {
            Mp::Map(v) => v.iter().find(|e| e.0 == k).map(|e| &e.1),
            _ => None,
        }
    }
    fn int(&self, k: &str) -> Result<i64, String> {
        match self.get(k) {
            Some(Mp::Int(i)) => Ok(*i),
            _ => Err(format!("ClassCastException: missing integer '{}'", k)),
        }
    }
    pub fn str_(&self, k: &str) -> Result<&str, String> {
        match self.get(k) {
            Some(Mp::Str(s)) => Ok(s),
            _ => Err(format!("ClassCastException: missing string '{}'", k)),
        }
    }
    pub fn arr(&self, k: &str) -> Result<&Vec<Mp>, String> {
        match self.get(k) {
            Some(Mp::Arr(a)) => Ok(a),
            _ => Err(format!("ClassCastException: missing array '{}'", k)),
        }
    }
    fn bin(&self, k: &str) -> Result<&Vec<u8>, String> {
        match self.get(k) {
            Some(Mp::Bin(b)) => Ok(b),
            _ => Err(format!("ClassCastException: missing binary '{}'", k)),
        }
    }
}

struct Rd<'a> {
    b: &'a [u8],
    i: usize,
}

impl<'a> Rd<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8], String> {
        if self.i + n > self.b.len() {
            return Err("EOFException".into());
        }
        let s = &self.b[self.i..self.i + n];
        self.i += n;
        Ok(s)
    }
    fn u8(&mut self) -> Result<u8, String> {
        Ok(self.take(1)?[0])
    }
    fn u16(&mut self) -> Result<u16, String> {
        Ok(u16::from_be_bytes(self.take(2)?.try_into().unwrap()))
    }
    fn u32(&mut self) -> Result<u32, String> {
        Ok(u32::from_be_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn s(&mut self, n: usize) -> Result<String, String> {
        // new String(bytes, US_ASCII): non-ASCII bytes become U+FFFD
        Ok(self.take(n)?.iter().map(|&c| if c < 128 { c as char } else { '\u{FFFD}' }).collect())
    }
    fn value(&mut self) -> Result<Mp, String> {
        let t = self.u8()?;
        if t & 0x80 == 0 {
            return Ok(Mp::Int(t as i64));
        }
        if t & 0xf0 == 0x80 {
            return self.map((t & 15) as usize);
        }
        if t & 0xf0 == 0x90 {
            return self.array((t & 15) as usize);
        }
        if t & 0xe0 == 0xa0 {
            return Ok(Mp::Str(self.s((t & 31) as usize)?));
        }
        if t & 0xe0 == 0xe0 {
            return Ok(Mp::Int(t as i8 as i64));
        }
        Ok(match t {
            0xc0 => Mp::Nil,
            0xc2 => Mp::Bool(false),
            0xc3 => Mp::Bool(true),
            0xc4 => {
                let n = self.u8()? as usize;
                Mp::Bin(self.take(n)?.to_vec())
            }
            0xc5 => {
                let n = self.u16()? as usize;
                Mp::Bin(self.take(n)?.to_vec())
            }
            0xc6 => {
                let n = self.u32()? as usize;
                Mp::Bin(self.take(n)?.to_vec())
            }
            0xca => Mp::Float(f32::from_bits(self.u32()?) as f64),
            0xcb => Mp::Float(f64::from_bits(u64::from_be_bytes(self.take(8)?.try_into().unwrap()))),
            0xcc => Mp::Int(self.u8()? as i64),
            0xcd => Mp::Int(self.u16()? as i64),
            // readUnsignedInt returns a (possibly negative) Java int
            0xce => Mp::Int(self.u32()? as i32 as i64),
            0xd0 => Mp::Int(self.u8()? as i8 as i64),
            0xd1 => Mp::Int(self.u16()? as i16 as i64),
            0xd2 => Mp::Int(self.u32()? as i32 as i64),
            0xd9 => {
                let n = self.u8()? as usize;
                Mp::Str(self.s(n)?)
            }
            0xda => {
                let n = self.u16()? as usize;
                Mp::Str(self.s(n)?)
            }
            0xdb => {
                let n = self.u32()? as usize;
                Mp::Str(self.s(n)?)
            }
            0xdc => {
                let n = self.u16()? as usize;
                return self.array(n);
            }
            0xdd => {
                let n = self.u32()? as usize;
                return self.array(n);
            }
            0xde => {
                let n = self.u16()? as usize;
                return self.map(n);
            }
            0xdf => {
                let n = self.u32()? as usize;
                return self.map(n);
            }
            _ => return Err(format!("IllegalArgumentException: Unknown MessagePack type 0x{:x}", t)),
        })
    }
    fn map(&mut self, n: usize) -> Result<Mp, String> {
        let mut v: Vec<(String, Mp)> = Vec::with_capacity(n.min(1 << 16));
        for _ in 0..n {
            let k = match self.value()? {
                Mp::Str(s) => s,
                _ => return Err("ClassCastException: map key is not a string".into()),
            };
            let val = self.value()?;
            // LinkedHashMap.put: replace existing key in place
            if let Some(e) = v.iter_mut().find(|e| e.0 == k) {
                e.1 = val;
            } else {
                v.push((k, val));
            }
        }
        Ok(Mp::Map(v))
    }
    fn array(&mut self, n: usize) -> Result<Mp, String> {
        let mut v = Vec::with_capacity(n.min(1 << 16));
        for _ in 0..n {
            v.push(self.value()?);
        }
        Ok(Mp::Arr(v))
    }
}

pub fn decode_msgpack(b: &[u8]) -> Result<Mp, String> {
    Rd { b, i: 0 }.value()
}

// ------------------------------------------------------------------ encodings

#[allow(dead_code)] // mirrors the Java data model
#[derive(Debug, Clone)]
pub enum Data {
    Bytes(Vec<u8>),
    /// int arrays with their type code (1..6 = int8..uint32)
    Ints(Vec<i32>, i64),
    /// float arrays with type code 32 / 33
    Floats(Vec<f64>, i64),
    Strs(Vec<String>),
}

fn bytes_to(b: &[u8], ty: i64) -> Result<Data, String> {
    Ok(match ty {
        1 => Data::Ints(b.iter().map(|&x| x as i8 as i32).collect(), 1),
        2 => Data::Ints(b.chunks_exact(2).map(|c| i16::from_le_bytes([c[0], c[1]]) as i32).collect(), 2),
        3 => Data::Ints(b.chunks_exact(4).map(|c| i32::from_le_bytes(c.try_into().unwrap())).collect(), 3),
        4 => Data::Ints(b.iter().map(|&x| x as i32).collect(), 4),
        5 => Data::Ints(b.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]]) as i32).collect(), 5),
        6 => Data::Ints(b.chunks_exact(4).map(|c| u32::from_le_bytes(c.try_into().unwrap()) as i32).collect(), 6),
        32 => Data::Floats(b.chunks_exact(4).map(|c| f32::from_le_bytes(c.try_into().unwrap()) as f64).collect(), 32),
        33 => Data::Floats(b.chunks_exact(8).map(|c| f64::from_le_bytes(c.try_into().unwrap())).collect(), 33),
        _ => return Err(format!("UnsupportedOperationException: cannot handle byte array encoding type {}", ty)),
    })
}

fn apply(enc: &Mp, data: Data) -> Result<Data, String> {
    let kind = enc.str_("kind")?;
    match kind {
        "ByteArray" => match data {
            Data::Bytes(b) => bytes_to(&b, enc.int("type")?),
            _ => Err("ClassCastException: ByteArray decoding of non-bytes".into()),
        },
        "FixedPoint" => {
            let factor = enc.int("factor")? as i32;
            let src = enc.int("srcType")?;
            match data {
                Data::Ints(v, _) => {
                    let f = 1.0f64 / factor as f64;
                    Ok(Data::Floats(v.iter().map(|&x| x as f64 * f).collect(), if src == 32 { 32 } else { 33 }))
                }
                _ => Err("ClassCastException: FixedPoint of non-int data".into()),
            }
        }
        "IntervalQuantization" => {
            let min = enc.int("min")? as i32 as f64;
            let max = enc.int("max")? as i32 as f64;
            let steps = enc.int("numSteps")? as i32;
            let src = enc.int("srcType")?;
            match data {
                Data::Ints(v, _) => {
                    let delta = (max - min) / (steps as f64 - 1.0);
                    Ok(Data::Floats(v.iter().map(|&x| min + delta * x as f64).collect(), if src == 32 { 32 } else { 33 }))
                }
                _ => Err("ClassCastException: IntervalQuantization of non-int data".into()),
            }
        }
        "RunLength" => {
            let src_type = enc.int("srcType")?;
            let src_size = enc.int("srcSize")? as usize;
            match data {
                Data::Ints(v, t) => {
                    if v.is_empty() {
                        return Ok(Data::Ints(v, t));
                    }
                    let mut out = vec![0i32; src_size];
                    let mut off = 0usize;
                    let mut i = 0;
                    while i < v.len() {
                        let value = v[i];
                        let len = *v.get(i + 1).ok_or("ArrayIndexOutOfBoundsException")?;
                        for _ in 0..len.max(0) {
                            *out.get_mut(off).ok_or("ArrayIndexOutOfBoundsException")? = value;
                            off += 1;
                        }
                        i += 2;
                    }
                    Ok(Data::Ints(out, src_type))
                }
                _ => Err("ClassCastException: RunLength of non-int data".into()),
            }
        }
        "Delta" => {
            let origin = enc.int("origin")? as i32;
            match data {
                Data::Ints(v, t) => {
                    let mut out = vec![0i32; v.len()];
                    if !v.is_empty() {
                        out[0] = v[0].wrapping_add(origin);
                        for i in 1..v.len() {
                            out[i] = v[i].wrapping_add(out[i - 1]);
                        }
                    }
                    Ok(Data::Ints(out, t))
                }
                _ => Err("ClassCastException: Delta of non-int data".into()),
            }
        }
        "IntegerPacking" => {
            let byte_count = enc.int("byteCount")?;
            let unsigned = matches!(enc.get("isUnsigned"), Some(Mp::Bool(true)));
            let src_size = enc.int("srcSize")? as usize;
            match data {
                Data::Ints(input, _) => {
                    if input.len() == src_size && byte_count == 4 {
                        return Ok(Data::Ints(input, 3));
                    }
                    let (lower, upper) = if unsigned {
                        (0, if byte_count == 1 { 255 } else { 65535 })
                    } else {
                        let u = if byte_count == 1 { 127 } else { 32767 };
                        (-u - 1, u)
                    };
                    let mut out = vec![0i32; src_size];
                    let mut i = 0;
                    let mut j = 0;
                    while i < input.len() {
                        let mut value = 0i32;
                        let mut t = input[i];
                        while if unsigned { t == upper } else { t == upper || t == lower } {
                            value = value.wrapping_add(t);
                            i += 1;
                            t = *input.get(i).ok_or("ArrayIndexOutOfBoundsException")?;
                        }
                        value = value.wrapping_add(t);
                        *out.get_mut(j).ok_or("ArrayIndexOutOfBoundsException")? = value;
                        i += 1;
                        j += 1;
                    }
                    Ok(Data::Ints(out, 3))
                }
                _ => Err("ClassCastException: IntegerPacking of non-int data".into()),
            }
        }
        "StringArray" => {
            let string_data = enc.str_("stringData")?;
            let offsets_raw = enc.bin("offsets")?;
            let offsets = match decode_chain(Data::Bytes(offsets_raw.clone()), enc.arr("offsetEncoding")?)? {
                Data::Ints(v, _) => v,
                _ => return Err("ClassCastException: string offsets".into()),
            };
            let indices = match data {
                Data::Bytes(b) => match decode_chain(Data::Bytes(b), enc.arr("dataEncoding")?)? {
                    Data::Ints(v, _) => v,
                    _ => return Err("ClassCastException: string indices".into()),
                },
                _ => return Err("ClassCastException: StringArray of non-bytes".into()),
            };
            let chars: Vec<char> = string_data.chars().collect();
            let mut strings: Vec<String> = Vec::with_capacity(offsets.len());
            strings.push(String::new());
            for i in 1..offsets.len() {
                let (a, b) = (offsets[i - 1] as usize, offsets[i] as usize);
                if a > b || b > chars.len() {
                    return Err("StringIndexOutOfBoundsException".into());
                }
                strings.push(chars[a..b].iter().collect());
            }
            let mut out = Vec::with_capacity(indices.len());
            for idx in indices {
                out.push(strings.get((idx + 1) as usize).cloned().ok_or("ArrayIndexOutOfBoundsException")?);
            }
            Ok(Data::Strs(out))
        }
        _ => Err(format!("IllegalArgumentException: Unsupported Encoding kind: {}", kind)),
    }
}

/// Encodings are applied in order of the deque in StringArrayEncoding (Deque.removeLast in decode)
fn decode_chain(data: Data, encodings: &[Mp]) -> Result<Data, String> {
    let mut cur = data;
    for e in encodings.iter().rev() {
        cur = apply(e, cur)?;
    }
    Ok(cur)
}

/// BinaryCifCodec.decode(Map): encodings applied from last to first
pub fn decode_column_data(encoded: &Mp) -> Result<Data, String> {
    let data = encoded.bin("data")?.clone();
    decode_chain(Data::Bytes(data), encoded.arr("encoding")?)
}

/// first data block: category name (without '_') -> (rowCount, encoded column maps by name)
pub struct BCategory {
    pub row_count: usize,
    pub columns: Vec<(String, Mp)>,
}

pub fn read_first_block(bytes: &[u8]) -> Result<Vec<(String, BCategory)>, String> {
    let root = decode_msgpack(bytes).map_err(|e| format!("ParsingException: Parsing failed. ({})", e))?;
    let version = root.str_("version")?;
    if !version.starts_with("0.3") {
        return Err(format!("ParsingException: Unsupported format version. Current {}, required 0.3.", version));
    }
    let blocks = root.arr("dataBlocks")?;
    let block = blocks.first().ok_or("IndexOutOfBoundsException: Index 0 out of bounds for length 0")?;
    let mut out: Vec<(String, BCategory)> = Vec::new();
    let cats = match block.arr("categories") {
        Ok(c) => c,
        Err(_) => return Ok(out),
    };
    let mut seen: BTreeMap<String, usize> = BTreeMap::new();
    for cat in cats {
        let name = cat.str_("name")?;
        let name = name.get(1..).unwrap_or("").to_string();
        let row_count = cat.int("rowCount")? as usize;
        let mut columns = Vec::new();
        for c in cat.arr("columns")? {
            columns.push((c.str_("name")?.to_string(), c.clone()));
        }
        let key = name.to_lowercase();
        let entry = BCategory { row_count, columns };
        if let Some(&i) = seen.get(&key) {
            out[i] = (name, entry);
        } else {
            seen.insert(key, out.len());
            out.push((name, entry));
        }
    }
    Ok(out)
}
