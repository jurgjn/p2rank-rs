//! Minimal reader for the Java Object Serialization stream protocol, sufficient to load a
//! FasterForest `LegacyFlatBinaryForest` from P2Rank's `model.zst` (zstd-compressed
//! `ObjectOutputStream` output) without a JVM.

use crate::forest::Forest;
use std::collections::HashMap;
use std::path::Path;

const TC_NULL: u8 = 0x70;
const TC_REFERENCE: u8 = 0x71;
const TC_CLASSDESC: u8 = 0x72;
const TC_OBJECT: u8 = 0x73;
const TC_STRING: u8 = 0x74;
const TC_ARRAY: u8 = 0x75;
const TC_CLASS: u8 = 0x76;
const TC_BLOCKDATA: u8 = 0x77;
const TC_ENDBLOCKDATA: u8 = 0x78;
const TC_RESET: u8 = 0x79;
const TC_BLOCKDATALONG: u8 = 0x7A;
const TC_LONGSTRING: u8 = 0x7C;
const TC_PROXYCLASSDESC: u8 = 0x7D;
const TC_ENUM: u8 = 0x7E;
const BASE_WIRE_HANDLE: u32 = 0x7E0000;

const SC_WRITE_METHOD: u8 = 0x01;
const SC_SERIALIZABLE: u8 = 0x02;
const SC_EXTERNALIZABLE: u8 = 0x04;
const SC_BLOCK_DATA: u8 = 0x08;

#[allow(dead_code)] // mirrors the Java data model
#[derive(Debug, Clone)]
pub enum Value {
    Null,
    Prim(f64, i64),
    Str(String),
    Obj(usize),
    Arr(usize),
    Class,
}

#[derive(Debug)]
struct ClassDesc {
    name: String,
    flags: u8,
    fields: Vec<(u8, String)>,
    super_desc: Option<usize>,
}

#[allow(dead_code)] // mirrors the Java data model
#[derive(Debug)]
pub enum Array {
    I32(Vec<i32>),
    F64(Vec<f64>),
    Other(Vec<Value>),
    Bytes(usize),
}

pub struct Object {
    pub class: String,
    pub fields: HashMap<String, Value>,
}

enum Handle {
    Desc(usize),
    Obj(usize),
    Arr(usize),
    Str(String),
    Class,
}

struct Reader<'a> {
    b: &'a [u8],
    i: usize,
    descs: Vec<ClassDesc>,
    objs: Vec<Object>,
    arrays: Vec<Array>,
    handles: Vec<Handle>,
}

type R<T> = Result<T, String>;

impl<'a> Reader<'a> {
    fn u8(&mut self) -> R<u8> {
        let v = *self.b.get(self.i).ok_or("unexpected end of stream")?;
        self.i += 1;
        Ok(v)
    }
    fn take(&mut self, n: usize) -> R<&'a [u8]> {
        if self.i + n > self.b.len() {
            return Err("unexpected end of stream".into());
        }
        let s = &self.b[self.i..self.i + n];
        self.i += n;
        Ok(s)
    }
    fn u16(&mut self) -> R<u16> {
        Ok(u16::from_be_bytes(self.take(2)?.try_into().unwrap()))
    }
    fn i32(&mut self) -> R<i32> {
        Ok(i32::from_be_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn i64(&mut self) -> R<i64> {
        Ok(i64::from_be_bytes(self.take(8)?.try_into().unwrap()))
    }
    fn utf(&mut self) -> R<String> {
        let n = self.u16()? as usize;
        Ok(String::from_utf8_lossy(self.take(n)?).to_string())
    }
    fn long_utf(&mut self) -> R<String> {
        let n = self.i64()? as usize;
        Ok(String::from_utf8_lossy(self.take(n)?).to_string())
    }

    fn new_handle(&mut self, h: Handle) -> usize {
        self.handles.push(h);
        self.handles.len() - 1
    }

    fn handle(&mut self) -> R<usize> {
        let h = self.i32()? as u32;
        h.checked_sub(BASE_WIRE_HANDLE).map(|x| x as usize).ok_or_else(|| "invalid handle".to_string())
    }

    fn string_value(&mut self) -> R<String> {
        match self.u8()? {
            TC_STRING => {
                let s = self.utf()?;
                self.new_handle(Handle::Str(s.clone()));
                Ok(s)
            }
            TC_LONGSTRING => {
                let s = self.long_utf()?;
                self.new_handle(Handle::Str(s.clone()));
                Ok(s)
            }
            TC_REFERENCE => {
                let h = self.handle()?;
                match self.handles.get(h) {
                    Some(Handle::Str(s)) => Ok(s.clone()),
                    _ => Err("reference is not a string".into()),
                }
            }
            t => Err(format!("expected string, got tag {:#x}", t)),
        }
    }

    fn class_desc(&mut self) -> R<Option<usize>> {
        match self.u8()? {
            TC_NULL => Ok(None),
            TC_REFERENCE => {
                let h = self.handle()?;
                match self.handles.get(h) {
                    Some(Handle::Desc(d)) => Ok(Some(*d)),
                    _ => Err("reference is not a class desc".into()),
                }
            }
            TC_CLASSDESC => {
                let name = self.utf()?;
                let _uid = self.i64()?;
                let idx = self.descs.len();
                self.descs.push(ClassDesc { name, flags: 0, fields: Vec::new(), super_desc: None });
                self.new_handle(Handle::Desc(idx));
                let flags = self.u8()?;
                let n = self.u16()?;
                let mut fields = Vec::new();
                for _ in 0..n {
                    let t = self.u8()?;
                    let fname = self.utf()?;
                    if t == b'L' || t == b'[' {
                        let _cn = self.string_value()?;
                    }
                    fields.push((t, fname));
                }
                self.skip_annotation()?;
                let sup = self.class_desc()?;
                let d = &mut self.descs[idx];
                d.flags = flags;
                d.fields = fields;
                d.super_desc = sup;
                Ok(Some(idx))
            }
            TC_PROXYCLASSDESC => Err("proxy classes are not supported".into()),
            t => Err(format!("expected class desc, got tag {:#x}", t)),
        }
    }

    fn skip_annotation(&mut self) -> R<()> {
        loop {
            match self.b.get(self.i).copied() {
                Some(TC_ENDBLOCKDATA) => {
                    self.i += 1;
                    return Ok(());
                }
                Some(TC_BLOCKDATA) => {
                    self.i += 1;
                    let n = self.u8()? as usize;
                    self.take(n)?;
                }
                Some(TC_BLOCKDATALONG) => {
                    self.i += 1;
                    let n = self.i32()? as usize;
                    self.take(n)?;
                }
                Some(_) => {
                    self.content()?;
                }
                None => return Err("unexpected end of stream".into()),
            }
        }
    }

    fn prim(&mut self, t: u8) -> R<Value> {
        Ok(match t {
            b'B' => Value::Prim(0.0, self.u8()? as i8 as i64),
            b'C' => Value::Prim(0.0, self.u16()? as i64),
            b'S' => Value::Prim(0.0, self.u16()? as i16 as i64),
            b'Z' => Value::Prim(0.0, self.u8()? as i64),
            b'I' => Value::Prim(0.0, self.i32()? as i64),
            b'J' => Value::Prim(0.0, self.i64()?),
            b'F' => Value::Prim(f32::from_bits(self.i32()? as u32) as f64, 0),
            b'D' => Value::Prim(f64::from_bits(self.i64()? as u64), 0),
            _ => return Err(format!("bad primitive type {}", t as char)),
        })
    }

    fn content(&mut self) -> R<Value> {
        let tag = self.u8()?;
        match tag {
            TC_NULL => Ok(Value::Null),
            TC_REFERENCE => {
                let h = self.handle()?;
                match self.handles.get(h) {
                    Some(Handle::Obj(o)) => Ok(Value::Obj(*o)),
                    Some(Handle::Arr(a)) => Ok(Value::Arr(*a)),
                    Some(Handle::Str(s)) => Ok(Value::Str(s.clone())),
                    Some(Handle::Class) => Ok(Value::Class),
                    Some(Handle::Desc(_)) => Ok(Value::Class),
                    None => Err("dangling handle".into()),
                }
            }
            TC_STRING | TC_LONGSTRING => {
                self.i -= 1;
                Ok(Value::Str(self.string_value()?))
            }
            TC_CLASS => {
                self.class_desc()?;
                self.new_handle(Handle::Class);
                Ok(Value::Class)
            }
            TC_CLASSDESC => {
                self.i -= 1;
                self.class_desc()?;
                Ok(Value::Class)
            }
            TC_ARRAY => {
                let d = self.class_desc()?.ok_or("array without class")?;
                let name = self.descs[d].name.clone();
                let ai = self.arrays.len();
                self.arrays.push(Array::Bytes(0));
                self.new_handle(Handle::Arr(ai));
                let n = self.i32()?;
                if n < 0 {
                    return Err("negative array size".into());
                }
                let n = n as usize;
                let arr = match name.as_str() {
                    "[I" => {
                        let raw = self.take(n * 4)?;
                        Array::I32(raw.chunks_exact(4).map(|c| i32::from_be_bytes(c.try_into().unwrap())).collect())
                    }
                    "[D" => {
                        let raw = self.take(n * 8)?;
                        Array::F64(raw.chunks_exact(8).map(|c| f64::from_bits(u64::from_be_bytes(c.try_into().unwrap()))).collect())
                    }
                    "[B" | "[Z" => {
                        self.take(n)?;
                        Array::Bytes(n)
                    }
                    "[S" | "[C" => {
                        self.take(n * 2)?;
                        Array::Bytes(n)
                    }
                    "[F" => {
                        self.take(n * 4)?;
                        Array::Bytes(n)
                    }
                    "[J" => {
                        self.take(n * 8)?;
                        Array::Bytes(n)
                    }
                    _ => {
                        let mut v = Vec::with_capacity(n);
                        for _ in 0..n {
                            v.push(self.content()?);
                        }
                        Array::Other(v)
                    }
                };
                self.arrays[ai] = arr;
                Ok(Value::Arr(ai))
            }
            TC_OBJECT => {
                let d = self.class_desc()?.ok_or("object without class")?;
                let oi = self.objs.len();
                self.objs.push(Object { class: self.descs[d].name.clone(), fields: HashMap::new() });
                self.new_handle(Handle::Obj(oi));
                // class hierarchy from the top-most superclass down
                let mut chain = vec![d];
                while let Some(s) = self.descs[*chain.last().unwrap()].super_desc {
                    chain.push(s);
                }
                chain.reverse();
                for cd in chain {
                    let flags = self.descs[cd].flags;
                    if flags & SC_EXTERNALIZABLE != 0 {
                        if flags & SC_BLOCK_DATA != 0 {
                            self.skip_annotation()?;
                            continue;
                        }
                        return Err("externalizable objects without block data are not supported".into());
                    }
                    if flags & SC_SERIALIZABLE == 0 {
                        continue;
                    }
                    let fields = self.descs[cd].fields.clone();
                    for (t, name) in fields {
                        let v = if t == b'L' || t == b'[' { self.content()? } else { self.prim(t)? };
                        self.objs[oi].fields.insert(name, v);
                    }
                    if flags & SC_WRITE_METHOD != 0 {
                        self.skip_annotation()?;
                    }
                }
                Ok(Value::Obj(oi))
            }
            TC_ENUM => {
                self.class_desc()?;
                self.new_handle(Handle::Class);
                self.string_value()?;
                Ok(Value::Class)
            }
            TC_RESET => {
                self.handles.clear();
                self.content()
            }
            TC_BLOCKDATA => {
                let n = self.u8()? as usize;
                self.take(n)?;
                self.content()
            }
            TC_BLOCKDATALONG => {
                let n = self.i32()? as usize;
                self.take(n)?;
                self.content()
            }
            t => Err(format!("unsupported serialization tag {:#x} at offset {}", t, self.i - 1)),
        }
    }
}

fn decompress_zstd(raw: &[u8]) -> R<Vec<u8>> {
    use std::io::Read;
    let mut dec = ruzstd::decoding::StreamingDecoder::new(raw).map_err(|e| format!("zstd: {:?}", e))?;
    let mut out = Vec::with_capacity(raw.len() * 4);
    dec.read_to_end(&mut out).map_err(|e| format!("zstd: {}", e))?;
    Ok(out)
}

/// Loads a `LegacyFlatBinaryForest` serialized by P2Rank into `model.zst`.
pub fn load_legacy_flat_forest_zst(path: &Path) -> R<Forest> {
    let raw = std::fs::read(path).map_err(|e| format!("cannot read model {}: {}", path.display(), e))?;
    let data = decompress_zstd(&raw)?;
    if data.len() < 4 || data[0..2] != [0xAC, 0xED] {
        return Err("model is not a Java serialization stream".into());
    }
    let mut r = Reader { b: &data, i: 4, descs: Vec::new(), objs: Vec::new(), arrays: Vec::new(), handles: Vec::new() };
    let root = r.content()?;
    let oi = match root {
        Value::Obj(o) => o,
        _ => return Err("model root is not an object".into()),
    };
    let obj = &r.objs[oi];
    if obj.class != "cz.siret.prank.fforest.api.LegacyFlatBinaryForest" {
        return Err(format!("unsupported model class {} (only LegacyFlatBinaryForest is supported)", obj.class));
    }
    let int_field = |name: &str| -> R<i64> {
        match obj.fields.get(name) {
            Some(Value::Prim(_, i)) => Ok(*i),
            _ => Err(format!("missing field {}", name)),
        }
    };
    let arr = |name: &str| -> R<&Array> {
        match obj.fields.get(name) {
            Some(Value::Arr(a)) => Ok(&r.arrays[*a]),
            _ => Err(format!("missing array field {}", name)),
        }
    };
    let ints = |name: &str| -> R<Vec<i32>> {
        match arr(name)? {
            Array::I32(v) => Ok(v.clone()),
            _ => Err(format!("field {} is not int[]", name)),
        }
    };
    let num_trees = int_field("numTrees")? as usize;
    let num_attributes = int_field("numAttributes")? as usize;
    let left = ints("childLeft")?;
    let right = ints("childRight")?;
    let attr = ints("attributeIndex")?;
    let split = match arr("splitPoint")? {
        Array::F64(v) => v.clone(),
        _ => return Err("splitPoint is not double[]".into()),
    };
    let leaves = match arr("classProbs")? {
        Array::Other(rows) => {
            let mut out = Vec::with_capacity(rows.len());
            for row in rows {
                match row {
                    Value::Arr(a) => match &r.arrays[*a] {
                        Array::F64(p) if p.len() == 2 => out.push([p[0], p[1]]),
                        _ => return Err("invalid classProbs row".into()),
                    },
                    Value::Null => out.push([f64::NAN, f64::NAN]),
                    _ => return Err("invalid classProbs row".into()),
                }
            }
            out
        }
        _ => return Err("classProbs is not double[][]".into()),
    };
    if left.len() != right.len() || left.len() != attr.len() || left.len() != split.len() {
        return Err("inconsistent forest arrays".into());
    }
    Forest::from_parts(num_trees, num_attributes, left, right, attr, split, leaves)
}
