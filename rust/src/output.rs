//! Output files: predictions/residues CSV, PyMOL and ChimeraX scripts, SAS points PDB (gzip).

use crate::jfmt::{decimal_format, format_f, format_wd, format_wf};
use crate::predict::{LabeledPoint, Pocket, ResidueRow};
use crate::protein::Protein;
use std::collections::BTreeSet;
use std::io::Write;

// ------------------------------------------------------------------ colors

/// Groovy `int / int` (BigDecimal division) converted to double, for 0 <= i < n
fn groovy_int_div_to_f64(i: i64, n: i64) -> f64 {
    decimal_to_f64(&groovy_int_div(i, n))
}

/// Groovy BigDecimalMath.divideImpl for non-negative ints: returns decimal string
fn groovy_int_div(i: i64, n: i64) -> String {
    if i == 0 {
        return "0".to_string();
    }
    let int_part = i / n;
    let mut rem = i % n;
    // terminating?
    let mut m = n / gcd(i, n);
    while m % 2 == 0 {
        m /= 2;
    }
    while m % 5 == 0 {
        m /= 5;
    }
    if m == 1 {
        let mut frac = String::new();
        while rem != 0 {
            rem *= 10;
            frac.push((b'0' + (rem / n) as u8) as char);
            rem %= n;
        }
        return if frac.is_empty() { int_part.to_string() } else { format!("{}.{}", int_part, frac) };
    }
    // non-terminating: round to (max(precision(i), precision(n)) + 10) significant digits HALF_UP,
    // then to scale 10 HALF_UP if the scale exceeds 10
    let prec = (i.to_string().len().max(n.to_string().len()) + 10) as i64;
    // generate digits of the full value (integer part + fraction) until prec significant digits + 1
    let mut digits: Vec<u8> = int_part.to_string().bytes().map(|b| b - b'0').collect();
    let int_len = digits.len() as i64;
    let mut sig_start: Option<i64> = if int_part != 0 { Some(0) } else { None };
    let mut pos = int_len; // number of digits generated so far
    loop {
        rem *= 10;
        let d = (rem / n) as u8;
        rem %= n;
        digits.push(d);
        pos += 1;
        if sig_start.is_none() && d != 0 {
            sig_start = Some(pos - 1);
        }
        if let Some(s) = sig_start {
            if pos - s >= prec + 1 {
                break;
            }
        }
    }
    let s = sig_start.unwrap();
    let keep = (s + prec) as usize; // digits to keep (int digits + fraction digits)
    let round_up = digits[keep] >= 5;
    digits.truncate(keep);
    let mut int_len = int_len as usize;
    if round_up && inc_digits(&mut digits) {
        digits.insert(0, 1);
        int_len += 1;
    }
    let mut scale = digits.len() - int_len;
    if scale > 10 {
        let round_up = digits[int_len + 10] >= 5;
        digits.truncate(int_len + 10);
        if round_up && inc_digits(&mut digits) {
            digits.insert(0, 1);
            int_len += 1;
        }
        scale = 10;
    }
    let ip: String = digits[..int_len].iter().map(|d| (b'0' + d) as char).collect();
    let fp: String = digits[int_len..int_len + scale].iter().map(|d| (b'0' + d) as char).collect();
    let ip = if ip.is_empty() { "0".to_string() } else { ip };
    format!("{}.{}", ip, fp)
}

fn inc_digits(d: &mut [u8]) -> bool {
    for i in (0..d.len()).rev() {
        if d[i] == 9 {
            d[i] = 0;
        } else {
            d[i] += 1;
            return false;
        }
    }
    true
}

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

fn decimal_to_f64(s: &str) -> f64 {
    s.parse::<f64>().unwrap()
}

/// java.awt.Color.HSBtoRGB
fn hsb_to_rgb(hue: f32, saturation: f32, brightness: f32) -> (i32, i32, i32) {
    let c = |v: f32| -> i32 { (v * 255.0f32 + 0.5f32) as i32 };
    if saturation == 0.0 {
        let v = c(brightness);
        return (v, v, v);
    }
    let h = (hue - (hue as f64).floor() as f32) * 6.0f32;
    let f = h - (h as f64).floor() as f32;
    let p = brightness * (1.0f32 - saturation);
    let q = brightness * (1.0f32 - saturation * f);
    let t = brightness * (1.0f32 - (saturation * (1.0f32 - f)));
    match h as i32 {
        0 => (c(brightness), c(t), c(p)),
        1 => (c(q), c(brightness), c(p)),
        2 => (c(p), c(brightness), c(t)),
        3 => (c(p), c(q), c(brightness)),
        4 => (c(t), c(p), c(brightness)),
        5 => (c(brightness), c(p), c(q)),
        _ => (0, 0, 0),
    }
}

/// ColorUtils.createSpectrum(n, 0.6, 0.6, 1.2)
pub fn pocket_colors(n: usize) -> Vec<(i32, i32, i32)> {
    if n == 1 {
        return vec![pocket_colors(2)[0]];
    }
    let (sat, hue_min, hue_max) = (0.6f64, 0.6f64, 1.20f64);
    let hue_range = hue_max - hue_min;
    (0..n)
        .map(|i| {
            let mut hue = hue_min + hue_range * groovy_int_div_to_f64(i as i64, n as i64);
            if hue > 1.0 {
                hue -= 1.0;
            }
            let bright = if i % 2 == 0 { 0.9f32 } else { 0.7f32 };
            hsb_to_rgb(hue as f32, sat as f32, bright)
        })
        .collect()
}

fn color_hex(c: (i32, i32, i32)) -> String {
    format!("#{:02x}{:02x}{:02x}", c.0 & 0xff, c.1 & 0xff, c.2 & 0xff)
}

/// Java Formatter "%5.3f" applied to a BigDecimal (HALF_UP on the decimal value)
fn format_bigdecimal_5_3(dec: &str) -> String {
    let (ip, fp) = dec.split_once('.').unwrap_or((dec, ""));
    let mut digits: Vec<u8> = ip.bytes().chain(fp.bytes()).map(|b| b - b'0').collect();
    let mut int_len = ip.len();
    let frac: Vec<u8> = fp.bytes().map(|b| b - b'0').collect();
    let mut kept: Vec<u8> = digits[..int_len].to_vec();
    for k in 0..3 {
        kept.push(*frac.get(k).unwrap_or(&0));
    }
    let round_up = frac.len() > 3 && frac[3] >= 5;
    if round_up && inc_digits(&mut kept) {
        kept.insert(0, 1);
        int_len += 1;
    }
    digits = kept;
    let ip: String = digits[..int_len].iter().map(|d| (b'0' + d) as char).collect();
    let ip = ip.trim_start_matches('0');
    let ip = if ip.is_empty() { "0" } else { ip };
    let fp: String = digits[int_len..].iter().map(|d| (b'0' + d) as char).collect();
    crate::jfmt::pad_left(&format!("{}.{}", ip, fp), 5)
}

fn py_color(c: (i32, i32, i32)) -> String {
    let f = |v: i32| format_bigdecimal_5_3(&groovy_int_div(v as i64, 255));
    format!("[{},{},{}]", f(c.0), f(c.1), f(c.2))
}

// ------------------------------------------------------------------ CSV

fn sorted_residue_ids(p: &Protein, pocket: &Pocket) -> String {
    let mut set: BTreeSet<String> = BTreeSet::new();
    for &a in &pocket.surface_atoms {
        let g = p.s.atoms[a].group;
        set.insert(p.s.groups[g].resnum.print_full());
    }
    // Java TreeSet<String> order (UTF-16 compareTo); BTreeSet byte order equals it for ASCII
    let mut v: Vec<String> = set.into_iter().collect();
    v.sort_by(|a, b| crate::jutil::java_str_cmp(a, b));
    v.join(" ")
}

pub fn predictions_csv(p: &Protein, pockets: &[Pocket]) -> String {
    let mut sb = String::with_capacity(8192);
    sb.push_str(
        "name, rank, score, probability, sas_points, surf_atoms, center_x, center_y, center_z, residue_ids, surf_atom_ids\n",
    );
    for pk in pockets {
        let mut serials: Vec<i32> = pk.surface_atoms.iter().map(|&a| p.s.atoms[a].serial).collect();
        serials.sort();
        let surf_ids: Vec<String> = serials.iter().map(|s| s.to_string()).collect();
        let row = [
            pk.name.clone(),
            pk.new_rank.to_string(),
            decimal_format(pk.score, 2),
            decimal_format(pk.proba, 3),
            pk.sas_points.len().to_string(),
            pk.surface_atoms.len().to_string(),
            decimal_format(pk.centroid[0], 4),
            decimal_format(pk.centroid[1], 4),
            decimal_format(pk.centroid[2], 4),
            sorted_residue_ids(p, pk),
            surf_ids.join(" "),
        ];
        sb.push_str(&row.join(", "));
        sb.push('\n');
    }
    sb
}

pub fn residues_csv(p: &Protein, rows: &[ResidueRow]) -> String {
    let mut s = String::with_capacity(rows.len() * 48);
    s.push_str("chain, residue_label, residue_name, score, zscore, probability, pocket\n");
    for (r, row) in p.residues.iter().zip(rows) {
        let g = &p.s.groups[r.group];
        s.push_str(&r.chain_author_id);
        s.push_str(", ");
        s.push_str(&g.resnum.to_jstring());
        s.push_str(", ");
        s.push_str(p.residue_code(r).unwrap_or("null"));
        s.push_str(", ");
        s.push_str(&format_f(row.score, 4));
        s.push_str(", ");
        s.push_str(&format_f(row.zscore, 4));
        s.push_str(", ");
        s.push_str(&format_f(row.proba, 4));
        s.push_str(", ");
        s.push_str(&row.pocket.to_string());
        s.push('\n');
    }
    s
}

// ------------------------------------------------------------------ visualizations

/// SAS points PDB text (formatted in parallel chunks, concatenated in order)
pub fn points_pdb(points: &[LabeledPoint], threads: usize) -> String {
    let n = points.len();
    let threads = threads.max(1).min(n.div_ceil(2048).max(1));
    if threads <= 1 {
        return points_pdb_range(points, 0, n);
    }
    let chunk = n.div_ceil(threads);
    let parts: Vec<String> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..threads).map(|t| sc.spawn(move || points_pdb_range(points, t * chunk, ((t + 1) * chunk).min(n)))).collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect()
    });
    parts.concat()
}

fn points_pdb_range(points: &[LabeledPoint], from: usize, to: usize) -> String {
    let mut s = String::with_capacity((to - from) * 68);
    for i in from..to {
        let lp = &points[i];
        s.push_str("HETATM");
        s.push_str(&format_wd(i as i64, 5));
        s.push_str(" H    STP 1  ");
        s.push_str(&format_wd(lp.pocket as i64, 2));
        s.push_str("    ");
        s.push_str(&format_wf(lp.pos[0], 8, 3));
        s.push_str(&format_wf(lp.pos[1], 8, 3));
        s.push_str(&format_wf(lp.pos[2], 8, 3));
        s.push_str("  0.50");
        s.push_str(&format_wf(lp.score, 6, 3));
        s.push('\n');
    }
    s
}

/// Writes gzip exactly like java.util.zip.GZIPOutputStream (default level, zlib deflate)
/// (`level` 6 = Deflater.DEFAULT_COMPRESSION, as in P2Rank)
pub fn write_gzip_java(path: &std::path::Path, data: &[u8], level: u32) -> std::io::Result<()> {
    let mut out: Vec<u8> = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff];
    let mut enc = flate2::write::DeflateEncoder::new(Vec::with_capacity(data.len() / 4), flate2::Compression::new(level));
    enc.write_all(data)?;
    let body = enc.finish()?;
    out.extend_from_slice(&body);
    let mut crc = flate2::Crc::new();
    crc.update(data);
    out.extend_from_slice(&crc.sum().to_le_bytes());
    out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    std::fs::write(path, out)
}

pub fn pymol_script(p: &Protein, pockets: &[Pocket], protein_file: &str, points_file: &str) -> String {
    let n = pockets.len();
    let colors = pocket_colors(n);
    let mut surf = String::new();
    for (k, pk) in pockets.iter().enumerate() {
        let i = k + 1;
        let ids: Vec<String> = pk.surface_atoms.iter().map(|&a| p.s.atoms[a].serial.to_string()).collect();
        surf.push_str(&format!("set_color pcol{} = {}\n", i, py_color(colors[k])));
        surf.push_str(&format!("select surf_pocket{}, protein and id [{}] \n", i, ids.join(",")));
        surf.push_str(&format!("set surface_color,  pcol{}, surf_pocket{} \n", i, i));
    }
    format!(
        r#"from pymol import cmd,stored

set depth_cue, 1
set fog_start, 0.4

set_color b_col, [36,36,85]
set_color t_col, [10,10,10]
set bg_rgb_bottom, b_col
set bg_rgb_top, t_col      
set bg_gradient

set  spec_power  =  200
set  spec_refl   =  0

load "{protein_file}", protein
create ligands, protein and organic
select xlig, protein and organic
delete xlig

hide everything, all

color white, elem c
color bluewhite, protein
#show_as cartoon, protein
show surface, protein
#set transparency, 0.15

show sticks, ligands
set stick_color, magenta






# SAS points

load "{points_file}", points
hide nonbonded, points
show nb_spheres, points
set sphere_scale, 0.2, points
cmd.spectrum("b", "green_red", selection="points", minimum=0, maximum=0.7)


stored.list=[]
cmd.iterate("(resn STP)","stored.list.append(resi)")    # read info about residues STP
lastSTP=stored.list[-1] # get the index of the last residue
hide lines, resn STP

cmd.select("rest", "resn STP and resi 0")

for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.4","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))



{surf}



deselect

orient
"#
    )
}

pub fn chimerax_script(p: &Protein, pockets: &[Pocket], protein_file: &str, points_file: &str) -> String {
    let n = pockets.len();
    let colors: Vec<String> = pocket_colors(n).into_iter().map(color_hex).collect();
    let mut res = String::new();
    for i in 0..n {
        res.push_str(&format!("color name color_pocket{} {}\n", i + 1, colors[i]));
    }
    res.push('\n');
    for (i, pk) in pockets.iter().enumerate() {
        let atoms: Vec<String> = pk.surface_atoms.iter().map(|&a| format!("@@serial_number={}", p.s.atoms[a].serial)).collect();
        res.push_str(&format!("name pocket{}_atoms {}\n", i + 1, atoms.join(" ")));
        res.push_str(&format!("color pocket{}_atoms color_pocket{}\n", i + 1, i + 1));
        res.push('\n');
    }
    format!(
        r#"open {protein_file}
surf
hide solvent
color protein #bdbdde
color ligand magenta


{res}

open {points_file}
color by bfactor #2 palette lime:red  range 0,0.7


set bgcolor #242455
transparency 35
graphics silhouettes true
lighting full
view
"#
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_groovy_div() {
        assert_eq!(groovy_int_div(1, 3), "0.3333333333");
        assert_eq!(groovy_int_div(1, 4), "0.25");
        assert_eq!(groovy_int_div(92, 255), "0.3607843137");
        assert_eq!(format_bigdecimal_5_3("0.3607843137"), "0.361");
        assert_eq!(format_bigdecimal_5_3("1"), "1.000");
    }

    #[test]
    fn test_colors() {
        let c = pocket_colors(4);
        let hex: Vec<String> = c.into_iter().map(color_hex).collect();
        assert_eq!(hex, vec!["#5c93e6", "#7d47b3", "#e65cae", "#b36847"]);
    }
}
