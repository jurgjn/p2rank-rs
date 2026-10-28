//! Number formatting reproducing Java semantics:
//! - `java.text.DecimalFormat("0.000")`: HALF_EVEN rounding of the exact binary value
//! - `java.util.Formatter` `%.Nf` / `%W.Nf`: HALF_UP rounding of the shortest repr digits
//!   (the `Double.toString` digits), `-` sign for negative values and -0.0

/// Shortest round-trip decimal digits of a finite non-negative double, with the decimal exponent
/// in Java's convention (value = 0.d1d2... * 10^exp).
fn shortest_digits(v: f64) -> (Vec<u8>, i32) {
    if v == 0.0 {
        return (vec![b'0'], 1);
    }
    // ryu: shortest round-trip digits, e.g. "1.0005", "123.0", "1e-7", "1.5e20"
    let mut buf = ryu::Buffer::new();
    let s = buf.format_finite(v).as_bytes();
    let mut digits: Vec<u8> = Vec::with_capacity(20);
    let mut point: Option<i32> = None;
    let mut exp: i32 = 0;
    let mut i = 0;
    while i < s.len() {
        match s[i] {
            b'.' => point = Some(digits.len() as i32),
            b'e' | b'E' => {
                exp = std::str::from_utf8(&s[i + 1..]).unwrap().parse().unwrap();
                break;
            }
            c => digits.push(c),
        }
        i += 1;
    }
    let int_len = point.unwrap_or(digits.len() as i32);
    // strip leading zeros (adjusting the decimal exponent) and trailing zeros
    let lead = digits.iter().take_while(|&&c| c == b'0').count();
    digits.drain(..lead);
    while digits.len() > 1 && *digits.last().unwrap() == b'0' {
        digits.pop();
    }
    (digits, int_len - lead as i32 + exp)
}

fn is_negative(x: f64) -> bool {
    x < 0.0 || (x == 0.0 && x.is_sign_negative())
}

/// Increment a decimal digit string by one unit in the last place; returns true on overflow (all 9s).
fn increment(d: &mut [u8]) -> bool {
    for i in (0..d.len()).rev() {
        if d[i] == b'9' {
            d[i] = b'0';
        } else {
            d[i] += 1;
            return false;
        }
    }
    true
}

/// DecimalFormat with pattern "0." + "0"*frac
pub fn decimal_format(x: f64, frac: usize) -> String {
    crate::jdtoa::decimal_format(x, frac)
}

/// java.util.Formatter "%.{prec}f" (no width)
pub fn format_f(x: f64, prec: usize) -> String {
    if x.is_nan() {
        return "NaN".to_string();
    }
    let neg = is_negative(x);
    let v = x.abs();
    if v.is_infinite() {
        return if neg { "-Infinity".to_string() } else { "Infinity".to_string() };
    }
    let (mut digits, mut dec_exp) = shortest_digits(v);
    let n = digits.len() as i32;
    let keep = dec_exp + prec as i32;
    if keep < 0 {
        digits = vec![b'0'];
        dec_exp = 1;
    } else if keep < n {
        if keep == 0 {
            if digits[0] >= b'5' {
                digits = vec![b'1'];
                dec_exp += 1;
            } else {
                digits = vec![b'0'];
                dec_exp = 1;
            }
        } else {
            let k = keep as usize;
            let round_up = digits[k] >= b'5';
            digits.truncate(k);
            if round_up && increment(&mut digits) {
                digits = vec![b'1'];
                dec_exp += 1;
            }
        }
    }
    // render digits * 10^(dec_exp - len) with prec fractional digits
    let mut int_part = String::new();
    let mut frac_part = String::new();
    let len = digits.len() as i32;
    for i in 0..dec_exp.max(0) {
        int_part.push(if i < len { digits[i as usize] as char } else { '0' });
    }
    if int_part.is_empty() {
        int_part.push('0');
    }
    for k in 0..prec as i32 {
        let pos = dec_exp + k; // index into digits
        frac_part.push(if pos >= 0 && pos < len { digits[pos as usize] as char } else { '0' });
    }
    let mut s = String::new();
    if neg {
        s.push('-');
    }
    // strip leading zeros of int part
    let ip = int_part.trim_start_matches('0');
    s.push_str(if ip.is_empty() { "0" } else { ip });
    if prec > 0 {
        s.push('.');
        s.push_str(&frac_part);
    }
    s
}

/// java.util.Formatter "%{width}.{prec}f"
pub fn format_wf(x: f64, width: usize, prec: usize) -> String {
    let s = format_f(x, prec);
    pad_left(&s, width)
}

pub fn pad_left(s: &str, width: usize) -> String {
    let n = s.chars().count();
    if n >= width {
        s.to_string()
    } else {
        format!("{}{}", " ".repeat(width - n), s)
    }
}

/// java.util.Formatter "%{width}d"
pub fn format_wd(x: i64, width: usize) -> String {
    pad_left(&x.to_string(), width)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_against_java_table() {
        // Java reference output (tools/java/FmtTest.java); P2RS_FMT_TABLE may point to the full table
        let text = match std::env::var("P2RS_FMT_TABLE") {
            Ok(p) => std::fs::read_to_string(p).unwrap(),
            Err(_) => {
                use std::io::Read;
                let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/java_number_format.txt.gz");
                let mut s = String::new();
                flate2::read::GzDecoder::new(std::fs::File::open(path).unwrap()).read_to_string(&mut s).unwrap();
                s
            }
        };
        let mut n = 0;
        for line in text.lines() {
            let f: Vec<&str> = line.split('|').collect();
            let x = f64::from_bits(u64::from_str_radix(f[0], 16).unwrap());
            assert_eq!(decimal_format(x, 2), f[1], "0.00 {}", f[8]);
            assert_eq!(decimal_format(x, 3), f[2], "0.000 {}", f[8]);
            assert_eq!(decimal_format(x, 4), f[3], "0.0000 {}", f[8]);
            assert_eq!(format_f(x, 4), f[4], "%.4f {}", f[8]);
            assert_eq!(format_wf(x, 8, 3), f[5], "%8.3f {}", f[8]);
            assert_eq!(format_wf(x, 6, 3), f[6], "%6.3f {}", f[8]);
            assert_eq!(format_f(x, 3), f[7], "%.3f {}", f[8]);
            n += 1;
        }
        assert!(n > 1000);
    }

    #[test]
    fn test_basic() {
        assert_eq!(decimal_format(0.125, 2), "0.12");
        assert_eq!(decimal_format(-0.0, 2), "-0.00");
        assert_eq!(format_f(1.0005, 3), "1.001");
        assert_eq!(format_f(9.9995, 3), "10.000");
        assert_eq!(format_wf(-0.0, 8, 3), "  -0.000");
    }
}
