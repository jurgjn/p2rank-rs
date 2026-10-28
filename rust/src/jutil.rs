//! Helpers reproducing Java standard library semantics (string slicing, number parsing).

/// A line of text with Java (UTF-16 code unit) indexing semantics.
pub struct JLine {
    s: String,
    /// UTF-16 code units, only materialized for non-ASCII lines
    utf16: Option<Vec<u16>>,
}

#[derive(Debug, Clone, Copy)]
pub struct IndexOutOfBounds;

impl JLine {
    pub fn new(s: String) -> JLine {
        let utf16 = if s.is_ascii() { None } else { Some(s.encode_utf16().collect()) };
        JLine { s, utf16 }
    }

    pub fn as_str(&self) -> &str {
        &self.s
    }

    /// Java String.length()
    pub fn len(&self) -> usize {
        match &self.utf16 {
            None => self.s.len(),
            Some(u) => u.len(),
        }
    }

    /// Java String.substring(begin, end)
    pub fn sub(&self, begin: usize, end: usize) -> Result<String, IndexOutOfBounds> {
        if begin > end || end > self.len() {
            return Err(IndexOutOfBounds);
        }
        match &self.utf16 {
            None => Ok(self.s[begin..end].to_string()),
            Some(u) => Ok(String::from_utf16_lossy(&u[begin..end])),
        }
    }

    /// Java String.substring(begin)
    pub fn sub_from(&self, begin: usize) -> Result<String, IndexOutOfBounds> {
        self.sub(begin, self.len())
    }

    pub fn starts_with(&self, p: &str) -> bool {
        self.s.starts_with(p)
    }
}

/// Java String.trim(): strips all chars <= ' ' from both ends.
pub fn jtrim(s: &str) -> &str {
    s.trim_matches(|c: char| c <= ' ')
}

/// Java Character.isWhitespace for chars used in split("\\s+") / StringTokenizer contexts.
pub fn is_java_space(c: char) -> bool {
    // regex \s in Java = [ \t\n\x0B\f\r]
    matches!(c, ' ' | '\t' | '\n' | '\u{0B}' | '\u{0C}' | '\r')
}

/// Splits a byte buffer into lines following java.io.BufferedReader.readLine() rules
/// (terminators: "\n", "\r", "\r\n"). Input is decoded as UTF-8 (lossy).
pub fn read_lines(bytes: &[u8]) -> Vec<String> {
    let text = String::from_utf8_lossy(bytes);
    let mut lines = Vec::new();
    let b = text.as_bytes();
    let mut start = 0;
    let mut i = 0;
    while i < b.len() {
        let c = b[i];
        if c == b'\n' || c == b'\r' {
            lines.push(text[start..i].to_string());
            if c == b'\r' && i + 1 < b.len() && b[i + 1] == b'\n' {
                i += 1;
            }
            start = i + 1;
        }
        i += 1;
    }
    if start < b.len() {
        lines.push(text[start..].to_string());
    }
    lines
}

fn is_decimal_float_literal(s: &str) -> bool {
    // [+-]? (digits [. digits?] | . digits) ([eE][+-]?digits)?
    let b = s.as_bytes();
    let mut i = 0;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let mut int_digits = 0;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
        int_digits += 1;
    }
    let mut frac_digits = 0;
    if i < b.len() && b[i] == b'.' {
        i += 1;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
            frac_digits += 1;
        }
    }
    if int_digits + frac_digits == 0 {
        return false;
    }
    if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
        i += 1;
        if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
            i += 1;
        }
        let mut exp_digits = 0;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
            exp_digits += 1;
        }
        if exp_digits == 0 {
            return false;
        }
    }
    i == b.len()
}

/// Prepares a string for Java-style float parsing: trims, handles the optional
/// [fFdD] suffix and the NaN/Infinity literals. Returns the plain literal for Rust parsing.
fn java_float_literal(s: &str) -> Option<String> {
    let t = jtrim(s);
    if t.is_empty() {
        return None;
    }
    let (sign, rest) = match t.as_bytes()[0] {
        b'+' => ("", &t[1..]),
        b'-' => ("-", &t[1..]),
        _ => ("", t),
    };
    if rest == "NaN" {
        return Some("NaN".to_string());
    }
    if rest == "Infinity" {
        return Some(format!("{}inf", sign));
    }
    let mut body = t;
    if let Some(last) = body.chars().last() {
        if matches!(last, 'f' | 'F' | 'd' | 'D') {
            body = &body[..body.len() - 1];
        }
    }
    if is_decimal_float_literal(body) {
        Some(body.to_string())
    } else {
        None
    }
}

/// java.lang.Double.parseDouble (decimal literals; hex literals are not supported and rejected)
pub fn parse_java_double(s: &str) -> Option<f64> {
    let lit = java_float_literal(s)?;
    lit.parse::<f64>().ok()
}

/// java.lang.Float.parseFloat (decimal literals; hex literals are not supported and rejected)
pub fn parse_java_float(s: &str) -> Option<f32> {
    let lit = java_float_literal(s)?;
    lit.parse::<f32>().ok()
}

/// java.lang.Integer.parseInt (ASCII digits only)
pub fn parse_java_int(s: &str) -> Option<i32> {
    let b = s.as_bytes();
    if b.is_empty() {
        return None;
    }
    let (neg, digits) = match b[0] {
        b'+' => (false, &s[1..]),
        b'-' => (true, &s[1..]),
        _ => (false, s),
    };
    if digits.is_empty() || !digits.bytes().all(|c| c.is_ascii_digit()) {
        return None;
    }
    let mut v: i64 = 0;
    for c in digits.bytes() {
        v = v * 10 + (c - b'0') as i64;
        if v > (i32::MAX as i64) + 1 {
            return None;
        }
    }
    let v = if neg { -v } else { v };
    if v < i32::MIN as i64 || v > i32::MAX as i64 {
        return None;
    }
    Some(v as i32)
}

/// Java's String.compareTo (UTF-16 code unit order)
pub fn java_str_cmp(a: &str, b: &str) -> std::cmp::Ordering {
    if a.is_ascii() && b.is_ascii() {
        return a.cmp(b);
    }
    a.encode_utf16().cmp(b.encode_utf16())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_double() {
        assert_eq!(parse_java_double(" 12.5 "), Some(12.5));
        assert_eq!(parse_java_double("1."), Some(1.0));
        assert_eq!(parse_java_double(".5"), Some(0.5));
        assert_eq!(parse_java_double("1e3"), Some(1000.0));
        assert_eq!(parse_java_double("1.5d"), Some(1.5));
        assert_eq!(parse_java_double("-Infinity"), Some(f64::NEG_INFINITY));
        assert!(parse_java_double("inf").is_none());
        assert!(parse_java_double(".").is_none());
        assert!(parse_java_double("").is_none());
        assert!(parse_java_double("1e").is_none());
        assert!(parse_java_double("abc").is_none());
    }

    #[test]
    fn test_parse_int() {
        assert_eq!(parse_java_int("+12"), Some(12));
        assert_eq!(parse_java_int("-2147483648"), Some(i32::MIN));
        assert_eq!(parse_java_int("2147483648"), None);
        assert_eq!(parse_java_int(" 1"), None);
    }

    #[test]
    fn test_lines() {
        assert_eq!(read_lines(b"a\r\nb\rc\n\nd"), vec!["a", "b", "c", "", "d"]);
    }
}
