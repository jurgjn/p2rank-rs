//! Port of the legacy `jdk.internal.math.FloatingDecimal.BinaryToASCIIBuffer.dtoa` digit
//! generation and of `java.text.DigitList` rounding, as used by `java.text.DecimalFormat`.
//! DecimalFormat output depends on the exact digits and on the "exact" / "rounded up" flags
//! this algorithm produces, so it is reproduced literally.

use std::cmp::Ordering;

const EXP_SHIFT: i32 = 52;
const FRACT_HOB: u64 = 1u64 << 52;
const EXP_ONE: u64 = 1023u64 << 52;
const SIGNIF_BIT_MASK: u64 = 0x000F_FFFF_FFFF_FFFF;
const EXP_BIT_MASK: u64 = 0x7FF0_0000_0000_0000;
const SIGN_BIT_MASK: u64 = 0x8000_0000_0000_0000;
const MAX_SMALL_BIN_EXP: i32 = 62;
const MIN_SMALL_BIN_EXP: i32 = -(63 / 3);

const N_5_BITS: [i32; 27] = [0, 3, 5, 7, 10, 12, 14, 17, 19, 21, 24, 26, 28, 31, 33, 35, 38, 40, 42, 45, 47, 49, 52, 54, 56, 59, 61];
const INSIGNIFICANT_DIGITS_NUMBER: [i32; 64] = [
    0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 11, 11, 11,
    12, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19,
];

fn long_5_pow(i: i32) -> i64 {
    5i64.pow(i as u32)
}

fn small_5_pow(i: i32) -> i32 {
    5i32.pow(i as u32)
}

const LONG_5_POW_LEN: i32 = 27;

// ------------------------------------------------------------------ minimal big unsigned integer

#[derive(Clone, Debug)]
struct Big(Vec<u32>); // little endian

impl Big {
    fn from_u64(v: u64) -> Big {
        let mut b = Big(vec![v as u32, (v >> 32) as u32]);
        b.trim();
        b
    }
    fn trim(&mut self) {
        while self.0.len() > 1 && *self.0.last().unwrap() == 0 {
            self.0.pop();
        }
        if self.0.is_empty() {
            self.0.push(0);
        }
    }
    fn is_zero(&self) -> bool {
        self.0.iter().all(|&w| w == 0)
    }
    fn mul_small(&mut self, m: u32) {
        let mut carry: u64 = 0;
        for w in self.0.iter_mut() {
            let v = (*w as u64) * (m as u64) + carry;
            *w = v as u32;
            carry = v >> 32;
        }
        if carry > 0 {
            self.0.push(carry as u32);
        }
    }
    fn shl(&mut self, bits: i32) {
        if bits <= 0 {
            return;
        }
        let words = (bits / 32) as usize;
        let rem = (bits % 32) as u32;
        let mut out = vec![0u32; words];
        let mut carry = 0u32;
        for &w in &self.0 {
            if rem == 0 {
                out.push(w);
            } else {
                out.push((w << rem) | carry);
                carry = w >> (32 - rem);
            }
        }
        if carry > 0 {
            out.push(carry);
        }
        self.0 = out;
        self.trim();
    }
    fn pow52(p5: i32, p2: i32) -> Big {
        let mut b = Big::from_u64(1);
        for _ in 0..p5 {
            b.mul_small(5);
        }
        b.shl(p2);
        b
    }
    fn mul_pow52(v: u64, p5: i32, p2: i32) -> Big {
        let mut b = Big::from_u64(v);
        for _ in 0..p5 {
            b.mul_small(5);
        }
        b.shl(p2);
        b
    }
    fn cmp(&self, o: &Big) -> Ordering {
        let (a, b) = (&self.0, &o.0);
        let la = a.iter().rposition(|&w| w != 0).map(|i| i + 1).unwrap_or(0);
        let lb = b.iter().rposition(|&w| w != 0).map(|i| i + 1).unwrap_or(0);
        if la != lb {
            return la.cmp(&lb);
        }
        for i in (0..la).rev() {
            if a[i] != b[i] {
                return a[i].cmp(&b[i]);
            }
        }
        Ordering::Equal
    }
    fn add(&self, o: &Big) -> Big {
        let n = self.0.len().max(o.0.len());
        let mut out = Vec::with_capacity(n + 1);
        let mut carry = 0u64;
        for i in 0..n {
            let v = *self.0.get(i).unwrap_or(&0) as u64 + *o.0.get(i).unwrap_or(&0) as u64 + carry;
            out.push(v as u32);
            carry = v >> 32;
        }
        if carry > 0 {
            out.push(carry as u32);
        }
        let mut b = Big(out);
        b.trim();
        b
    }
    fn sub_assign(&mut self, o: &Big) {
        let mut borrow = 0i64;
        for i in 0..self.0.len() {
            let v = self.0[i] as i64 - *o.0.get(i).unwrap_or(&0) as i64 - borrow;
            if v < 0 {
                self.0[i] = (v + (1i64 << 32)) as u32;
                borrow = 1;
            } else {
                self.0[i] = v as u32;
                borrow = 0;
            }
        }
        self.trim();
    }
    /// q = self / s (q < 10 expected); self = (self % s) * 10
    fn quo_rem_iteration(&mut self, s: &Big) -> i32 {
        let mut q = 0;
        while self.cmp(s) != Ordering::Less {
            self.sub_assign(s);
            q += 1;
        }
        self.mul_small(10);
        q
    }
}

// ------------------------------------------------------------------ dtoa

pub struct Converted {
    pub digits: Vec<u8>,
    pub dec_exponent: i32,
    pub exact: bool,
    pub rounded_up: bool,
}

fn estimate_dec_exp(fract_bits: u64, bin_exp: i32) -> i32 {
    let d2 = f64::from_bits(EXP_ONE | (fract_bits & SIGNIF_BIT_MASK));
    let d = (d2 - 1.5f64) * 0.289529654f64 + 0.176091259 + (bin_exp as f64) * 0.301029995663981;
    let d_bits = d.to_bits();
    let exponent = (((d_bits & EXP_BIT_MASK) >> EXP_SHIFT) as i32) - 1023;
    let is_negative = (d_bits & SIGN_BIT_MASK) != 0;
    if (0..52).contains(&exponent) {
        let mask = SIGNIF_BIT_MASK >> exponent;
        let r = (((d_bits & SIGNIF_BIT_MASK) | FRACT_HOB) >> (EXP_SHIFT - exponent)) as i32;
        if is_negative {
            if (mask & d_bits) == 0 {
                -r
            } else {
                -r - 1
            }
        } else {
            r
        }
    } else if exponent < 0 {
        if (d_bits & !SIGN_BIT_MASK) == 0 {
            0
        } else if is_negative {
            -1
        } else {
            0
        }
    } else {
        d as i32
    }
}

fn insignificant_digits_for_pow2(p2: i32) -> i32 {
    if p2 > 1 && (p2 as usize) < INSIGNIFICANT_DIGITS_NUMBER.len() {
        INSIGNIFICANT_DIGITS_NUMBER[p2 as usize]
    } else {
        0
    }
}

struct Buf {
    digits: Vec<u8>,
    dec_exponent: i32,
    exact: bool,
    rounded_up: bool,
}

impl Buf {
    fn develop_long_digits(&mut self, mut dec_exponent: i32, mut lvalue: i64, insignificant: i32) {
        if insignificant != 0 {
            let pow10 = long_5_pow(insignificant) << insignificant;
            let residue = lvalue % pow10;
            lvalue /= pow10;
            dec_exponent += insignificant;
            if residue >= (pow10 >> 1) {
                lvalue += 1;
            }
        }
        // digits are produced from the least significant end
        let mut rev: Vec<u8> = Vec::new();
        let mut c = (lvalue % 10) as u8;
        lvalue /= 10;
        while c == 0 {
            dec_exponent += 1;
            c = (lvalue % 10) as u8;
            lvalue /= 10;
        }
        while lvalue != 0 {
            rev.push(c + b'0');
            dec_exponent += 1;
            c = (lvalue % 10) as u8;
            lvalue /= 10;
        }
        rev.push(c + b'0');
        rev.reverse();
        self.digits = rev;
        self.dec_exponent = dec_exponent + 1;
    }

    fn roundup(&mut self) {
        let mut i = self.digits.len() - 1;
        let mut q = self.digits[i];
        if q == b'9' {
            while q == b'9' && i > 0 {
                self.digits[i] = b'0';
                i -= 1;
                q = self.digits[i];
            }
            if q == b'9' {
                self.dec_exponent += 1;
                self.digits[0] = b'1';
                return;
            }
        }
        self.digits[i] = q + 1;
        self.rounded_up = true;
    }

    fn dtoa(&mut self, bin_exp: i32, mut fract_bits: u64, n_significant_bits: i32, is_compatible_format: bool) {
        let tail_zeros = fract_bits.trailing_zeros() as i32;
        let n_fract_bits = EXP_SHIFT + 1 - tail_zeros;
        self.rounded_up = false;
        self.exact = false;
        let n_tiny_bits = 0.max(n_fract_bits - bin_exp - 1);
        if (MIN_SMALL_BIN_EXP..=MAX_SMALL_BIN_EXP).contains(&bin_exp)
            && n_tiny_bits < LONG_5_POW_LEN
            && (n_fract_bits + N_5_BITS[n_tiny_bits as usize]) < 64
            && n_tiny_bits == 0
        {
            let insignificant =
                if bin_exp > n_significant_bits { insignificant_digits_for_pow2(bin_exp - n_significant_bits - 1) } else { 0 };
            if bin_exp >= EXP_SHIFT {
                fract_bits <<= bin_exp - EXP_SHIFT;
            } else {
                fract_bits >>= EXP_SHIFT - bin_exp;
            }
            self.develop_long_digits(0, fract_bits as i64, insignificant);
            return;
        }

        let mut dec_exp = estimate_dec_exp(fract_bits, bin_exp);
        let b5 = 0.max(-dec_exp);
        let mut b2 = b5 + n_tiny_bits + bin_exp;
        let s5 = 0.max(dec_exp);
        let mut s2 = s5 + n_tiny_bits;
        let m5 = b5;
        let mut m2 = b2 - n_significant_bits;

        fract_bits >>= tail_zeros;
        b2 -= n_fract_bits - 1;
        let common2 = b2.min(s2);
        b2 -= common2;
        s2 -= common2;
        m2 -= common2;
        if n_fract_bits == 1 {
            m2 -= 1;
        }
        if m2 < 0 {
            b2 -= m2;
            s2 -= m2;
            m2 = 0;
        }

        let mut digits: Vec<u8> = Vec::new();
        let mut low;
        let mut high;

        let b_bits = n_fract_bits + b2 + if (b5 as usize) < N_5_BITS.len() { N_5_BITS[b5 as usize] } else { b5 * 3 };
        let ten_s_bits = s2 + 1 + if ((s5 + 1) as usize) < N_5_BITS.len() { N_5_BITS[(s5 + 1) as usize] } else { (s5 + 1) * 3 };
        if b_bits < 64 && ten_s_bits < 64 {
            if b_bits < 32 && ten_s_bits < 32 {
                let mut b: i32 = ((fract_bits as i32).wrapping_mul(small_5_pow(b5))) << b2;
                let s: i32 = small_5_pow(s5) << s2;
                let mut m: i32 = small_5_pow(m5) << m2;
                let tens: i32 = s * 10;
                let mut q = b / s;
                b = 10i32.wrapping_mul(b % s);
                m = m.wrapping_mul(10);
                low = b < m;
                high = b.wrapping_add(m) > tens;
                if q == 0 && !high {
                    dec_exp -= 1;
                } else {
                    digits.push(b'0' + q as u8);
                }
                if !is_compatible_format || dec_exp < -3 || dec_exp >= 8 {
                    high = false;
                    low = false;
                }
                while !low && !high {
                    q = b / s;
                    b = 10i32.wrapping_mul(b % s);
                    m = m.wrapping_mul(10);
                    if m > 0 {
                        low = b < m;
                        high = b.wrapping_add(m) > tens;
                    } else {
                        low = true;
                        high = true;
                    }
                    digits.push(b'0' + q as u8);
                }
                // Java: (b<<1) - tens evaluated in int arithmetic
                let ldd = (b.wrapping_shl(1)).wrapping_sub(tens) as i64;
                self.exact = b == 0;
                return self.finish(digits, dec_exp, high, low, ldd);
            } else {
                let mut b: i64 = ((fract_bits as i64).wrapping_mul(long_5_pow(b5))) << b2;
                let s: i64 = long_5_pow(s5) << s2;
                let mut m: i64 = long_5_pow(m5) << m2;
                let tens: i64 = s * 10;
                let mut q = (b / s) as i32;
                b = 10i64.wrapping_mul(b % s);
                m = m.wrapping_mul(10);
                low = b < m;
                high = b.wrapping_add(m) > tens;
                if q == 0 && !high {
                    dec_exp -= 1;
                } else {
                    digits.push(b'0' + q as u8);
                }
                if !is_compatible_format || dec_exp < -3 || dec_exp >= 8 {
                    high = false;
                    low = false;
                }
                while !low && !high {
                    q = (b / s) as i32;
                    b = 10i64.wrapping_mul(b % s);
                    m = m.wrapping_mul(10);
                    if m > 0 {
                        low = b < m;
                        high = b.wrapping_add(m) > tens;
                    } else {
                        low = true;
                        high = true;
                    }
                    digits.push(b'0' + q as u8);
                }
                let ldd = (b.wrapping_shl(1)).wrapping_sub(tens);
                self.exact = b == 0;
                return self.finish(digits, dec_exp, high, low, ldd);
            }
        }

        let sval = Big::pow52(s5, s2);
        let mut bval = Big::mul_pow52(fract_bits, b5, b2);
        let mut mval = Big::pow52(m5 + 1, m2 + 1);
        let ten_sval = Big::pow52(s5 + 1, s2 + 1);
        let mut q = bval.quo_rem_iteration(&sval);
        low = bval.cmp(&mval) == Ordering::Less;
        high = ten_sval.cmp(&bval.add(&mval)) != Ordering::Greater;
        if q == 0 && !high {
            dec_exp -= 1;
        } else {
            digits.push(b'0' + q as u8);
        }
        if !is_compatible_format || dec_exp < -3 || dec_exp >= 8 {
            high = false;
            low = false;
        }
        while !low && !high {
            q = bval.quo_rem_iteration(&sval);
            mval.mul_small(10);
            low = bval.cmp(&mval) == Ordering::Less;
            high = ten_sval.cmp(&bval.add(&mval)) != Ordering::Greater;
            digits.push(b'0' + q as u8);
        }
        let ldd: i64 = if high && low {
            bval.shl(1);
            match bval.cmp(&ten_sval) {
                Ordering::Less => -1,
                Ordering::Equal => 0,
                Ordering::Greater => 1,
            }
        } else {
            0
        };
        self.exact = bval.is_zero();
        self.finish(digits, dec_exp, high, low, ldd)
    }

    fn finish(&mut self, digits: Vec<u8>, dec_exp: i32, high: bool, low: bool, low_digit_difference: i64) {
        self.digits = digits;
        self.dec_exponent = dec_exp + 1;
        if high {
            if low {
                if low_digit_difference == 0 {
                    if (self.digits[self.digits.len() - 1] & 1) != 0 {
                        self.roundup();
                    }
                } else if low_digit_difference > 0 {
                    self.roundup();
                }
            } else {
                self.roundup();
            }
        }
    }
}

/// FloatingDecimal.getBinaryToASCIIConverter(d) for a finite, non-zero, non-negative d
pub fn binary_to_ascii(d: f64) -> Converted {
    let d_bits = d.to_bits();
    let mut fract_bits = d_bits & SIGNIF_BIT_MASK;
    let mut bin_exp = ((d_bits & EXP_BIT_MASK) >> EXP_SHIFT) as i32;
    if bin_exp == 0 && fract_bits == 0 {
        return Converted { digits: vec![b'0'], dec_exponent: 0, exact: true, rounded_up: false };
    }
    let n_significant_bits;
    if bin_exp == 0 {
        let leading_zeros = fract_bits.leading_zeros() as i32;
        let shift = leading_zeros - (63 - EXP_SHIFT);
        fract_bits <<= shift;
        bin_exp = 1 - shift;
        n_significant_bits = 64 - leading_zeros;
    } else {
        fract_bits |= FRACT_HOB;
        n_significant_bits = EXP_SHIFT + 1;
    }
    bin_exp -= 1023;
    let mut buf = Buf { digits: Vec::new(), dec_exponent: 0, exact: false, rounded_up: false };
    buf.dtoa(bin_exp, fract_bits, n_significant_bits, true);
    Converted { digits: buf.digits, dec_exponent: buf.dec_exponent, exact: buf.exact, rounded_up: buf.rounded_up }
}

/// BinaryToASCIIBuffer.getChars (Java format string, e.g. "0.0050", "1.0E-5")
fn java_format_string(c: &Converted) -> String {
    let digits = &c.digits;
    let n = digits.len() as i32;
    let de = c.dec_exponent;
    let mut r = String::new();
    if de > 0 && de < 8 {
        let char_len = n.min(de);
        r.push_str(std::str::from_utf8(&digits[..char_len as usize]).unwrap());
        if char_len < de {
            for _ in 0..(de - char_len) {
                r.push('0');
            }
            r.push_str(".0");
        } else {
            r.push('.');
            if char_len < n {
                r.push_str(std::str::from_utf8(&digits[char_len as usize..]).unwrap());
            } else {
                r.push('0');
            }
        }
    } else if de <= 0 && de > -3 {
        r.push_str("0.");
        for _ in 0..(-de) {
            r.push('0');
        }
        r.push_str(std::str::from_utf8(digits).unwrap());
    } else {
        r.push(digits[0] as char);
        r.push('.');
        if n > 1 {
            r.push_str(std::str::from_utf8(&digits[1..]).unwrap());
        } else {
            r.push('0');
        }
        r.push('E');
        let e = if de <= 0 {
            r.push('-');
            -de + 1
        } else {
            de - 1
        };
        r.push_str(&e.to_string());
    }
    r
}

/// java.text.DigitList after set(isNegative, double, maximumDigits, fixedPoint=true)
struct DigitList {
    digits: Vec<u8>,
    count: usize,
    decimal_at: i32,
}

fn should_round_up(dl: &DigitList, maximum_digits: usize, already_rounded: bool, exact: bool) -> bool {
    if maximum_digits < dl.count {
        let d = dl.digits[maximum_digits];
        if d > b'5' {
            return true;
        } else if d == b'5' {
            if maximum_digits == dl.count - 1 {
                if exact {
                    return maximum_digits > 0 && (dl.digits[maximum_digits - 1] % 2 != 0);
                } else {
                    return !already_rounded;
                }
            } else {
                return dl.digits[maximum_digits + 1..dl.count].iter().any(|&c| c != b'0');
            }
        }
    }
    false
}

fn digit_list_set(x: f64, maximum_digits: i32) -> DigitList {
    let conv = binary_to_ascii(x);
    let s = java_format_string(&conv);
    let src = s.as_bytes();
    let mut digits: Vec<u8> = Vec::new();
    let mut decimal_at: i32 = -1;
    let mut exponent: i32 = 0;
    let mut leading_zeros_after_decimal = 0;
    let mut non_zero_seen = false;
    let mut i = 0;
    while i < src.len() {
        let c = src[i];
        i += 1;
        if c == b'.' {
            decimal_at = digits.len() as i32;
        } else if c == b'e' || c == b'E' {
            let rest = std::str::from_utf8(&src[i..]).unwrap();
            exponent = rest.parse().unwrap();
            break;
        } else {
            if !non_zero_seen {
                non_zero_seen = c != b'0';
                if !non_zero_seen && decimal_at != -1 {
                    leading_zeros_after_decimal += 1;
                }
            }
            if non_zero_seen {
                digits.push(c);
            }
        }
    }
    if decimal_at == -1 {
        decimal_at = digits.len() as i32;
    }
    if non_zero_seen {
        decimal_at += exponent - leading_zeros_after_decimal;
    }
    let mut dl = DigitList { count: digits.len(), digits, decimal_at };

    if -dl.decimal_at > maximum_digits {
        dl.count = 0;
        return dl;
    } else if -dl.decimal_at == maximum_digits {
        if should_round_up(&dl, 0, conv.rounded_up, conv.exact) {
            dl.count = 1;
            dl.decimal_at += 1;
            dl.digits[0] = b'1';
        } else {
            dl.count = 0;
        }
        return dl;
    }
    while dl.count > 1 && dl.digits[dl.count - 1] == b'0' {
        dl.count -= 1;
    }
    let max = maximum_digits + dl.decimal_at;
    if max >= 0 && (max as usize) < dl.count {
        let mut maxd = max as usize;
        if should_round_up(&dl, maxd, conv.rounded_up, conv.exact) {
            // roundUp
            let mut m = maxd as i32;
            loop {
                m -= 1;
                if m < 0 {
                    dl.digits[0] = b'1';
                    dl.decimal_at += 1;
                    m = 0;
                    break;
                }
                dl.digits[m as usize] += 1;
                if dl.digits[m as usize] <= b'9' {
                    break;
                }
            }
            maxd = (m + 1) as usize;
        }
        dl.count = maxd;
        while dl.count > 1 && dl.digits[dl.count - 1] == b'0' {
            dl.count -= 1;
        }
    }
    dl
}

/// java.text.DecimalFormat("0." + "0"*frac).format(x) with '.' decimal separator
pub fn decimal_format(x: f64, frac: usize) -> String {
    let neg = x < 0.0 || (x == 0.0 && x.is_sign_negative());
    let v = x.abs();
    let mut s = String::new();
    if neg {
        s.push('-');
    }
    if v.is_nan() {
        return "\u{FFFD}".to_string();
    }
    if v.is_infinite() {
        s.push('\u{221E}');
        return s;
    }
    let dl = if v == 0.0 { DigitList { digits: vec![], count: 0, decimal_at: 0 } } else { digit_list_set(v, frac as i32) };
    let digit = |pos: i32| -> char {
        if pos >= 0 && (pos as usize) < dl.count {
            dl.digits[pos as usize] as char
        } else {
            '0'
        }
    };
    if dl.count == 0 || dl.decimal_at <= 0 {
        s.push('0');
    } else {
        for i in 0..dl.decimal_at {
            s.push(digit(i));
        }
    }
    if frac > 0 {
        s.push('.');
        for i in 0..frac as i32 {
            let pos = if dl.count == 0 { -1 } else { dl.decimal_at + i };
            s.push(digit(pos));
        }
    }
    s
}

/// Double.toString digits are not needed here; exposed for tests.
#[cfg(test)]
pub fn java_double_digits_string(x: f64) -> String {
    java_format_string(&binary_to_ascii(x))
}
