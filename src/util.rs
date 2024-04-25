pub trait AbsMax {
    type Output;
    fn abs(self) -> Self::Output;
    const MAX: Self::Output;
}

impl AbsMax for i8 {
    type Output = i8;
    fn abs(self) -> Self::Output {
        i8::abs(self)
    }
    const MAX: Self::Output = i8::MAX;
}

impl AbsMax for i16 {
    type Output = i16;
    fn abs(self) -> Self::Output {
        i16::abs(self)
    }
    const MAX: Self::Output = i16::MAX;
}

impl AbsMax for i32 {
    type Output = i32;
    fn abs(self) -> Self::Output {
        i32::abs(self)
    }
    const MAX: Self::Output = i32::MAX;
}

impl AbsMax for i64 {
    type Output = i64;
    fn abs(self) -> Self::Output {
        i64::abs(self)
    }
    const MAX: Self::Output = i64::MAX;
}

impl AbsMax for i128 {
    type Output = i128;
    fn abs(self) -> Self::Output {
        i128::abs(self)
    }
    const MAX: Self::Output = i128::MAX;
}

impl AbsMax for u8 {
    type Output = u8;
    fn abs(self) -> Self::Output {
        u8::abs(self)
    }
    const MAX: Self::Output = u8::MAX;
}

impl AbsMax for u16 {
    type Output = u16;
    fn abs(self) -> Self::Output {
        u16::abs(self)
    }
    const MAX: Self::Output = u16::MAX;
}

impl AbsMax for u32 {
    type Output = u32;
    fn abs(self) -> Self::Output {
        u32::abs(self)
    }
    const MAX: Self::Output = u32::MAX;
}

impl AbsMax for u64 {
    type Output = u64;
    fn abs(self) -> Self::Output {
        u64::abs(self)
    }
    const MAX: Self::Output = u64::MAX;
}

impl AbsMax for u128 {
    type Output = u128;
    fn abs(self) -> Self::Output {
        u128::abs(self)
    }
    const MAX: Self::Output = u128::MAX;
}

pub(crate) trait NumericElement:
    From<i8>
    + std::ops::Rem<Output = Self>
    + std::ops::Sub<Output = Self>
    + std::ops::Add<Output = Self>
    + std::ops::Div<Output = Self>
    + std::cmp::Ord
    + std::ops::Mul<Output = Self>
    + std::fmt::Display
    + std::ops::RemAssign
    + std::ops::AddAssign
    + Copy
    + AbsMax<Output = Self>
{
}

impl<T> NumericElement for T where
    T: From<i8>
        + std::ops::Rem<Output = Self>
        + std::ops::Sub<Output = Self>
        + std::ops::Add<Output = Self>
        + std::ops::Div<Output = Self>
        + std::cmp::Ord
        + std::ops::Mul<Output = Self>
        + std::fmt::Display
        + std::ops::RemAssign
        + std::ops::AddAssign
        + Copy
        + AbsMax<Output = Self>
{
}

/// Find the greatest common divisor.
fn gcd<T>(mut n: T, mut m: T) -> Result<T, &'static str>
where
    T: NumericElement,
{
    if n <= T::from(0) || m <= T::from(0) {
        return Err("zero or negative values not supported");
    }
    while m != T::from(0) {
        if m < n {
            std::mem::swap(&mut m, &mut n);
        }
        m = m % n;
    }
    Ok(n)
}

/// This is a brute-force implementation of modular inverse. The Extended Euclidian Algorithm might be a better choice.
fn meziriac<T>(a: T, b: T) -> Result<T, &'static str>
where
    T: NumericElement,
{
    let mut g = T::from(1);
    if b == T::from(1) {
        g = T::from(1);
    } else if a == b {
        g = T::from(0);
    } else {
        while g < T::MAX {
            if ((g * a) % b) == T::from(1) {
                break;
            }
            g += T::from(1);
        }
    }
    Ok(g)
}

/// Core implementation of intersection of two residual classes.
pub(crate) fn intersection<T>(m1: T, m2: T, mut s1: T, mut s2: T) -> Result<(T, T), &'static str>
where
    T: NumericElement,
{
    if m1 == T::from(0) || m2 == T::from(0) {
        // intersection of null and anything is null
        return Ok((T::from(0), T::from(0)));
    }
    // normalize shifts
    s1 %= m1;
    s2 %= m2;

    // use common divisor
    let d = gcd(m1, m2)?;
    let md1 = m1 / d;
    let md2 = m2 / d;
    let span: T = (s2 - s1).abs();

    if d != T::from(1) && (span % d != T::from(0)) {
        return Ok((T::from(0), T::from(0))); // no intersection
    }
    // NOTE: though this case was specified, it seems impossible to replicate
    // if d != 1 && (span % d == 0) && (s1 != s2) && (md1 == md2) {
    //     return Ok((d, s1));
    // }

    // d might be 1
    let m = md1 * md2 * d;
    Ok((m, (s1 + (meziriac(md1, md2).unwrap() * span * md1)) % m))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gcd_a() {
        assert_eq!(gcd(14, 15).unwrap(), 1);
    }

    #[test]
    fn test_gcd_b() {
        assert_eq!(gcd(12, 8).unwrap(), 4);
    }

    #[test]
    fn test_gcd_c() {
        let a = 2 * 3 * 5 * 11 * 17;
        let b = 3 * 7 * 11 * 13 * 19;
        assert_eq!(gcd(a, b).unwrap(), 3 * 11);
    }

    #[test]
    fn test_gcd_d() {
        assert_eq!(gcd(12, 0).is_err(), true);
    }

    #[test]
    fn test_gcd_e() {
        assert_eq!(gcd(0, 3).is_err(), true);
    }

    #[test]
    fn test_intersection_a() {
        assert_eq!(intersection(0, 0, 2, 3).unwrap(), (0, 0));
    }

    #[test]
    fn test_intersection_b() {
        assert_eq!(intersection(45, 40, 11, 1).unwrap(), (360, 101));
    }

    #[test]
    fn test_meziriac_a() {
        assert_eq!(meziriac(1, 1).unwrap(), 1);
        assert_eq!(meziriac(10, 1).unwrap(), 1);
        assert_eq!(meziriac(10, 10).unwrap(), 0);
        assert_eq!(meziriac(12, 12).unwrap(), 0);
        assert_eq!(meziriac(3, 11).unwrap(), 4);
        assert_eq!(meziriac(20, 9).unwrap(), 5);
        assert_eq!(meziriac(101, 13).unwrap(), 4);
    }
}
