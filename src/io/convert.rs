pub const AUTOAA: f64 = 0.52917726;
pub const AATOAU: f64 = 1.0 / AUTOAA;

pub fn autoaa(x: f64) -> f64 {
    x * (0.52917726)
}

pub fn aatoau(x: f64) -> f64 {
    x * (1.0 / 0.52917726)
}

pub const AMUTOKG: f64 = 1.660539040e-27;
pub const KGTOAMU: f64 = 1.0 / AMUTOKG;
pub const METOKG: f64 = 9.10938356e-31;
pub const KGTOME: f64 = 1.0 / METOKG;
pub const AMUTOAU: f64 = AMUTOKG * KGTOME;
