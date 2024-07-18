use super::{cdata::CintBasis, libccint::CINTgto_norm};

#[derive(Debug, Clone)]
pub struct CGTO {
    pub kappa_of: i8,
    pub angl: Vec<u8>,
    pub exp: Vec<f64>,
    pub coeff: Vec<Vec<f64>>,
}

impl CGTO {
    pub fn norm(mut self) -> Self {
        // pub fn CINTgto_norm(n: ::std::os::raw::c_int, a: f64) -> f64;
        self.coeff = self
            .coeff
            .iter()
            .enumerate()
            .map(|(i, coeff)| {
                coeff
                    .iter()
                    .enumerate()
                    .map(|(ic, c)| c * unsafe { CINTgto_norm(self.angl[i].into(), self.exp[ic]) })
                    .collect()
            })
            .collect();
        self
    }

    pub fn gen_bas(&self, iangl: usize, ptr_exp: i32, ptr_coeff: i32) -> CintBasis {
        CintBasis {
            atom_of: -1,
            ang_of: self.angl[iangl].into(),
            nprim_of: self.exp.len() as i32,
            nctr_of: 1,
            kappa_of: self.kappa_of as i32,
            ptr_exp,
            ptr_coeff,
            reserve_baslot: 0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct CintAtomGroup {
    basis: Option<Vec<CGTO>>,
    charge_of: u8,
    nuc_mod_of: u8,
    zeta: f64,
    frac_charge: f64,
    coordinates: Vec<[f64; 3]>,
}

impl CintAtomGroup {
    pub fn new(
        basis: Option<Vec<CGTO>>,
        charge_of: u8,
        nuc_mod_of: u8,
        zeta: f64,
        frac_charge: f64,
        coordinates: Vec<[f64; 3]>,
    ) -> Self {
        Self {
            basis,
            charge_of,
            nuc_mod_of,
            zeta,
            frac_charge,
            coordinates,
        }
    }
    pub fn from_other_group(atom_group: &impl AtomGroup) -> Self {
        Self {
            basis: atom_group.basis().clone(),
            charge_of: *atom_group.charge_of(),
            nuc_mod_of: *atom_group.nuc_mod_of(),
            zeta: *atom_group.zeta(),
            frac_charge: *atom_group.frac_charge(),
            coordinates: atom_group.coordinates().clone(),
        }
    }
}

impl AtomGroup for CintAtomGroup {
    fn basis_mut(&mut self) -> &mut Option<Vec<CGTO>> {
        &mut self.basis
    }

    fn charge_of_mut(&mut self) -> &mut u8 {
        &mut self.charge_of
    }

    fn nuc_mod_of_mut(&mut self) -> &mut u8 {
        &mut self.nuc_mod_of
    }

    fn zeta_mut(&mut self) -> &mut f64 {
        &mut self.zeta
    }

    fn frac_charge_mut(&mut self) -> &mut f64 {
        &mut self.frac_charge
    }

    fn coordinates_mut(&mut self) -> &mut Vec<[f64; 3]> {
        &mut self.coordinates
    }

    fn basis(&self) -> &Option<Vec<CGTO>> {
        &self.basis
    }

    fn charge_of(&self) -> &u8 {
        &self.charge_of
    }

    fn nuc_mod_of(&self) -> &u8 {
        &self.nuc_mod_of
    }

    fn zeta(&self) -> &f64 {
        &self.zeta
    }

    fn frac_charge(&self) -> &f64 {
        &self.frac_charge
    }

    fn coordinates(&self) -> &Vec<[f64; 3]> {
        &self.coordinates
    }
}

pub trait AtomGroup {
    fn basis_mut(&mut self) -> &mut Option<Vec<CGTO>>;
    fn charge_of_mut(&mut self) -> &mut u8;
    fn nuc_mod_of_mut(&mut self) -> &mut u8;
    fn zeta_mut(&mut self) -> &mut f64;
    fn frac_charge_mut(&mut self) -> &mut f64;
    fn coordinates_mut(&mut self) -> &mut Vec<[f64; 3]>;

    fn basis(&self) -> &Option<Vec<CGTO>>;
    fn charge_of(&self) -> &u8;
    fn nuc_mod_of(&self) -> &u8;
    fn zeta(&self) -> &f64;
    fn frac_charge(&self) -> &f64;
    fn coordinates(&self) -> &Vec<[f64; 3]>;
}
