use core::{io::{bse::{ElectronShells, Elements}, PERIODIC_TABLE}, num::convert::aa2au};
use std::collections::BTreeMap;

use crate::cint::NAtom;

use super::{cdata::CintBasis, libcint::CINTgto_norm};

pub trait GetCGTO {
    type OutPut;
    fn to_cgto(&self, kappa_of: i8, norm: bool) -> Self::OutPut;
}

impl GetCGTO for ElectronShells {
    type OutPut = CGTO;
    fn to_cgto(&self, kappa_of: i8, norm: bool) -> Self::OutPut {
        let cgto = CGTO {
            kappa_of,
            angl: self.angular_momentum().clone(),
            exp: self.exponents().iter().map(|x| x.parse().unwrap()).collect(),
            coeff: self
                .coefficients()
                .iter()
                .map(|y| y.iter().map(|x| x.parse().unwrap()).collect())
                .collect(),
        };
        if norm {
            cgto.norm()
        } else {
            cgto
        }
    }
}

impl GetCGTO for Elements {
    type OutPut = Vec<CGTO>;
    fn to_cgto(&self, kappa_of: i8, norm: bool) -> Self::OutPut {
        self.electron_shells()
            .iter()
            .map(|eshl| eshl.to_cgto(kappa_of, norm))
            .collect()
    }
}

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

    pub fn from_xyz(xyz_str: &str, basis: Option<BTreeMap<u8, Vec<CGTO>>>) -> Vec<Self> {
        let mut natm: NAtom = 0;
        let xyz_list: Vec<&str> = xyz_str.trim_end().lines().collect();
        match xyz_list[0].parse() {
            Err(why) => panic!("{:?}", why),
            Ok(v) => natm = v,
        }

        let mut atoms_map: BTreeMap<u8, Vec<[f64; 3]>> = BTreeMap::new();

        let mut counter = 0;
        xyz_list[2..xyz_list.len()].iter().for_each(|line| {
            let mut split_s = line.trim().split_whitespace();
            let nuc = match split_s.next() {
                Some(symbol) => match PERIODIC_TABLE.iter().position(|ele| *ele == symbol) {
                    None => panic!(""),
                    Some(nuc) => (nuc + 1) as u8,
                },
                None => panic!(""),
            };

            let mut coor = [0.0; 3];
            for (i, x_str) in split_s.enumerate() {
                match x_str.parse() {
                    Ok(x) => coor[i] = aa2au(x),
                    Err(why) => panic!(""),
                }
            }

            match atoms_map.get_mut(&nuc) {
                Some(atom_coors) => atom_coors.push(coor),
                None => {
                    atoms_map.insert(nuc, vec![coor]);
                    ()
                }
            }

            counter += 1;
        });
        assert_eq!(natm, counter);

        atoms_map
            .iter()
            .map(|(nuc, coors)| {
                Self::new(
                    match &basis {
                        Some(bas) => match bas.get(&nuc) {
                            Some(b) => Some(b.to_vec()),
                            None => panic!(""),
                        },
                        None => None,
                    },
                    *nuc,
                    0,
                    0.0,
                    0.0,
                    coors.to_vec(),
                )
            })
            .collect()
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
