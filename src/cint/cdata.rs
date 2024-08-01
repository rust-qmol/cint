use core::io::bse::JsonBasis;
use std::{collections::BTreeMap, fs};


use super::{
    intor::Intor,
    libcint::{
        ANG_MAX, NGRIDS, PTR_COMMON_ORIG, PTR_ENV_START, PTR_EXPCUTOFF, PTR_F12_ZETA, PTR_GRIDS,
        PTR_GTG_ZETA, PTR_RANGE_OMEGA, PTR_RINV_ORIG, PTR_RINV_ZETA,
    },
    rawdata::{AtomGroup, CintAtomGroup, GetCGTO},
    AtomIndex, BasisIndex, NAtom,
};

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct CintAtom {
    charge_of: i32,
    coord: i32,
    nuc_mod_of: i32,
    zeta: i32,
    frac_charge: i32,
    reserve_atmslot: i32,
}

impl CintAtom {
    pub fn empty() -> Self {
        CintAtom {
            charge_of: 0,
            coord: 0,
            nuc_mod_of: 0,
            zeta: 0,
            frac_charge: 0,
            reserve_atmslot: 0,
        }
    }
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct CintBasis {
    pub(super) atom_of: i32,
    pub(super) ang_of: i32,
    pub(super) nprim_of: i32,
    pub(super) nctr_of: i32,
    pub(super) kappa_of: i32,
    pub(super) ptr_exp: i32,
    pub(super) ptr_coeff: i32,
    pub(super) reserve_baslot: i32,
}

impl CintBasis {
    pub fn empty() -> Self {
        Self {
            atom_of: 0,
            ang_of: 0,
            nprim_of: 0,
            nctr_of: 0,
            kappa_of: 0,
            ptr_exp: 0,
            ptr_coeff: 0,
            reserve_baslot: 0,
        }
    }
}

const ATM_OFFSET: usize = PTR_ENV_START as usize;
const ATM_SLOT: usize = 4;

#[repr(C)]
#[derive(Debug, Clone)]
pub struct CintEnv {
    pub(self) data: Vec<f64>,
}

impl CintEnv {
    pub fn expcutoff(&mut self) -> &mut f64 {
        &mut self.data[PTR_EXPCUTOFF as usize]
    }
    pub fn common_orig(&mut self) -> &mut [f64] {
        &mut self.data[(PTR_COMMON_ORIG as usize)..(PTR_RINV_ORIG as usize)]
    }
    pub fn rinv_orig(&mut self) -> &mut [f64] {
        &mut self.data[(PTR_RINV_ORIG as usize)..(PTR_RINV_ZETA as usize)]
    }
    pub fn rinv_zeta(&mut self) -> &mut f64 {
        &mut self.data[PTR_RANGE_OMEGA as usize]
    }
    pub fn range_omega(&mut self) -> &mut f64 {
        &mut self.data[PTR_RANGE_OMEGA as usize]
    }
    pub fn f12_zeta(&mut self) -> &mut f64 {
        &mut self.data[PTR_F12_ZETA as usize]
    }
    pub fn gtg_zeta(&mut self) -> &mut f64 {
        &mut self.data[PTR_GTG_ZETA as usize]
    }
    pub fn ngrids(&mut self) -> &mut f64 {
        &mut self.data[NGRIDS as usize]
    }
    pub fn grids(&mut self) -> &mut [f64] {
        &mut self.data[(PTR_GRIDS as usize)..(PTR_ENV_START as usize)]
    }
    pub fn atom(&mut self, iatm: usize) -> &mut [f64] {
        let start = ATM_OFFSET + ATM_SLOT * iatm;
        let end = ATM_OFFSET + ATM_SLOT * (iatm + 1);
        &mut self.data[start..end]
    }
    pub fn as_ptr(&self) -> *const f64 {
        self.data.as_ptr()
    }
}

#[derive(Debug, Clone)]
pub struct CintDate {
    atom_groups: Vec<CintAtomGroup>,
    basis_template: Vec<Vec<CintBasis>>,
    atmbas_index: Vec<AtomIndex>,
    atm: Vec<CintAtom>,
    pub env: CintEnv,
}

impl CintDate {
    pub fn new(atom_groups: Vec<impl AtomGroup>) -> Self {
        atom_groups.iter().for_each(|atoms| match atoms.basis() {
            Some(_) => {}
            None => panic!(),
        });

        let natm: NAtom = atom_groups
            .iter()
            .map(|atom| atom.coordinates().len())
            .sum();

        //
        let mut ptr_exp = (PTR_ENV_START as usize) + natm * ATM_SLOT;
        let mut ptr_coeff = ptr_exp;
        let basis_template: Vec<Vec<CintBasis>> = atom_groups
            .iter()
            .map(|atoms| match &atoms.basis() {
                Some(basis) => basis
                    .iter()
                    .flat_map(|cgto| {
                        ptr_coeff += cgto.exp.len();
                        let atom_basis: Vec<CintBasis> = cgto
                            .angl
                            .iter()
                            .enumerate()
                            .map(|(i, angl)| {
                                assert!(*angl < ANG_MAX as u8);
                                let cint_basis = cgto.gen_bas(i, ptr_exp as i32, ptr_coeff as i32);
                                ptr_coeff += cgto.exp.len();
                                cint_basis
                            })
                            .collect();
                        ptr_exp += cgto.exp.len() * (cgto.coeff.len() + 1);
                        atom_basis
                    })
                    .collect(),
                None => panic!(""),
            })
            .collect();

        //
        let atmbas_index = atom_groups
            .iter()
            .enumerate()
            .flat_map(|(i, atoms)| vec![i; atoms.coordinates().len()])
            .collect();

        //
        let mut iatm = 0;
        let atm: Vec<CintAtom> = atom_groups
            .iter()
            .flat_map(|atoms| {
                atoms
                    .coordinates()
                    .iter()
                    .map(|coor| {
                        let coord_i = (PTR_ENV_START as usize) + iatm * ATM_SLOT;
                        let zeta_i = coord_i + 3;
                        iatm += 1;
                        CintAtom {
                            charge_of: (*atoms.charge_of()).into(),
                            coord: coord_i as i32,
                            nuc_mod_of: (*atoms.nuc_mod_of()).into(),
                            zeta: zeta_i as i32,
                            frac_charge: (*atoms.frac_charge()) as i32,
                            reserve_atmslot: 0,
                        }
                    })
                    .collect::<Vec<CintAtom>>()
            })
            .collect();

        //
        let env_atm = atom_groups
            .iter()
            .flat_map(|atoms| {
                atoms
                    .coordinates()
                    .iter()
                    .flat_map(|coor| vec![coor[0], coor[1], coor[2], *atoms.zeta()])
            })
            .collect();

        let env_bas = atom_groups
            .iter()
            .flat_map(|atoms| match &atoms.basis() {
                Some(basis) => basis.iter().flat_map(|cgto| {
                    cgto.exp
                        .iter()
                        .map(|x| *x)
                        .chain(cgto.coeff.iter().flat_map(|coeff| coeff.iter().map(|x| *x)))
                }),
                None => panic!(""),
            })
            .collect();

        let env = CintEnv {
            data: vec![vec![0.0_f64; PTR_ENV_START as usize], env_atm, env_bas].concat(),
        };

        //
        Self {
            atom_groups: atom_groups
                .iter()
                .map(|atom| CintAtomGroup::from_other_group(atom))
                .collect(),
            basis_template,
            atmbas_index,
            atm,
            env,
        }
    }

    pub fn fron_xyz(xyz_str: &str, basis_path: &str) -> Self {
        let basis_str = fs::read_to_string(basis_path).expect("Error in reading the basis file");
        let json: JsonBasis =
            serde_json::from_str(&basis_str).expect("read basis from json failed");

        let mut atom_group = CintAtomGroup::from_xyz(xyz_str, None);
        atom_group
            .iter_mut()
            .for_each(|atoms| match json.get_elements(*atoms.charge_of()) {
                Some(ele) => (*atoms.basis_mut()) = Some(ele.to_cgto(0, true)),
                None => panic!(""),
            });

        CintDate::new(atom_group)
    }

    fn gen_bas_select(&self, which_bas: BTreeMap<AtomIndex, Vec<BasisIndex>>) -> Vec<CintBasis> {
        which_bas
            .iter()
            .flat_map(|(iatm, ibas_all)| {
                ibas_all.iter().flat_map(|ibas| {
                    self.basis_template[*ibas]
                        .iter()
                        .map(|bas| {
                            let mut rbas = bas.clone();
                            rbas.atom_of = *iatm as i32;
                            rbas
                        })
                        .collect::<Vec<CintBasis>>()
                })
            })
            .collect()
    }

    fn gen_bas(&self, iatm: Vec<AtomIndex>) -> Vec<CintBasis> {
        iatm.iter()
            .map(|i| self.atmbas_index[*i])
            .enumerate()
            .flat_map(|(iatm, ibas)| {
                let i = &iatm;
                self.basis_template[ibas]
                    .iter()
                    .map(|bas| {
                        let mut rbas = bas.clone();
                        rbas.atom_of = *i as i32;
                        rbas
                    })
                    .collect::<Vec<CintBasis>>()
            })
            .collect()
    }

    fn gen_bas_all(&self) -> Vec<CintBasis> {
        self.atmbas_index
            .iter()
            .enumerate()
            .flat_map(|(iatm, ibas)| {
                let i = &iatm;
                self.basis_template[*ibas]
                    .iter()
                    .map(|bas| {
                        let mut rbas = bas.clone();
                        rbas.atom_of = *i as i32;
                        rbas
                    })
                    .collect::<Vec<CintBasis>>()
            })
            .collect()
    }

    pub fn gen_intor_select<'a>(
        &'a self,
        which_bas: BTreeMap<AtomIndex, Vec<BasisIndex>>,
    ) -> Intor {
        Intor::new(&self.atm, self.gen_bas_all(), &self.env)
    }

    pub fn gen_intor<'a>(&'a self, iatm: Vec<AtomIndex>) -> Intor {
        Intor::new(&self.atm, self.gen_bas_all(), &self.env)
    }

    pub fn gen_intor_all<'a>(&'a self) -> Intor {
        Intor::new(&self.atm, self.gen_bas_all(), &self.env)
    }
}
