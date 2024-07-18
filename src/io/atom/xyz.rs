use std::collections::BTreeMap;

use super::super::convert::aatoau;
use super::super::PERIODIC_TABLE;
use crate::cint::{
    rawdata::{CintAtomGroup, CGTO},
    NAtom,
};

impl CintAtomGroup {
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
                    Ok(x) => coor[i] = aatoau(x),
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
