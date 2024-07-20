#[allow(unused)]
pub mod cint;
#[allow(unused)]
pub mod io;

#[cfg(test)]
mod tests {

    #[test]
    fn test_rawdata() {
        use crate::cint::{cdata::CintDate, libcint::int1e_ovlp_cart};
        use std::collections::BTreeMap;

        //
        let work_path = std::env::current_dir().unwrap();

        let basis_path = format!(
            "{}/basis_set_exchange/basis_set_exchange/data/sto/STO-3G.1.json",
            work_path.to_str().unwrap()
        );

        let xyz_str = "3

        H -0.5  0.0  0.0
        O  0.0  0.0  0.0
        H  0.5  0.0  0.0
        ";

        let cint_data = CintDate::fron_xyz(xyz_str, &basis_path);

        //
        let intor_all = cint_data.gen_intor_all();
        let out = unsafe { intor_all.int_cart([0, 2], Some(int1e_ovlp_cart)) };
        println!("{:?}", out);

        //
        let intor_atm = cint_data.gen_intor(vec![0, 2]);
        let out = unsafe { intor_atm.int_cart([0, 1], Some(int1e_ovlp_cart)) };
        println!("{:?}", out);

        //
        let mut some_bas = BTreeMap::new();
        some_bas.insert(0, vec![0]);
        some_bas.insert(2, vec![0]);
        let intor_bas = cint_data.gen_intor_select(some_bas);
        let out = unsafe { intor_bas.int_cart([0, 1], Some(int1e_ovlp_cart)) };
        println!("{:?}", out);
    }
}
