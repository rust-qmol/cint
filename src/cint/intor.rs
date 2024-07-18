use std::ops::Index;

use super::{
    cdata::{CintAtom, CintBasis, CintEnv},
    libccint::{CINTIntegralFunction, CINTcgtos_cart, CINTcgtos_spheric, CINTcgtos_spinor},
    NAtom,
};

#[derive(Debug, Clone)]
pub struct IntorResult<const N: NAtom> {
    dims: Vec<i32>,
    out: Vec<f64>,
}

impl<const N: NAtom> IntorResult<N> {
    pub fn dims(&self) -> &Vec<i32> {
        &self.dims
    }

    pub fn out(&self) -> &Vec<f64> {
        &self.out
    }
}

impl<const N: NAtom> Index<Vec<usize>> for IntorResult<N> {
    type Output = f64;

    fn index(&self, index: Vec<usize>) -> &Self::Output {
        assert_eq!(index.len(), N);
        let out_i = index.last().unwrap()
            + self.dims[1..N]
                .iter()
                .rev()
                .enumerate()
                .fold(0, |i, (p, x)| (p + index[i]) * (*x as usize));
        &self.out[out_i]
    }
}

#[derive(Debug, Clone)]
pub struct Intor<'a> {
    dims_cart: Vec<i32>,
    dims_spheric: Vec<i32>,
    dims_spinor: Vec<i32>,
    natm: i32,
    atm: &'a Vec<CintAtom>,
    nbas: i32,
    bas: Vec<CintBasis>,
    env: &'a CintEnv,
}

impl<'a> Intor<'a> {
    pub fn new(atm: &'a Vec<CintAtom>, bas: Vec<CintBasis>, env: &'a CintEnv) -> Self {
        let dims_cart = bas
            .iter()
            .enumerate()
            .map(|(i, b)| unsafe { CINTcgtos_cart(i as i32, bas.as_ptr() as *const i32) })
            .collect();
        let dims_spheric = bas
            .iter()
            .enumerate()
            .map(|(i, b)| unsafe { CINTcgtos_spheric(i as i32, bas.as_ptr() as *const i32) })
            .collect();
        let dims_spinor = bas
            .iter()
            .enumerate()
            .map(|(i, b)| unsafe { CINTcgtos_spinor(i as i32, bas.as_ptr() as *const i32) })
            .collect();
        Self {
            dims_cart,
            dims_spheric,
            dims_spinor,
            natm: atm.len() as i32,
            atm,
            nbas: bas.len() as i32,
            bas,
            env,
        }
    }

    pub fn natm(&self) -> i32 {
        self.natm
    }

    pub fn nbas(&self) -> i32 {
        self.nbas
    }

    pub unsafe fn int_cart<const N: NAtom>(
        &self,
        shls: [i32; N],
        int_func: CINTIntegralFunction,
    ) -> IntorResult<N> {
        let dims: Vec<i32> = shls.iter().map(|l| self.dims_cart[*l as usize]).collect();
        let mut out = vec![0.0; dims.iter().fold(1, |p, x| p * x) as usize];
        match int_func {
            Some(func) => {
                unsafe {
                    func(
                        out.as_mut_ptr(),
                        dims.as_ptr(),
                        shls.as_ptr(),
                        self.atm.as_ptr() as *const i32,
                        self.natm,
                        self.bas.as_ptr() as *const i32,
                        self.nbas,
                        self.env.as_ptr(),
                        std::ptr::null(),
                        std::ptr::null_mut(),
                    );
                }
                ()
            }
            None => panic!(""),
        }
        IntorResult { dims, out }
    }
}
