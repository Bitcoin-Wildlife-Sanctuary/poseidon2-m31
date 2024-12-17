pub mod params;

use crate::params::{FIRST_FOUR_ROUND_RC, LAST_FOUR_ROUNDS_RC, MAT_DIAG16_M_1, PARTIAL_ROUNDS_RC};
use sha2::Digest;
use std::cmp::PartialEq;

pub fn modp(a: u32) -> u32 {
    a % ((1 << 31) - 1)
}

pub fn addmod(a: u32, b: u32) -> u32 {
    (a + b) % ((1 << 31) - 1)
}

pub fn apply_4x4_mds_matrix(x0: &mut u32, x1: &mut u32, x2: &mut u32, x3: &mut u32) {
    let t0 = addmod(*x0, *x1);
    let t1 = addmod(*x2, *x3);

    let x1_times_2 = addmod(*x1, *x1);
    let x3_times_2 = addmod(*x3, *x3);

    let t2 = addmod(x1_times_2, t1);
    let t3 = addmod(x3_times_2, t0);

    let t1_times_2 = addmod(t1, t1);
    let t1_times_4 = addmod(t1_times_2, t1_times_2);

    let t0_times_2 = addmod(t0, t0);
    let t0_times_4 = addmod(t0_times_2, t0_times_2);

    let t4 = addmod(t1_times_4, t3);
    let t5 = addmod(t0_times_4, t2);

    let t6 = addmod(t3, t5);
    let t7 = addmod(t2, t4);

    *x0 = t6;
    *x1 = t5;
    *x2 = t7;
    *x3 = t4;
}

pub fn apply_16x16_mds_matrix(state: &mut [u32; 16]) {
    let mut t0 = state[0];
    let mut t1 = state[1];
    let mut t2 = state[2];
    let mut t3 = state[3];
    let mut t4 = state[4];
    let mut t5 = state[5];
    let mut t6 = state[6];
    let mut t7 = state[7];
    let mut t8 = state[8];
    let mut t9 = state[9];
    let mut t10 = state[10];
    let mut t11 = state[11];
    let mut t12 = state[12];
    let mut t13 = state[13];
    let mut t14 = state[14];
    let mut t15 = state[15];

    apply_4x4_mds_matrix(&mut t0, &mut t1, &mut t2, &mut t3);
    apply_4x4_mds_matrix(&mut t4, &mut t5, &mut t6, &mut t7);
    apply_4x4_mds_matrix(&mut t8, &mut t9, &mut t10, &mut t11);
    apply_4x4_mds_matrix(&mut t12, &mut t13, &mut t14, &mut t15);

    let t = [
        t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15,
    ];

    for i in 0..4 {
        state[i] = addmod(t[i], t[i]);
        state[i] = addmod(state[i], t[i + 4]);
        state[i] = addmod(state[i], t[i + 8]);
        state[i] = addmod(state[i], t[i + 12]);
    }

    for i in 4..8 {
        state[i] = addmod(t[i], t[i]);
        state[i] = addmod(state[i], t[i - 4]);
        state[i] = addmod(state[i], t[i + 4]);
        state[i] = addmod(state[i], t[i + 8]);
    }

    for i in 8..12 {
        state[i] = addmod(t[i], t[i]);
        state[i] = addmod(state[i], t[i - 8]);
        state[i] = addmod(state[i], t[i - 4]);
        state[i] = addmod(state[i], t[i + 4]);
    }

    for i in 12..16 {
        state[i] = addmod(t[i], t[i]);
        state[i] = addmod(state[i], t[i - 12]);
        state[i] = addmod(state[i], t[i - 8]);
        state[i] = addmod(state[i], t[i - 4]);
    }
}

pub fn pow5mod(a: u32) -> u32 {
    let a = a as u64;

    let pow2mod = a * a % ((1 << 31) - 1);
    let pow4mod = pow2mod * pow2mod % ((1 << 31) - 1);
    (pow4mod * a % ((1 << 31) - 1)) as u32
}

pub fn mulmod(a: u32, b: u32) -> u32 {
    (((a as u64) * (b as u64)) % ((1 << 31) - 1)) as u32
}

pub fn compute_iv_values(domain_separator: &[u8]) -> [u32; 8] {
    let mut sha256 = sha2::Sha256::new();
    Digest::update(&mut sha256, domain_separator);

    let bytes = sha256.finalize().to_vec();

    [
        modp(u32::from_be_bytes(bytes[0..4].try_into().unwrap())),
        modp(u32::from_be_bytes(bytes[4..8].try_into().unwrap())),
        modp(u32::from_be_bytes(bytes[8..12].try_into().unwrap())),
        modp(u32::from_be_bytes(bytes[12..16].try_into().unwrap())),
        modp(u32::from_be_bytes(bytes[16..20].try_into().unwrap())),
        modp(u32::from_be_bytes(bytes[20..24].try_into().unwrap())),
        modp(u32::from_be_bytes(bytes[24..28].try_into().unwrap())),
        modp(u32::from_be_bytes(bytes[28..32].try_into().unwrap())),
    ]
}

pub fn poseidon2_permute(state: &mut [u32; 16]) {
    apply_16x16_mds_matrix(state);

    for r in 0..4 {
        for i in 0..16 {
            state[i] = addmod(state[i], FIRST_FOUR_ROUND_RC[r][i]);
        }

        for i in 0..16 {
            state[i] = pow5mod(state[i]);
        }

        apply_16x16_mds_matrix(state);
    }

    for r in 0..14 {
        state[0] = addmod(state[0], PARTIAL_ROUNDS_RC[r]);
        state[0] = pow5mod(state[0]);

        let mut sum = 0;
        for i in 0..16 {
            sum = addmod(sum, state[i]);
        }

        for i in 0..16 {
            state[i] = addmod(sum, mulmod(state[i], MAT_DIAG16_M_1[i]));
        }
    }

    for r in 0..4 {
        for i in 0..16 {
            state[i] = addmod(state[i], LAST_FOUR_ROUNDS_RC[r][i]);
        }

        for i in 0..16 {
            state[i] = pow5mod(state[i]);
        }

        apply_16x16_mds_matrix(state);
    }
}

pub struct Poseidon31Hasher {
    pub buffer: Vec<u32>,
}

impl Poseidon31Hasher {
    pub fn new() -> Self {
        Self { buffer: Vec::new() }
    }

    pub fn update(&mut self, v: impl AsRef<[u32]>) {
        self.buffer.extend_from_slice(v.as_ref())
    }

    pub fn finalize(&self, description: impl ToString) -> [u32; 8] {
        let mut input = self.buffer.clone();
        let l = input.len();

        let iv = compute_iv_values(
            format!(
                "Poseidon2 M31 hashing {} M31 elements for {}",
                l,
                description.to_string()
            )
            .as_bytes(),
        );

        let mut state = [
            0, 0, 0, 0, 0, 0, 0, 0, iv[0], iv[1], iv[2], iv[3], iv[4], iv[5], iv[6], iv[7],
        ];
        input.resize(l.div_ceil(8) * 8, 0);

        for chunk in input.chunks(8) {
            state[0..8].copy_from_slice(chunk);
            poseidon2_permute(&mut state);
        }

        [
            state[0], state[1], state[2], state[3], state[4], state[5], state[6], state[7],
        ]
    }

    pub fn finalize_reset(&mut self, description: impl ToString) -> [u32; 8] {
        let res = self.finalize(description);
        self.buffer.clear();
        res
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub enum Poseidon31Mode {
    ABSORB,
    SQUEEZE,
}

#[derive(Clone)]
pub struct Poseidon31Sponge {
    pub state: [u32; 16],
    pub buffer: Vec<u32>,
    pub mode: Poseidon31Mode,
    pub squeeze_index: usize,
}

impl Default for Poseidon31Sponge {
    fn default() -> Self {
        Self::new("Poseidon2 M31 sponge")
    }
}

impl Poseidon31Sponge {
    pub fn new(description: impl ToString) -> Self {
        let iv = compute_iv_values(
            format!("Poseidon2 M31 sponge for {}", description.to_string()).as_bytes(),
        );
        Self {
            state: [
                0, 0, 0, 0, 0, 0, 0, 0, iv[0], iv[1], iv[2], iv[3], iv[4], iv[5], iv[6], iv[7],
            ],
            buffer: Vec::new(),
            mode: Poseidon31Mode::ABSORB,
            squeeze_index: 0,
        }
    }

    pub fn absorb(&mut self, data: &[u32]) {
        if data.is_empty() {
            return;
        }

        self.mode = Poseidon31Mode::ABSORB;
        self.squeeze_index = 0;

        let bl = self.buffer.len();
        let l = bl + data.len();

        if l >= 8 {
            self.state[0..bl].copy_from_slice(self.buffer.as_slice());
            self.state[bl..8].copy_from_slice(&data[0..8 - bl]);
            self.buffer.clear();
            poseidon2_permute(&mut self.state);

            self.absorb(&data[8 - bl..]);
        } else {
            self.buffer.extend_from_slice(data);
        }
    }

    pub fn squeeze(&mut self, size: usize) -> Vec<u32> {
        assert!(size > 0);

        if self.mode == Poseidon31Mode::ABSORB && !self.buffer.is_empty() {
            for i in 0..8 {
                self.state[i] = 0;
            }
            self.state[0..self.buffer.len()].copy_from_slice(self.buffer.as_slice());
            self.buffer.clear();
            poseidon2_permute(&mut self.state);

            self.squeeze_index = 0;
        }

        self.mode = Poseidon31Mode::SQUEEZE;

        let mut res = vec![];

        while res.len() < size {
            let need = size - res.len();
            let ready = 8 - self.squeeze_index;

            if ready > need {
                res.extend_from_slice(&self.state[self.squeeze_index..self.squeeze_index + need]);
                self.squeeze_index += need;
            } else {
                res.extend_from_slice(&self.state[self.squeeze_index..8]);

                for i in 0..8 {
                    self.state[i] = 0;
                }
                poseidon2_permute(&mut self.state);
                self.squeeze_index = 0;
            }
        }

        res
    }
}

pub struct Poseidon31CRH;

impl Poseidon31CRH {
    pub fn hash_fixed_length(data: &[u32]) -> [u32; 8] {
        let mut cur_layer = vec![];
        for chunk in data.chunks(16) {
            let mut cur = [0u32; 16];
            cur[0..chunk.len()].copy_from_slice(chunk);

            poseidon2_permute(&mut cur);

            for i in 0..std::cmp::min(8, chunk.len()) {
                cur[i] = addmod(cur[i], chunk[i]);
            }

            cur_layer.push([
                cur[0], cur[1], cur[2], cur[3], cur[4], cur[5], cur[6], cur[7],
            ]);
        }

        while cur_layer.len() > 1 {
            let mut new_layer = vec![];
            for chunk in cur_layer.chunks(2) {
                if chunk.len() == 1 {
                    new_layer.push(chunk[0]);
                } else {
                    let mut cur = [0u32; 16];
                    for i in 0..8 {
                        cur[i] = chunk[0][i];
                        cur[i + 8] = chunk[1][i];
                    }

                    poseidon2_permute(&mut cur);

                    for i in 0..8 {
                        cur[i] = addmod(cur[i], chunk[0][i]);
                    }

                    new_layer.push([
                        cur[0], cur[1], cur[2], cur[3], cur[4], cur[5], cur[6], cur[7],
                    ]);
                }
            }

            cur_layer = new_layer;
        }

        cur_layer[0]
    }
}

#[cfg(test)]
mod tests {
    use crate::{modp, poseidon2_permute, Poseidon31CRH, Poseidon31Sponge};
    use rand::Rng;
    use rand_chacha::rand_core::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_poseidon2_permute() {
        let mut state = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
        poseidon2_permute(&mut state);

        let expected = [
            1348310665, 996460804, 2044919169, 1269301599, 615961333, 595876573, 1377780500,
            1776267289, 715842585, 1823756332, 1870636634, 1979645732, 311256455, 1364752356,
            58674647, 323699327,
        ];
        assert_eq!(state, expected);
    }

    #[test]
    fn test_absorb_and_squeeze() {
        let mut sponge = Poseidon31Sponge::new("test");
        let c1 = sponge.state[8..16].to_vec();

        sponge.absorb(&[1, 2, 3, 4, 5, 6, 7, 8, 9]);

        let c2 = {
            let mut state = [
                1, 2, 3, 4, 5, 6, 7, 8, c1[0], c1[1], c1[2], c1[3], c1[4], c1[5], c1[6], c1[7],
            ];
            poseidon2_permute(&mut state);
            state[8..16].to_vec()
        };

        assert_eq!(sponge.state[8..16], c2);

        let a = sponge.squeeze(1);

        let c3_full = {
            let mut state = [
                9, 0, 0, 0, 0, 0, 0, 0, c2[0], c2[1], c2[2], c2[3], c2[4], c2[5], c2[6], c2[7],
            ];
            poseidon2_permute(&mut state);
            state
        };

        assert_eq!(a, vec![c3_full[0]]);

        let a = sponge.squeeze(2);
        assert_eq!(a, vec![c3_full[1], c3_full[2]]);

        sponge.absorb(&[10, 11, 12]);
        sponge.absorb(&[13]);

        let a = sponge.squeeze(3);

        let c3 = c3_full[8..16].to_vec();

        let c4_full = {
            let mut state = [
                10, 11, 12, 13, 0, 0, 0, 0, c3[0], c3[1], c3[2], c3[3], c3[4], c3[5], c3[6], c3[7],
            ];
            poseidon2_permute(&mut state);
            state
        };
        assert_eq!(a, vec![c4_full[0], c4_full[1], c4_full[2]]);

        let _ = sponge.squeeze(8);
        let a = sponge.squeeze(1);

        let c4 = c4_full[8..16].to_vec();

        let c5_full = {
            let mut state = [
                0, 0, 0, 0, 0, 0, 0, 0, c4[0], c4[1], c4[2], c4[3], c4[4], c4[5], c4[6], c4[7],
            ];
            poseidon2_permute(&mut state);
            state
        };
        assert_eq!(a, vec![c5_full[3]]);

        sponge.absorb(&[13]);
        let a = sponge.squeeze(2);

        let c5 = c5_full[8..16].to_vec();
        let c6_full = {
            let mut state = [
                13, 0, 0, 0, 0, 0, 0, 0, c5[0], c5[1], c5[2], c5[3], c5[4], c5[5], c5[6], c5[7],
            ];
            poseidon2_permute(&mut state);
            state
        };

        assert_eq!(a, vec![c6_full[0], c6_full[1]]);
    }

    #[test]
    fn test_commit() {
        let mut prng = ChaCha20Rng::seed_from_u64(0);

        let mut a = [0u32; 5];
        let mut b = [0u32; 17];

        for i in 0..5 {
            a[i] = modp(prng.gen::<u32>());
        }

        for i in 0..17 {
            b[i] = modp(prng.gen::<u32>());
        }

        let _ = Poseidon31CRH::hash_fixed_length(&a);
        let _ = Poseidon31CRH::hash_fixed_length(&b);
    }
}
