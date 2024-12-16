pub mod params;
use crate::params::{FIRST_FOUR_ROUND_RC, LAST_FOUR_ROUNDS_RC, MAT_DIAG16_M_1, PARTIAL_ROUNDS_RC};
use sha2::Digest;

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

#[cfg(test)]
mod tests {
    use crate::poseidon2_permute;

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
}
