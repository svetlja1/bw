#define ECMULT_BIG_TABLE

#ifndef WIN64
//todo - check why with 64bit performance is lower
//#define ARCH 64
typedef unsigned __int128 uint128_t;
#define ARCH 32
#else
#define ARCH 32
#endif

#ifdef _DEBUG
#define ECMULT_GEN_PREC_BITS 4
//#define ECMULT_WINDOW_SIZE 4
#else
#define ECMULT_GEN_PREC_BITS 16
//#define ECMULT_WINDOW_SIZE 22
#endif 

#define ECMULT_GEN_PREC_B ECMULT_GEN_PREC_BITS

#ifdef ECMULT_BIG_TABLE
#define ECMULT_GEN_PREC_G (1 << ECMULT_GEN_PREC_B)
#define ECMULT_GEN_PREC_N (256 / ECMULT_GEN_PREC_B)
#else
#define ECMULT_GEN_PREC_G (1 << ECMULT_GEN_PREC_B)
#define ECMULT_GEN_PREC_N (256 / ECMULT_GEN_PREC_B)
#endif

__device__ __constant__ unsigned int ECMULT_WINDOW_SIZE_CONST[1];
__device__ __constant__ unsigned int WINDOWS_SIZE_CONST[1];

__device__ unsigned int ECMULT_WINDOW_SIZE;
__device__  size_t WINDOWS;
__device__  size_t WINDOW_SIZE;
#define ECMULT_TABLE_SIZE(w) ((256 / ECMULT_WINDOW_SIZE) * WINDOW_SIZE + (1 << (256 % ECMULT_WINDOW_SIZE)))


#define RIPEMD160_BLOCK_LENGTH 64
#define RIPEMD160_DIGEST_LENGTH 20
#define SECP256K1_FLAGS_TYPE_MASK ((1 << 8) - 1)
#define SECP256K1_FLAGS_TYPE_CONTEXT (1 << 0)
#define SECP256K1_FLAGS_TYPE_COMPRESSION (1 << 1)
#define SECP256K1_FLAGS_BIT_CONTEXT_VERIFY (1 << 8)
#define SECP256K1_FLAGS_BIT_CONTEXT_SIGN (1 << 9)
#define SECP256K1_FLAGS_BIT_CONTEXT_DECLASSIFY (1 << 10)
#define SECP256K1_FLAGS_BIT_COMPRESSION (1 << 8)
#define SECP256K1_EC_COMPRESSED (SECP256K1_FLAGS_TYPE_COMPRESSION | SECP256K1_FLAGS_BIT_COMPRESSION)
#define SECP256K1_EC_UNCOMPRESSED (SECP256K1_FLAGS_TYPE_COMPRESSION)
#define SECP256K1_TAG_PUBKEY_EVEN 0x02
#define SECP256K1_TAG_PUBKEY_ODD 0x03
#define SECP256K1_TAG_PUBKEY_UNCOMPRESSED 0x04
#define SECP256K1_TAG_PUBKEY_HYBRID_EVEN 0x06
#define SECP256K1_TAG_PUBKEY_HYBRID_ODD 0x07

#define SECP256K1_FE_CONST(d7, d6, d5, d4, d3, d2, d1, d0) {SECP256K1_FE_CONST_INNER((d7), (d6), (d5), (d4), (d3), (d2), (d1), (d0))}
#define SECP256K1_GE_CONST(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p) {SECP256K1_FE_CONST((a),(b),(c),(d),(e),(f),(g),(h)), SECP256K1_FE_CONST((i),(j),(k),(l),(m),(n),(o),(p)), 0}
#define SECP256K1_GE_CONST_INFINITY {SECP256K1_FE_CONST(0, 0, 0, 0, 0, 0, 0, 0), SECP256K1_FE_CONST(0, 0, 0, 0, 0, 0, 0, 0), 1}


#define VERIFY_CHECK(cond) do { (void)(cond); } while(0)

#if ARCH==32
typedef struct {
    uint32_t n[10];
} secp256k1_fe;

typedef struct {
    uint32_t n[8];
} secp256k1_fe_storage;

typedef struct {
    uint32_t d[8];
} secp256k1_scalar;
#define SECP256K1_SCALAR_CONST(d7, d6, d5, d4, d3, d2, d1, d0) {{(d0), (d1), (d2), (d3), (d4), (d5), (d6), (d7)}}
#define SECP256K1_N_0 ((uint32_t)0xD0364141UL)
#define SECP256K1_N_1 ((uint32_t)0xBFD25E8CUL)
#define SECP256K1_N_2 ((uint32_t)0xAF48A03BUL)
#define SECP256K1_N_3 ((uint32_t)0xBAAEDCE6UL)
#define SECP256K1_N_4 ((uint32_t)0xFFFFFFFEUL)
#define SECP256K1_N_5 ((uint32_t)0xFFFFFFFFUL)
#define SECP256K1_N_6 ((uint32_t)0xFFFFFFFFUL)
#define SECP256K1_N_7 ((uint32_t)0xFFFFFFFFUL)
#define SECP256K1_N_C_0 (~SECP256K1_N_0 + 1)
#define SECP256K1_N_C_1 (~SECP256K1_N_1)
#define SECP256K1_N_C_2 (~SECP256K1_N_2)
#define SECP256K1_N_C_3 (~SECP256K1_N_3)
#define SECP256K1_N_C_4 (1)
#define SECP256K1_FE_CONST_INNER(d7, d6, d5, d4, d3, d2, d1, d0) { (d0) & 0x3FFFFFFUL, (((uint32_t)d0) >> 26) | (((uint32_t)(d1) & 0xFFFFFUL) << 6), (((uint32_t)d1) >> 20) | (((uint32_t)(d2) & 0x3FFFUL) << 12), (((uint32_t)d2) >> 14) | (((uint32_t)(d3) & 0xFFUL) << 18), (((uint32_t)d3) >> 8) | (((uint32_t)(d4) & 0x3UL) << 24), (((uint32_t)d4) >> 2) & 0x3FFFFFFUL,(((uint32_t)d4) >> 28) | (((uint32_t)(d5) & 0x3FFFFFUL) << 4), (((uint32_t)d5) >> 22) | (((uint32_t)(d6) & 0xFFFFUL) << 10), (((uint32_t)d6) >> 16) | (((uint32_t)(d7) & 0x3FFUL) << 16), (((uint32_t)d7) >> 10) }
#define SECP256K1_FE_STORAGE_CONST(d7, d6, d5, d4, d3, d2, d1, d0) {{ (d0), (d1), (d2), (d3), (d4), (d5), (d6), (d7) }}
#else
typedef struct {
    /* X = sum(i=0..4, elem[i]*2^52) mod n */
    uint64_t n[5];
} secp256k1_fe;

typedef struct {
    uint64_t n[4];
} secp256k1_fe_storage;

typedef struct {
    uint64_t d[4];
} secp256k1_scalar;

#define SECP256K1_SCALAR_CONST(d7, d6, d5, d4, d3, d2, d1, d0) {{((uint64_t)(d1)) << 32 | (d0), ((uint64_t)(d3)) << 32 | (d2), ((uint64_t)(d5)) << 32 | (d4), ((uint64_t)(d7)) << 32 | (d6)}}
#define SECP256K1_N_0 ((uint64_t)0xBFD25E8CD0364141ULL)
#define SECP256K1_N_1 ((uint64_t)0xBAAEDCE6AF48A03BULL)
#define SECP256K1_N_2 ((uint64_t)0xFFFFFFFFFFFFFFFEULL)
#define SECP256K1_N_3 ((uint64_t)0xFFFFFFFFFFFFFFFFULL)
#define SECP256K1_N_C_0 (~SECP256K1_N_0 + 1)
#define SECP256K1_N_C_1 (~SECP256K1_N_1)
#define SECP256K1_N_C_2 (1)
#define SECP256K1_N_H_0 ((uint64_t)0xDFE92F46681B20A0ULL)
#define SECP256K1_N_H_1 ((uint64_t)0x5D576E7357A4501DULL)
#define SECP256K1_N_H_2 ((uint64_t)0xFFFFFFFFFFFFFFFFULL)
#define SECP256K1_N_H_3 ((uint64_t)0x7FFFFFFFFFFFFFFFULL)

#define SECP256K1_FE_CONST_INNER(d7, d6, d5, d4, d3, d2, d1, d0) { \
    (d0) | (((uint64_t)(d1) & 0xFFFFFUL) << 32), \
    ((uint64_t)(d1) >> 20) | (((uint64_t)(d2)) << 12) | (((uint64_t)(d3) & 0xFFUL) << 44), \
    ((uint64_t)(d3) >> 8) | (((uint64_t)(d4) & 0xFFFFFFFUL) << 24), \
    ((uint64_t)(d4) >> 28) | (((uint64_t)(d5)) << 4) | (((uint64_t)(d6) & 0xFFFFUL) << 36), \
    ((uint64_t)(d6) >> 16) | (((uint64_t)(d7)) << 16) \
}

#define SECP256K1_FE_STORAGE_CONST(d7, d6, d5, d4, d3, d2, d1, d0) {{ \
    (d0) | (((uint64_t)(d1)) << 32), \
    (d2) | (((uint64_t)(d3)) << 32), \
    (d4) | (((uint64_t)(d5)) << 32), \
    (d6) | (((uint64_t)(d7)) << 32) \
}}

#endif

typedef struct {
    secp256k1_fe x;
    secp256k1_fe y;
    int infinity;
} secp256k1_ge;

typedef struct {
    secp256k1_fe x;
    secp256k1_fe y;
    secp256k1_fe z;
    int infinity;
} secp256k1_gej;

typedef struct {
    secp256k1_fe_storage x;
    secp256k1_fe_storage y;
} secp256k1_ge_storage;

typedef struct {
    unsigned char data[64];
} secp256k1_pubkey;

typedef struct {
    int32_t v[9];
} secp256k1_modinv32_signed30;

typedef struct {
    // The modulus in signed30 notation, must be odd and in [3, 2^256]. 
    secp256k1_modinv32_signed30 modulus;

    // modulus^{-1} mod 2^30 
    uint32_t modulus_inv30;
} secp256k1_modinv32_modinfo;

__device__ __constant__ secp256k1_modinv32_modinfo secp256k1_const_modinfo_fe = {
    {{-0x3D1, -4, 0, 0, 0, 0, 0, 0, 65536}},
    0x2DDACACFL
};

#define SECP256K1_GE_STORAGE_CONST(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p) {SECP256K1_FE_STORAGE_CONST((a),(b),(c),(d),(e),(f),(g),(h)), SECP256K1_FE_STORAGE_CONST((i),(j),(k),(l),(m),(n),(o),(p))}
#define SC SECP256K1_GE_STORAGE_CONST

#define SECP256K1_G SECP256K1_GE_CONST(\
    0x79BE667EUL, 0xF9DCBBACUL, 0x55A06295UL, 0xCE870B07UL,\
    0x029BFCDBUL, 0x2DCE28D9UL, 0x59F2815BUL, 0x16F81798UL,\
    0x483ADA77UL, 0x26A3C465UL, 0x5DA4FBFCUL, 0x0E1108A8UL,\
    0xFD17B448UL, 0xA6855419UL, 0x9C47D08FUL, 0xFB10D4B8UL\
)

__device__ __constant__ secp256k1_ge secp256k1_ge_const_g = SECP256K1_GE_CONST(
    0x79BE667EUL, 0xF9DCBBACUL, 0x55A06295UL, 0xCE870B07UL,
    0x029BFCDBUL, 0x2DCE28D9UL, 0x59F2815BUL, 0x16F81798UL,
    0x483ADA77UL, 0x26A3C465UL, 0x5DA4FBFCUL, 0x0E1108A8UL,
    0xFD17B448UL, 0xA6855419UL, 0x9C47D08FUL, 0xFB10D4B8UL
);


__device__ __constant__ uint32_t M = 0x3FFFFFFUL, R0 = 0x3D10UL, R1 = 0x400UL;

__device__ __constant__ uint32_t M30u = UINT32_MAX >> 2;
__device__ __constant__ uint32_t M26u = UINT32_MAX >> 6;
__device__ __constant__ int32_t M30 = (int32_t)(UINT32_MAX >> 2);


__device__ __constant__ secp256k1_scalar SCALAR_ONE = SECP256K1_SCALAR_CONST(0, 0, 0, 0, 0, 0, 0, 1);
__device__ __constant__ secp256k1_fe FE_ONE = SECP256K1_FE_CONST(0, 0, 0, 0, 0, 0, 0, 1);

__device__ __inline__ int secp256k1_ctz32_var(const uint32_t x) {
    return __clz(__brev(x));
    //return secp256k1_ctz32_var_debruijn(x);
}

/* inv256[i] = -(2*i+1)^-1 (mod 256) */
__device__ __constant__ uint8_t inv256[128] = {
    0xFF, 0x55, 0x33, 0x49, 0xC7, 0x5D, 0x3B, 0x11, 0x0F, 0xE5, 0xC3, 0x59,
    0xD7, 0xED, 0xCB, 0x21, 0x1F, 0x75, 0x53, 0x69, 0xE7, 0x7D, 0x5B, 0x31,
    0x2F, 0x05, 0xE3, 0x79, 0xF7, 0x0D, 0xEB, 0x41, 0x3F, 0x95, 0x73, 0x89,
    0x07, 0x9D, 0x7B, 0x51, 0x4F, 0x25, 0x03, 0x99, 0x17, 0x2D, 0x0B, 0x61,
    0x5F, 0xB5, 0x93, 0xA9, 0x27, 0xBD, 0x9B, 0x71, 0x6F, 0x45, 0x23, 0xB9,
    0x37, 0x4D, 0x2B, 0x81, 0x7F, 0xD5, 0xB3, 0xC9, 0x47, 0xDD, 0xBB, 0x91,
    0x8F, 0x65, 0x43, 0xD9, 0x57, 0x6D, 0x4B, 0xA1, 0x9F, 0xF5, 0xD3, 0xE9,
    0x67, 0xFD, 0xDB, 0xB1, 0xAF, 0x85, 0x63, 0xF9, 0x77, 0x8D, 0x6B, 0xC1,
    0xBF, 0x15, 0xF3, 0x09, 0x87, 0x1D, 0xFB, 0xD1, 0xCF, 0xA5, 0x83, 0x19,
    0x97, 0xAD, 0x8B, 0xE1, 0xDF, 0x35, 0x13, 0x29, 0xA7, 0x3D, 0x1B, 0xF1,
    0xEF, 0xC5, 0xA3, 0x39, 0xB7, 0xCD, 0xAB, 0x01
};