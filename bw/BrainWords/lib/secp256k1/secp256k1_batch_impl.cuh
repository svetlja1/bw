/**********************************************************************
 * Copyright (c) 2016 Llamasoft                                       *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _SECP256K1_BATCH_IMPL_H_
#define _SECP256K1_BATCH_IMPL_H_

#include <stddef.h>

typedef struct secp256k1_scratch_struct secp256k1_scratch;
typedef struct secp256k1_scratch_struct2 secp256k1_scratch2;
/* Scratch space for secp256k1_ec_pubkey_create_batch's temporary results. */
struct secp256k1_scratch_struct {
    /* Output from individual secp256k1_ecmult_gen. */
    secp256k1_gej gej[THREAD_STEPS];

    /* Input and output buffers for secp256k1_fe_inv_all_var. */
    secp256k1_fe  fe_in[THREAD_STEPS];
    secp256k1_fe  fe_out[THREAD_STEPS];
};
struct secp256k1_scratch_struct2 {
    /* Output from individual secp256k1_ecmult_gen. */
    secp256k1_gej gej[THREAD_STEPS_BRUTE];

    /* Input and output buffers for secp256k1_fe_inv_all_var. */
    secp256k1_fe  fe_in[THREAD_STEPS_BRUTE];
    secp256k1_fe  fe_out[THREAD_STEPS_BRUTE];
};

__device__ __inline__ size_t secp256k1_ec_pubkey_create_serialized_batch_myunsafe_brute(unsigned char* __restrict__  pubkeys, const unsigned char* __restrict__ privkeys, const int keyLen, const secp256k1_ge_storage* __restrict__ precPtr, const size_t precPitch) {
    secp256k1_scalar s_privkey;
    secp256k1_ge ge_pubkey;
    size_t dummy;
    int i, out_keys;
    secp256k1_scratch2 scr;
    out_keys = 0;
    int windowLimit = WINDOWS_SIZE_CONST[0];
    int windowMultLimit = ECMULT_WINDOW_SIZE_CONST[0];
#pragma unroll
    for (i = 0; i < keyLen; i++) {
        /* Convert private key to scalar form. */
        secp256k1_scalar_set_b32(&s_privkey, &(privkeys[32 * i]), NULL);
#ifdef ECMULT_BIG_TABLE
        secp256k1_ecmult_big(&(scr.gej[i]), &s_privkey, precPtr, precPitch, windowLimit, windowMultLimit);
#else
        secp256k1_ecmult_gen_fast(&(scr.gej[i]), &s_privkey, prec);
#endif // !ECM        

        if (scr.gej[i].infinity) { continue; }
        scr.fe_in[out_keys] = scr.gej[i].z;
        out_keys++;
    }
    if (out_keys > 0) {
        secp256k1_fe_inv_all_var(out_keys, scr.fe_out, scr.fe_in);
    }
    out_keys = 0;
    for (i = 0; i < keyLen; i++) {
        if (scr.gej[i].infinity) {
            continue;
        }
        secp256k1_ge_set_gej_zinv(&ge_pubkey, &(scr.gej[i]), &(scr.fe_out[out_keys]));
        secp256k1_eckey_pubkey_serialize(&ge_pubkey, &(pubkeys[65 * i]), &dummy, false);
        out_keys++;
    }
    return out_keys;
}

__device__ __inline__ size_t secp256k1_ec_pubkey_create_serialized_batch_myunsafe(unsigned char* __restrict__  pubkeys, const unsigned char* __restrict__ privkeys,  const int keyLen, const secp256k1_ge_storage* __restrict__ precPtr, const size_t precPitch) {
    secp256k1_scalar s_privkey;
    secp256k1_ge ge_pubkey;
    size_t dummy;
    int i, out_keys;
    secp256k1_scratch scr;
    /* Blank all of the output, regardless of what happens.                 */
    /* This marks all output keys as invalid until successfully created.    */
    //memset(pubkeys, 0, sizeof(*pubkeys) * pubkey_size * key_count);

    out_keys = 0;
    int windowLimit = WINDOWS_SIZE_CONST[0];
    int windowMultLimit = ECMULT_WINDOW_SIZE_CONST[0];
#pragma unroll
    for (i = 0; i < keyLen; i++) {
        /* Convert private key to scalar form. */
        secp256k1_scalar_set_b32(&s_privkey, &(privkeys[32 * i]), NULL);
#ifdef ECMULT_BIG_TABLE
        secp256k1_ecmult_big(&(scr.gej[i]), &s_privkey, precPtr, precPitch, windowLimit, windowMultLimit);
#else
        secp256k1_ecmult_gen_fast(&(scr.gej[i]), &s_privkey, prec);      
#endif // !ECM        

        /* If the result is the point at infinity, the pubkey is invalid. */
        if (scr.gej[i].infinity) { continue; }


        /* Save the Jacobian pubkey's Z coordinate for batch inversion. */
        scr.fe_in[out_keys] = scr.gej[i].z;
        out_keys++;
    }


    /* Assuming we have at least one non-infinite Jacobian pubkey. */
    if (out_keys > 0) {
        /* Invert all Jacobian public keys' Z values in one go. */
        secp256k1_fe_inv_all_var(out_keys, scr.fe_out, scr.fe_in);
    }


    /* Using the inverted Z values, convert each Jacobian public key to affine, */
    /*   then serialize the affine version to the pubkey buffer.                */
    out_keys = 0;

    for (i = 0; i < keyLen; i++) {
        /* Skip inverting infinite values. */
        /* The corresponding pubkey is already filled with \0 bytes from earlier. */
        if (scr.gej[i].infinity) {
            continue;
        }

        /* Otherwise, load the next inverted Z value and convert the pubkey to affine coordinates. */
        secp256k1_ge_set_gej_zinv(&ge_pubkey, &(scr.gej[i]), &(scr.fe_out[out_keys]));

        /* Serialize the public key into the requested format. */
        secp256k1_eckey_pubkey_serialize(&ge_pubkey, &(pubkeys[65 * i]), &dummy, false);
        out_keys++;
    }


    /* Returning the number of successfully converted private keys. */
    return out_keys;
}

#endif
