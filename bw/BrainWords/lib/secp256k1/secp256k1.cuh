/*
__device__ __inline__ void memczero(void* s, size_t len, int flag) {
    unsigned char* p = (unsigned char*)s;
    volatile int vflag = flag;
    unsigned char mask = -(unsigned char)vflag;
    while (len) {
        *p &= ~mask;
        p++;
        len--;
    }
}
*/

__device__ __inline__ uint64_t secp256k1_scalar_shr_any(secp256k1_scalar* __restrict__ s, unsigned int n) {
    unsigned int cur_shift = 0, offset = 0;
    uint64_t rtn = 0;

    //VERIFY_CHECK(s != NULL);
    //VERIFY_CHECK(n > 0);
    //VERIFY_CHECK(n <= 64);

    
    while (n > 0) {
        // Shift up to 15 bits at a time, or N bits, whichever is smaller.  
        // secp256k1_scalar_shr_int() is hard limited to (0 < n < 16).      
        cur_shift = (n > 15 ? 15 : n);

        rtn |= ((uint64_t)secp256k1_scalar_shr_int(s, cur_shift) << (uint64_t)offset);

        offset += cur_shift;
        n -= cur_shift;
    }

    return rtn;
}


__device__ __inline__ int64_t secp256k1_scalar_sdigit_single(secp256k1_scalar* __restrict__  s, const unsigned int w) {
    int64_t sdigit = 0;

    /* Represents a 1 bit in the next window's least significant bit.       */
    /* VERIFY_CHECK verifies that (1 << w) won't touch int64_t's sign bit.  */
    int64_t overflow_bit = (int64_t)(1 << w);

    /* Represents the maximum positive value in a w-bit precomp table.  */
    /* Values greater than this are converted to negative values and    */
    /*   will "reverse borrow" a bit from the next window.              */
    int64_t precomp_max = (int64_t)(1 << (w - 1));

    //VERIFY_CHECK(s != NULL);
    //VERIFY_CHECK(w >= 1);
    //VERIFY_CHECK(w <= 62);

    sdigit = (int64_t)secp256k1_scalar_shr_any(s, w);

    if (sdigit <= precomp_max) {
        /* A w-bit precomp table has this digit as a positive value, return as-is.  */
        return sdigit;

    }
    else {
        secp256k1_scalar one;
        secp256k1_scalar_set_int(&one, 1);

        /* Convert this digit to a negative value, but balance s by adding it's value.  */
        /* Subtracting our sdigit value carries over into a 1 bit of the next digit.    */
        /* Since s has been shifted down w bits, s += 1 does the same thing.            */
        sdigit -= overflow_bit;

        secp256k1_scalar_add(s, s, &one);

        return sdigit;
    }
}

__device__ __inline__ void secp256k1_ecmult_gen_fast(secp256k1_gej* r, secp256k1_scalar* gn, const secp256k1_ge_storage _prec[ECMULT_GEN_PREC_N][ECMULT_GEN_PREC_G]) {
    secp256k1_ge add;
    int bits;
    secp256k1_gej_set_infinity(r);

    add.infinity = 0;
    for (int j = 0; j < ECMULT_GEN_PREC_N; j++) {
        bits = secp256k1_scalar_get_bits(gn, j * ECMULT_GEN_PREC_B, ECMULT_GEN_PREC_B);
        secp256k1_ge_from_storage(&add, &_prec[j][bits]);
        secp256k1_gej_add_ge(r, r, &add);
    }
}
/*
__device__ __inline__ void secp256k1_ecmult_gen(secp256k1_gej *r,  secp256k1_scalar *gn) {
    secp256k1_ge add;
    secp256k1_ge_storage adds;
    int bits;
    int i, j;

    memset(&adds, 0, sizeof(adds));
    secp256k1_gej_set_infinity(r);

    add.infinity = 0;
    for (j = 0; j < ECMULT_GEN_PREC_N; j++) {
        bits = secp256k1_scalar_get_bits(gn, j * ECMULT_GEN_PREC_B, ECMULT_GEN_PREC_B);
        for (i = 0; i < ECMULT_GEN_PREC_G; i++) {
            uint32_t mask0, mask1;
            mask0 = (i == bits) + ~((uint32_t)0);
            mask1 = ~mask0;
            
            adds.x.n[0] = (adds.x.n[0] & mask0) | (prec[j][i].x.n[0] & mask1);
            adds.x.n[1] = (adds.x.n[1] & mask0) | (prec[j][i].x.n[1] & mask1);
            adds.x.n[2] = (adds.x.n[2] & mask0) | (prec[j][i].x.n[2] & mask1);
            adds.x.n[3] = (adds.x.n[3] & mask0) | (prec[j][i].x.n[3] & mask1);
            adds.x.n[4] = (adds.x.n[4] & mask0) | (prec[j][i].x.n[4] & mask1);
            adds.x.n[5] = (adds.x.n[5] & mask0) | (prec[j][i].x.n[5] & mask1);
            adds.x.n[6] = (adds.x.n[6] & mask0) | (prec[j][i].x.n[6] & mask1);
            adds.x.n[7] = (adds.x.n[7] & mask0) | (prec[j][i].x.n[7] & mask1);

            adds.y.n[0] = (adds.y.n[0] & mask0) | (prec[j][i].y.n[0] & mask1);
            adds.y.n[1] = (adds.y.n[1] & mask0) | (prec[j][i].y.n[1] & mask1);
            adds.y.n[2] = (adds.y.n[2] & mask0) | (prec[j][i].y.n[2] & mask1);
            adds.y.n[3] = (adds.y.n[3] & mask0) | (prec[j][i].y.n[3] & mask1);
            adds.y.n[4] = (adds.y.n[4] & mask0) | (prec[j][i].y.n[4] & mask1);
            adds.y.n[5] = (adds.y.n[5] & mask0) | (prec[j][i].y.n[5] & mask1);
            adds.y.n[6] = (adds.y.n[6] & mask0) | (prec[j][i].y.n[6] & mask1);
            adds.y.n[7] = (adds.y.n[7] & mask0) | (prec[j][i].y.n[7] & mask1);
        }
        secp256k1_ge_from_storage(&add, &adds);
        
        secp256k1_ge_from_storage(&add, &prec[j][bits]);
        secp256k1_gej_add_ge(r, r, &add);
    }
    bits = 0;
    secp256k1_ge_clear(&add);
}
*/
/*
__device__ __inline__ void secp256k1_pubkey_save(secp256k1_pubkey* pubkey, secp256k1_ge* ge) {
  secp256k1_fe_normalize_var(&ge->x);
  secp256k1_fe_normalize_var(&ge->y);
  secp256k1_fe_get_b32(pubkey->data, &ge->x);
  secp256k1_fe_get_b32(pubkey->data + 32, &ge->y);
}
*/
/*
__device__ __inline__ int secp256k1_ec_pubkey_xyz(secp256k1_gej* pj, const unsigned char* seckey, secp256k1_ge_storage _prec[ECMULT_GEN_PREC_N][ECMULT_GEN_PREC_G]) {
    //secp256k1_gej pj;
    secp256k1_scalar sec;
    secp256k1_scalar secp256k1_scalar_one = SECP256K1_SCALAR_CONST(0, 0, 0, 0, 0, 0, 0, 1);
    int ret = 0;
    ret = secp256k1_scalar_set_b32_seckey(&sec, seckey);
    secp256k1_scalar_cmov(&sec, &secp256k1_scalar_one, !ret);
    secp256k1_ecmult_gen_fast(pj, &sec, _prec);
    //secp256k1_ge_set_gej(p, &pj);

    //secp256k1_pubkey_save(pubkey, p);
    //memczero(pubkey, sizeof(*pubkey), !ret);

    secp256k1_scalar_clear(&sec);
    return ret;
}*/


/** Multiply with the generator: R = a*G.
 *
 *  Args:   bmul:   pointer to an ecmult_big_context (cannot be NULL)
 *  Out:    r:      set to a*G where G is the generator (cannot be NULL)
 *  In:     a:      the scalar to multiply the generator by (cannot be NULL)
 */
__device__ __inline__ void secp256k1_ecmult_big(secp256k1_gej* __restrict__ r, const secp256k1_scalar* __restrict__ a, const secp256k1_ge_storage* __restrict__ precPtr, const size_t precPitch, const int windowLimit, const unsigned int windowEcmultLimit) {
    unsigned int window = 0;
    int64_t sdigit = 0;
    secp256k1_ge window_value;
    /* Copy of the input scalar which secp256k1_scalar_sdigit_single will destroy. */
    secp256k1_scalar privkey = *a;

    //VERIFY_CHECK(r != NULL);
    //VERIFY_CHECK(a != NULL);

    /* Until we hit a non-zero window, the value of r is undefined. */
    secp256k1_gej_set_infinity(r);

    /* If the privkey is zero, bail. */
    if (secp256k1_scalar_is_zero(&privkey)) { return; }
    
    /* Incrementally convert the privkey into signed digit form, one window at a time. */
    while (window < windowLimit && !secp256k1_scalar_is_zero(&privkey)) {
        sdigit = secp256k1_scalar_sdigit_single(&privkey, windowEcmultLimit);

        /* Zero windows have no representation in our precomputed table. */
      
        if (sdigit != 0) {

            secp256k1_ge_storage* ROW_PREC = (secp256k1_ge_storage*)((char*)precPtr + window * precPitch) + (llabs(sdigit) - 1);
            secp256k1_ge_from_storage(&window_value, ROW_PREC);

            if (sdigit < 0) {
                /* Use the positive precomp index and negate the result. */
                
                secp256k1_ge_neg(&window_value, &window_value);
            }
            else {
                /* Use the precomp index and result as-is.  */
                //secp256k1_ge_storage* ROW_PREC = (secp256k1_ge_storage*)((char*)precPtr + window * precPitch) + (+(sdigit)-1);
                //secp256k1_ge_from_storage(&window_value, ROW_PREC);
            }

            /* The first addition is automatically replaced by a load when r = inf. */
            secp256k1_gej_add_ge_var(r, r, &window_value, NULL);
        }

        window++;
    }

    /* If privkey isn't zero, something broke.  */
    //VERIFY_CHECK(secp256k1_scalar_is_zero(&privkey));
}



/*
__device__ __inline__ int secp256k1_eckey_privkey_tweak_add(secp256k1_scalar *key, const secp256k1_scalar *tweak) {
    secp256k1_scalar_add(key, key, tweak);
    return !secp256k1_scalar_is_zero(key);
}
*/
/*
__device__ __inline__ int secp256k1_ec_seckey_tweak_add(unsigned char *seckey, const unsigned char *tweak) {
  secp256k1_scalar term;
  secp256k1_scalar sec;
  int ret = 0;
  int overflow = 0;
  secp256k1_scalar_set_b32(&term, tweak, &overflow);
  ret = secp256k1_scalar_set_b32_seckey(&sec, seckey);

  ret &= (!overflow) & secp256k1_eckey_privkey_tweak_add(&sec, &term);
  secp256k1_scalar secp256k1_scalar_zero = SECP256K1_SCALAR_CONST(0, 0, 0, 0, 0, 0, 0, 0);
  secp256k1_scalar_cmov(&sec, &secp256k1_scalar_zero, !ret);
  secp256k1_scalar_get_b32(seckey, &sec);

  secp256k1_scalar_clear(&sec);
  secp256k1_scalar_clear(&term);
  return ret;
}
*/
/*
__device__ __inline__ int secp256k1_pubkey_load(secp256k1_ge* ge, const secp256k1_pubkey* pubkey) {
  secp256k1_fe x, y;
  secp256k1_fe_set_b32(&x, pubkey->data);
  secp256k1_fe_set_b32(&y, pubkey->data + 32);
  secp256k1_ge_set_xy(ge, &x, &y);
   
  return 1;
}
*/
__device__ __inline__ int secp256k1_eckey_pubkey_serialize(secp256k1_ge *elem, unsigned char *pub, size_t *size, const int compressed) {
    if (secp256k1_ge_is_infinity(elem)) {
        return 0;
    }
    secp256k1_fe_normalize_var(&elem->x);
    secp256k1_fe_normalize_var(&elem->y);
    secp256k1_fe_get_b32(&pub[1], &elem->x);
    if (compressed) {
        *size = 33;
        pub[0] = secp256k1_fe_is_odd(&elem->y) ? SECP256K1_TAG_PUBKEY_ODD : SECP256K1_TAG_PUBKEY_EVEN;
    } else {
        *size = 65;
        pub[0] = SECP256K1_TAG_PUBKEY_UNCOMPRESSED;
        secp256k1_fe_get_b32(&pub[33], &elem->y);
    }
    return 1;
}
/*
__device__ __inline__ int secp256k1_ec_pubkey_serialize(unsigned char *output, size_t outputlen, const secp256k1_pubkey* pubkey, unsigned int flags) {
    secp256k1_ge Q;
    int ret = 0;
    memset(output, 0, outputlen);
    if (secp256k1_pubkey_load(&Q, pubkey)) {
        ret = secp256k1_eckey_pubkey_serialize(&Q, output, &outputlen, flags & SECP256K1_FLAGS_BIT_COMPRESSION);
    }
    return ret;
}
*/
/*
__device__ __inline__ int secp256k1_ec_pubkey_create(secp256k1_pubkey* pubkey, const unsigned char* seckey) {
    secp256k1_gej pj;
    secp256k1_ge p;
    secp256k1_scalar sec;
    secp256k1_scalar secp256k1_scalar_one = SECP256K1_SCALAR_CONST(0, 0, 0, 0, 0, 0, 0, 1);
    int ret = 0;

    memset(pubkey, 0, sizeof(*pubkey));

    ret = secp256k1_scalar_set_b32_seckey(&sec, seckey);

    secp256k1_scalar_cmov(&sec, &secp256k1_scalar_one, !ret);

    //secp256k1_ecmult_gen(&pj, &sec);
    secp256k1_ecmult_big(&pj, &sec);
    secp256k1_ge_set_gej(&p, &pj);
    secp256k1_pubkey_save(pubkey, &p);

    memczero(pubkey, sizeof(*pubkey), !ret);

    secp256k1_scalar_clear(&sec);
    return ret;
}*/