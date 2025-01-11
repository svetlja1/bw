#if ARCH==32
#include "secp256k1_field.cuh"
#else
#include "secp256k1_field64.cuh"
#endif

__device__ __inline__ void secp256k1_ge_from_storage(secp256k1_ge * __restrict__ r, const secp256k1_ge_storage * __restrict__ a) {
    secp256k1_fe_from_storage(&r->x, &a->x);
    secp256k1_fe_from_storage(&r->y, &a->y);
    r->infinity = 0;
}

__device__ __inline__ int secp256k1_ge_set_xquad(secp256k1_ge * __restrict__ r, const secp256k1_fe * __restrict__ x) {
    secp256k1_fe x2, x3, c;
    r->x = *x;
    secp256k1_fe_sqr(&x2, x);
    secp256k1_fe_mul(&x3, x, &x2);
    r->infinity = 0;
    secp256k1_fe_set_int(&c, 7);
    secp256k1_fe_add(&c, &x3);
    return secp256k1_fe_sqrt(&r->y, &c);
}

__device__ __inline__ int secp256k1_ge_set_xo_var(secp256k1_ge * __restrict__ r, const secp256k1_fe * __restrict__ x, int odd) {
    if (!secp256k1_ge_set_xquad(r, x)) {
        return 0;
    }
    secp256k1_fe_normalize_var(&r->y);
    if (secp256k1_fe_is_odd(&r->y) != odd) {
        secp256k1_fe_negate(&r->y, &r->y, 1);
    }
    return 1;
}

__device__ __inline__ void secp256k1_gej_set_ge(secp256k1_gej * __restrict__ r, const secp256k1_ge * __restrict__ a) {
   r->infinity = a->infinity;
   r->x = a->x;
   r->y = a->y;
   secp256k1_fe_set_int(&r->z, 1);
}

__device__ __inline__ void secp256k1_gej_double_nonzero(secp256k1_gej * __restrict__ r, const secp256k1_gej * __restrict__ a) {
    /* Operations: 3 mul, 4 sqr, 8 add/half/mul_int/negate */
    secp256k1_fe l, s, t;

    r->infinity = a->infinity;

    /* Formula used:
     * L = (3/2) * X1^2
     * S = Y1^2
     * T = -X1*S
     * X3 = L^2 + 2*T
     * Y3 = -(L*(X3 + T) + S^2)
     * Z3 = Y1*Z1
     */

    secp256k1_fe_mul(&r->z, &a->z, &a->y); /* Z3 = Y1*Z1 (1) */
    secp256k1_fe_sqr(&s, &a->y);           /* S = Y1^2 (1) */
    secp256k1_fe_sqr(&l, &a->x);           /* L = X1^2 (1) */
    secp256k1_fe_mul_int(&l, 3);           /* L = 3*X1^2 (3) */
    secp256k1_fe_half(&l);                 /* L = 3/2*X1^2 (2) */
    secp256k1_fe_negate(&t, &s, 1);        /* T = -S (2) */
    secp256k1_fe_mul(&t, &t, &a->x);       /* T = -X1*S (1) */
    secp256k1_fe_sqr(&r->x, &l);           /* X3 = L^2 (1) */
    secp256k1_fe_add(&r->x, &t);           /* X3 = L^2 + T (2) */
    secp256k1_fe_add(&r->x, &t);           /* X3 = L^2 + 2*T (3) */
    secp256k1_fe_sqr(&s, &s);              /* S' = S^2 (1) */
    secp256k1_fe_add(&t, &r->x);           /* T' = X3 + T (4) */
    secp256k1_fe_mul(&r->y, &t, &l);       /* Y3 = L*(X3 + T) (1) */
    secp256k1_fe_add(&r->y, &s);           /* Y3 = L*(X3 + T) + S^2 (2) */
    secp256k1_fe_negate(&r->y, &r->y, 2);  /* Y3 = -(L*(X3 + T) + S^2) (3) */
}

__device__ __inline__ void secp256k1_gej_neg(secp256k1_gej * __restrict__ r, const secp256k1_gej * __restrict__ a) {
    r->infinity = a->infinity;
    r->x = a->x;
    r->y = a->y;
    r->z = a->z;
    secp256k1_fe_normalize_weak(&r->y);
    secp256k1_fe_negate(&r->y, &r->y, 1);
}

__device__ __inline__ void secp256k1_gej_double_var(secp256k1_gej * __restrict__ r, const secp256k1_gej * __restrict__ a, secp256k1_fe *rzr) {
    /** For secp256k1, 2Q is infinity if and only if Q is infinity. This is because if 2Q = infinity,
     *  Q must equal -Q, or that Q.y == -(Q.y), or Q.y is 0. For a point on y^2 = x^3 + 7 to have
     *  y=0, x^3 must be -7 mod p. However, -7 has no cube root mod p.
     *
     *  Having said this, if this function receives a point on a sextic twist, e.g. by
     *  a fault attack, it is possible for y to be 0. This happens for y^2 = x^3 + 6,
     *  since -6 does have a cube root mod p. For this point, this function will not set
     *  the infinity flag even though the point doubles to infinity, and the result
     *  point will be gibberish (z = 0 but infinity = 0).
     */
    if (a->infinity) {
        r->infinity = 1;
        if (rzr == NULL) {

        } else {
          secp256k1_fe_set_int(rzr, 1);
        }
        return;
    }

    if (rzr != NULL) {
        *rzr = a->y;
        secp256k1_fe_normalize_weak(rzr);
        //secp256k1_fe_mul_int(rzr, 2);
    }

    secp256k1_gej_double_nonzero(r, a);
}

__device__ __inline__ void secp256k1_gej_add_var(secp256k1_gej * __restrict__ r, const secp256k1_gej * __restrict__ a, const secp256k1_gej * __restrict__ b, secp256k1_fe *rzr) {
    /* Operations: 12 mul, 4 sqr, 2 normalize, 12 mul_int/add/negate */
    secp256k1_fe z22, z12, u1, u2, s1, s2, h, i, i2, h2, h3, t;

    if (a->infinity) {
        *r = *b;
        return;
    }

    if (b->infinity) {
        if (rzr != NULL) {
            secp256k1_fe_set_int(rzr, 1);
        }
        *r = *a;
        return;
    }

    r->infinity = 0;
    secp256k1_fe_sqr(&z22, &b->z);
    secp256k1_fe_sqr(&z12, &a->z);
    secp256k1_fe_mul(&u1, &a->x, &z22);
    secp256k1_fe_mul(&u2, &b->x, &z12);
    secp256k1_fe_mul(&s1, &a->y, &z22); 
    secp256k1_fe_mul(&s1, &s1, &b->z);
    secp256k1_fe_mul(&s2, &b->y, &z12); 
    secp256k1_fe_mul(&s2, &s2, &a->z);
    secp256k1_fe_negate(&h, &u1, 1); 
    secp256k1_fe_add(&h, &u2);
    secp256k1_fe_negate(&i, &s1, 1); 
    secp256k1_fe_add(&i, &s2);
    if (secp256k1_fe_normalizes_to_zero_var(&h)) {
        if (secp256k1_fe_normalizes_to_zero_var(&i)) {
            secp256k1_gej_double_var(r, a, rzr);
        } else {
            if (rzr != NULL) {
                secp256k1_fe_set_int(rzr, 0);
            }
            r->infinity = 1;
        }
        return;
    }
    secp256k1_fe_sqr(&i2, &i);
    secp256k1_fe_sqr(&h2, &h);
    secp256k1_fe_mul(&h3, &h, &h2);
    secp256k1_fe_mul(&h, &h, &b->z);
    if (rzr != NULL) {
        *rzr = h;
    }
    secp256k1_fe_mul(&r->z, &a->z, &h);
    secp256k1_fe_mul(&t, &u1, &h2);
    r->x = t; secp256k1_fe_mul_int(&r->x, 2); secp256k1_fe_add(&r->x, &h3); secp256k1_fe_negate(&r->x, &r->x, 3); secp256k1_fe_add(&r->x, &i2);
    secp256k1_fe_negate(&r->y, &r->x, 5); secp256k1_fe_add(&r->y, &t); secp256k1_fe_mul(&r->y, &r->y, &i);
    secp256k1_fe_mul(&h3, &h3, &s1); secp256k1_fe_negate(&h3, &h3, 1);
    secp256k1_fe_add(&r->y, &h3);
}

__device__ __inline__ void secp256k1_gej_set_infinity(secp256k1_gej* __restrict__ r) {
    r->infinity = 1;
    secp256k1_fe_clear(&r->x);
    secp256k1_fe_clear(&r->y);
    secp256k1_fe_clear(&r->z);
}

__device__ __inline__ void secp256k1_gej_add_ge_var(secp256k1_gej * __restrict__ r, const secp256k1_gej * __restrict__ a, const secp256k1_ge * __restrict__ b, secp256k1_fe *rzr) {
    /* 8 mul, 3 sqr, 13 add/negate/normalize_weak/normalizes_to_zero (ignoring special cases) */
    secp256k1_fe z12, u1, u2, s1, s2, h, i, h2, h3, t;
    if (a->infinity) {
        //VERIFY_CHECK(rzr == NULL);
        secp256k1_gej_set_ge(r, b);
        return;
    }
    if (b->infinity) {
        if (rzr != NULL) {
            secp256k1_fe_set_int(rzr, 1);
        }
        *r = *a;
        return;
    }

    secp256k1_fe_sqr(&z12, &a->z);
    u1 = a->x; secp256k1_fe_normalize_weak(&u1);
    secp256k1_fe_mul(&u2, &b->x, &z12);
    s1 = a->y; secp256k1_fe_normalize_weak(&s1);
    secp256k1_fe_mul(&s2, &b->y, &z12); secp256k1_fe_mul(&s2, &s2, &a->z);
    secp256k1_fe_negate(&h, &u1, 1); secp256k1_fe_add(&h, &u2);
    secp256k1_fe_negate(&i, &s2, 1); secp256k1_fe_add(&i, &s1);
    if (secp256k1_fe_normalizes_to_zero_var(&h)) {
        if (secp256k1_fe_normalizes_to_zero_var(&i)) {
            secp256k1_gej_double_var(r, a, rzr);
        } else {
            if (rzr != NULL) {
                secp256k1_fe_set_int(rzr, 0);
            }
            secp256k1_gej_set_infinity(r);
        }
        return;
    }

    r->infinity = 0;
    if (rzr != NULL) {
        *rzr = h;
    }
    secp256k1_fe_mul(&r->z, &a->z, &h);

    secp256k1_fe_sqr(&h2, &h);
    secp256k1_fe_negate(&h2, &h2, 1);
    secp256k1_fe_mul(&h3, &h2, &h);
    secp256k1_fe_mul(&t, &u1, &h2);

    secp256k1_fe_sqr(&r->x, &i);
    secp256k1_fe_add(&r->x, &h3);
    secp256k1_fe_add(&r->x, &t);
    secp256k1_fe_add(&r->x, &t);

    secp256k1_fe_add(&t, &r->x);
    secp256k1_fe_mul(&r->y, &t, &i);
    secp256k1_fe_mul(&h3, &h3, &s1);
    secp256k1_fe_add(&r->y, &h3);
}

__device__ __inline__ void secp256k1_gej_add_ge(secp256k1_gej * __restrict__ r, const secp256k1_gej * __restrict__ a, const secp256k1_ge * __restrict__ b) {
    /* Operations: 7 mul, 5 sqr, 4 normalize, 21 mul_int/add/negate/cmov */
    secp256k1_fe fe_1 = SECP256K1_FE_CONST(0, 0, 0, 0, 0, 0, 0, 1);
    secp256k1_fe zz, u1, u2, s1, s2, t, tt, m, n, q, rr;
    secp256k1_fe m_alt, rr_alt;
    int infinity, degenerate;
    
    secp256k1_fe_sqr(&zz, &a->z);                       /* z = Z1^2 */
    u1 = a->x; secp256k1_fe_normalize_weak(&u1);        /* u1 = U1 = X1*Z2^2 (1) */
    secp256k1_fe_mul(&u2, &b->x, &zz);                  /* u2 = U2 = X2*Z1^2 (1) */
    s1 = a->y; secp256k1_fe_normalize_weak(&s1);        /* s1 = S1 = Y1*Z2^3 (1) */
    secp256k1_fe_mul(&s2, &b->y, &zz);                  /* s2 = Y2*Z1^2 (1) */
    secp256k1_fe_mul(&s2, &s2, &a->z);                  /* s2 = S2 = Y2*Z1^3 (1) */
    t = u1; secp256k1_fe_add(&t, &u2);                  /* t = T = U1+U2 (2) */
    m = s1; secp256k1_fe_add(&m, &s2);                  /* m = M = S1+S2 (2) */
    secp256k1_fe_sqr(&rr, &t);                          /* rr = T^2 (1) */
    secp256k1_fe_negate(&m_alt, &u2, 1);                /* Malt = -X2*Z1^2 */
    secp256k1_fe_mul(&tt, &u1, &m_alt);                 /* tt = -U1*U2 (2) */
    secp256k1_fe_add(&rr, &tt);                         /* rr = R = T^2-U1*U2 (3) */
    /** If lambda = R/M = 0/0 we have a problem (except in the "trivial"
     *  case that Z = z1z2 = 0, and this is special-cased later on). */
    degenerate = secp256k1_fe_normalizes_to_zero(&m) &
                 secp256k1_fe_normalizes_to_zero(&rr);
    /* This only occurs when y1 == -y2 and x1^3 == x2^3, but x1 != x2.
     * This means either x1 == beta*x2 or beta*x1 == x2, where beta is
     * a nontrivial cube root of one. In either case, an alternate
     * non-indeterminate expression for lambda is (y1 - y2)/(x1 - x2),
     * so we set R/M equal to this. */
    rr_alt = s1;
    secp256k1_fe_mul_int(&rr_alt, 2);       /* rr = Y1*Z2^3 - Y2*Z1^3 (2) */
    secp256k1_fe_add(&m_alt, &u1);          /* Malt = X1*Z2^2 - X2*Z1^2 */

    secp256k1_fe_cmov(&rr_alt, &rr, !degenerate);
    secp256k1_fe_cmov(&m_alt, &m, !degenerate);
    /* Now Ralt / Malt = lambda and is guaranteed not to be 0/0.
     * From here on out Ralt and Malt represent the numerator
     * and denominator of lambda; R and M represent the explicit
     * expressions x1^2 + x2^2 + x1x2 and y1 + y2. */
    secp256k1_fe_sqr(&n, &m_alt);                       /* n = Malt^2 (1) */
    secp256k1_fe_mul(&q, &n, &t);                       /* q = Q = T*Malt^2 (1) */
    /* These two lines use the observation that either M == Malt or M == 0,
     * so M^3 * Malt is either Malt^4 (which is computed by squaring), or
     * zero (which is "computed" by cmov). So the cost is one squaring
     * versus two multiplications. */
    secp256k1_fe_sqr(&n, &n);
    secp256k1_fe_cmov(&n, &m, degenerate);              /* n = M^3 * Malt (2) */
    secp256k1_fe_sqr(&t, &rr_alt);                      /* t = Ralt^2 (1) */
    secp256k1_fe_mul(&r->z, &a->z, &m_alt);             /* r->z = Malt*Z (1) */
    infinity = secp256k1_fe_normalizes_to_zero(&r->z) * (1 - a->infinity);
    secp256k1_fe_mul_int(&r->z, 2);                     /* r->z = Z3 = 2*Malt*Z (2) */
    secp256k1_fe_negate(&q, &q, 1);                     /* q = -Q (2) */
    secp256k1_fe_add(&t, &q);                           /* t = Ralt^2-Q (3) */
    secp256k1_fe_normalize_weak(&t);
    r->x = t;                                           /* r->x = Ralt^2-Q (1) */
    secp256k1_fe_mul_int(&t, 2);                        /* t = 2*x3 (2) */
    secp256k1_fe_add(&t, &q);                           /* t = 2*x3 - Q: (4) */
    secp256k1_fe_mul(&t, &t, &rr_alt);                  /* t = Ralt*(2*x3 - Q) (1) */
    secp256k1_fe_add(&t, &n);                           /* t = Ralt*(2*x3 - Q) + M^3*Malt (3) */
    secp256k1_fe_negate(&r->y, &t, 3);                  /* r->y = Ralt*(Q - 2x3) - M^3*Malt (4) */
    secp256k1_fe_normalize_weak(&r->y);
    secp256k1_fe_mul_int(&r->x, 4);                     /* r->x = X3 = 4*(Ralt^2-Q) */
    secp256k1_fe_mul_int(&r->y, 4);                     /* r->y = Y3 = 4*Ralt*(Q - 2x3) - 4*M^3*Malt (4) */

    /** In case a->infinity == 1, replace r with (b->x, b->y, 1). */
    secp256k1_fe_cmov(&r->x, &b->x, a->infinity);
    secp256k1_fe_cmov(&r->y, &b->y, a->infinity);
    secp256k1_fe_cmov(&r->z, &fe_1, a->infinity);
    r->infinity = infinity;
}

__device__ __inline__ void secp256k1_ge_clear(secp256k1_ge *r) {
    r->infinity = 0;
    secp256k1_fe_clear(&r->x);
    secp256k1_fe_clear(&r->y);
}



__device__ __inline__ void secp256k1_ge_set_gej(secp256k1_ge * __restrict__ r, secp256k1_gej * __restrict__ a) {
    secp256k1_fe z2, z3;
    r->infinity = a->infinity;
    secp256k1_fe_inv(&a->z, &a->z);
    secp256k1_fe_sqr(&z2, &a->z);
    secp256k1_fe_mul(&z3, &a->z, &z2);
    secp256k1_fe_mul(&a->x, &a->x, &z2);
    secp256k1_fe_mul(&a->y, &a->y, &z3);
    secp256k1_fe_set_int(&a->z, 1);
    r->x = a->x;
    r->y = a->y;
}

__device__ __inline__ void secp256k1_ge_set_gej_zinv(secp256k1_ge * __restrict__ r, const secp256k1_gej * __restrict__ a, const secp256k1_fe * __restrict__ zi) {
  secp256k1_fe zi2;
  secp256k1_fe zi3;
  secp256k1_fe_sqr(&zi2, zi);
  secp256k1_fe_mul(&zi3, &zi2, zi);
  secp256k1_fe_mul(&r->x, &a->x, &zi2);
  secp256k1_fe_mul(&r->y, &a->y, &zi3);
  r->infinity = a->infinity;
}



__device__ __inline__ void secp256k1_ge_to_storage(secp256k1_ge_storage * __restrict__ r, const secp256k1_ge * __restrict__ a) {
    secp256k1_fe x, y;
    x = a->x;
    secp256k1_fe_normalize(&x);
    y = a->y;
    secp256k1_fe_normalize(&y);
    secp256k1_fe_to_storage(&r->x, &x);
    secp256k1_fe_to_storage(&r->y, &y);
}

__device__ __inline__ void secp256k1_ge_set_xy(secp256k1_ge * __restrict__ r, const secp256k1_fe * __restrict__ x, const secp256k1_fe * __restrict__ y) {
    r->infinity = 0;
    r->x = *x;
    r->y = *y;
}

__device__ __inline__ int secp256k1_ge_is_infinity(const secp256k1_ge * __restrict__ a) {
    return a->infinity;
}

__device__ __inline__ void secp256k1_ge_set_infinity(secp256k1_ge* __restrict__  r) {
    r->infinity = 1;
    secp256k1_fe_clear(&r->x);
    secp256k1_fe_clear(&r->y);
}

__device__ __inline__ void secp256k1_ge_set_all_gej_var(secp256k1_ge* __restrict__ r, const secp256k1_gej* __restrict__ a, size_t len) {
    secp256k1_fe u;
    size_t i;
    size_t last_i = SIZE_MAX;

    for (i = 0; i < len; i++) {
        if (a[i].infinity) {
            secp256k1_ge_set_infinity(&r[i]);
        }
        else {
            /* Use destination's x coordinates as scratch space */
            if (last_i == SIZE_MAX) {
                r[i].x = a[i].z;
            }
            else {
                secp256k1_fe_mul(&r[i].x, &r[last_i].x, &a[i].z);
            }
            last_i = i;
        }
    }
    if (last_i == SIZE_MAX) {
        return;
    }
    secp256k1_fe_inv_var(&u, &r[last_i].x);

    i = last_i;
    while (i > 0) {
        i--;
        if (!a[i].infinity) {
            secp256k1_fe_mul(&r[last_i].x, &r[i].x, &u);
            secp256k1_fe_mul(&u, &u, &a[last_i].z);
            last_i = i;
        }
    }
    //VERIFY_CHECK(!a[last_i].infinity);
    r[last_i].x = u;

    for (i = 0; i < len; i++) {
        if (!a[i].infinity) {
            secp256k1_ge_set_gej_zinv(&r[i], &a[i], &r[i].x);
        }
    }
}


__device__ __inline__  void secp256k1_ge_neg(secp256k1_ge* __restrict__ r, const secp256k1_ge* __restrict__ a) {
    *r = *a;
    secp256k1_fe_normalize_weak(&r->y);
    secp256k1_fe_negate(&r->y, &r->y, 1);
}