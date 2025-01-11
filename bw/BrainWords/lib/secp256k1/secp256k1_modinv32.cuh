/***********************************************************************
 * Copyright (c) 2020 Peter Dettman                                    *
 * Distributed under the MIT software license, see the accompanying    *
 * file COPYING or https://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef SECP256K1_MODINV32_IMPL_H
#define SECP256K1_MODINV32_IMPL_H

#include <stdlib.h>

 /* This file implements modular inversion based on the paper "Fast constant-time gcd computation and
  * modular inversion" by Daniel J. Bernstein and Bo-Yin Yang.
  *
  * For an explanation of the algorithm, see doc/safegcd_implementation.md. This file contains an
  * implementation for N=30, using 30-bit signed limbs represented as int32_t.
  */

/* Take as input a signed30 number in range (-2*modulus,modulus), and add a multiple of the modulus
 * to it to bring it to range [0,modulus). If sign < 0, the input will also be negated in the
 * process. The input must have limbs in range (-2^30,2^30). The output will have limbs in range
 * [0,2^30). */
__device__ __inline__ void secp256k1_modinv32_normalize_30(secp256k1_modinv32_signed30* __restrict__  r, int32_t sign, const secp256k1_modinv32_modinfo* __restrict__ modinfo) {
    int32_t r0 = r->v[0], r1 = r->v[1], r2 = r->v[2], r3 = r->v[3], r4 = r->v[4],
        r5 = r->v[5], r6 = r->v[6], r7 = r->v[7], r8 = r->v[8];
    int32_t cond_add, cond_negate;

    /* In a first step, add the modulus if the input is negative, and then negate if requested.
     * This brings r from range (-2*modulus,modulus) to range (-modulus,modulus). As all input
     * limbs are in range (-2^30,2^30), this cannot overflow an int32_t. Note that the right
     * shifts below are signed sign-extending shifts (see assumptions.h for tests that that is
     * indeed the behavior of the right shift operator). */
    cond_add = r8 >> 31;
    r0 += modinfo->modulus.v[0] & cond_add;
    r1 += modinfo->modulus.v[1] & cond_add;
    r2 += modinfo->modulus.v[2] & cond_add;
    r3 += modinfo->modulus.v[3] & cond_add;
    r4 += modinfo->modulus.v[4] & cond_add;
    r5 += modinfo->modulus.v[5] & cond_add;
    r6 += modinfo->modulus.v[6] & cond_add;
    r7 += modinfo->modulus.v[7] & cond_add;
    r8 += modinfo->modulus.v[8] & cond_add;
    cond_negate = sign >> 31;
    r0 = (r0 ^ cond_negate) - cond_negate;
    r1 = (r1 ^ cond_negate) - cond_negate;
    r2 = (r2 ^ cond_negate) - cond_negate;
    r3 = (r3 ^ cond_negate) - cond_negate;
    r4 = (r4 ^ cond_negate) - cond_negate;
    r5 = (r5 ^ cond_negate) - cond_negate;
    r6 = (r6 ^ cond_negate) - cond_negate;
    r7 = (r7 ^ cond_negate) - cond_negate;
    r8 = (r8 ^ cond_negate) - cond_negate;
    /* Propagate the top bits, to bring limbs back to range (-2^30,2^30). */
    r1 += r0 >> 30; r0 &= M30;
    r2 += r1 >> 30; r1 &= M30;
    r3 += r2 >> 30; r2 &= M30;
    r4 += r3 >> 30; r3 &= M30;
    r5 += r4 >> 30; r4 &= M30;
    r6 += r5 >> 30; r5 &= M30;
    r7 += r6 >> 30; r6 &= M30;
    r8 += r7 >> 30; r7 &= M30;

    /* In a second step add the modulus again if the result is still negative, bringing r to range
     * [0,modulus). */
    cond_add = r8 >> 31;
    r0 += modinfo->modulus.v[0] & cond_add;
    r1 += modinfo->modulus.v[1] & cond_add;
    r2 += modinfo->modulus.v[2] & cond_add;
    r3 += modinfo->modulus.v[3] & cond_add;
    r4 += modinfo->modulus.v[4] & cond_add;
    r5 += modinfo->modulus.v[5] & cond_add;
    r6 += modinfo->modulus.v[6] & cond_add;
    r7 += modinfo->modulus.v[7] & cond_add;
    r8 += modinfo->modulus.v[8] & cond_add;
    /* And propagate again. */
    r1 += r0 >> 30; r0 &= M30;
    r2 += r1 >> 30; r1 &= M30;
    r3 += r2 >> 30; r2 &= M30;
    r4 += r3 >> 30; r3 &= M30;
    r5 += r4 >> 30; r4 &= M30;
    r6 += r5 >> 30; r5 &= M30;
    r7 += r6 >> 30; r6 &= M30;
    r8 += r7 >> 30; r7 &= M30;

    r->v[0] = r0;
    r->v[1] = r1;
    r->v[2] = r2;
    r->v[3] = r3;
    r->v[4] = r4;
    r->v[5] = r5;
    r->v[6] = r6;
    r->v[7] = r7;
    r->v[8] = r8;


}

/* Data type for transition matrices (see section 3 of explanation).
 *
 * t = [ u  v ]
 *     [ q  r ]
 */
typedef struct {
    int32_t u, v, q, r;
} secp256k1_modinv32_trans2x2;

/* Compute the transition matrix and zeta for 30 divsteps.
 *
 * Input:  zeta: initial zeta
 *         f0:   bottom limb of initial f
 *         g0:   bottom limb of initial g
 * Output: t: transition matrix
 * Return: final zeta
 *
 * Implements the divsteps_n_matrix function from the explanation.
 */
__device__ __inline__ int32_t secp256k1_modinv32_divsteps_30(int32_t zeta, uint32_t f0, uint32_t g0, secp256k1_modinv32_trans2x2* __restrict__  t) {
    /* u,v,q,r are the elements of the transformation matrix being built up,
     * starting with the identity matrix. Semantically they are signed integers
     * in range [-2^30,2^30], but here represented as unsigned mod 2^32. This
     * permits left shifting (which is UB for negative numbers). The range
     * being inside [-2^31,2^31) means that casting to signed works correctly.
     */
    uint32_t u = 1, v = 0, q = 0, r = 1;
    uint32_t c1, c2, f = f0, g = g0, x, y, z;
    int i;
#pragma unroll
    for (i = 0; i < 30; ++i) {
        /* Compute conditional masks for (zeta < 0) and for (g & 1). */
        c1 = zeta >> 31;
        c2 = -(g & 1);
        /* Compute x,y,z, conditionally negated versions of f,u,v. */
        x = (f ^ c1) - c1;
        y = (u ^ c1) - c1;
        z = (v ^ c1) - c1;
        /* Conditionally add x,y,z to g,q,r. */
        g += x & c2;
        q += y & c2;
        r += z & c2;
        /* In what follows, c1 is a condition mask for (zeta < 0) and (g & 1). */
        c1 &= c2;
        /* Conditionally change zeta into -zeta-2 or zeta-1. */
        zeta = (zeta ^ c1) - 1;
        /* Conditionally add g,q,r to f,u,v. */
        f += g & c1;
        u += q & c1;
        v += r & c1;
        /* Shifts */
        g >>= 1;
        u <<= 1;
        v <<= 1;
    }
    /* Return data in t and return value. */
    t->u = (int32_t)u;
    t->v = (int32_t)v;
    t->q = (int32_t)q;
    t->r = (int32_t)r;
    return zeta;
}

/* Compute the transition matrix and eta for 30 divsteps (variable time).
 *
 * Input:  eta: initial eta
 *         f0:  bottom limb of initial f
 *         g0:  bottom limb of initial g
 * Output: t: transition matrix
 * Return: final eta
 *
 * Implements the divsteps_n_matrix_var function from the explanation.
 */
__device__ __inline__ int32_t secp256k1_modinv32_divsteps_30_var(int32_t eta, uint32_t f0, uint32_t g0, secp256k1_modinv32_trans2x2* __restrict__   t) {
    /* Transformation matrix; see comments in secp256k1_modinv32_divsteps_30. */
    uint32_t u = 1, v = 0, q = 0, r = 1;
    uint32_t f = f0, g = g0, m;
    uint16_t w;
    int i = 30, limit, zeros;

    for (;;) {
        /* Use a sentinel bit to count zeros only up to i. */
        zeros = secp256k1_ctz32_var(g | (UINT32_MAX << i));
        /* Perform zeros divsteps at once; they all just divide g by two. */
        g >>= zeros;
        u <<= zeros;
        v <<= zeros;
        eta -= zeros;
        i -= zeros;
        /* We're done once we've done 30 divsteps. */
        if (i == 0) break;
        /* If eta is negative, negate it and replace f,g with g,-f. */
        if (eta < 0) {
            uint32_t tmp;
            eta = -eta;
            tmp = f; f = g; g = -tmp;
            tmp = u; u = q; q = -tmp;
            tmp = v; v = r; r = -tmp;
        }
        /* eta is now >= 0. In what follows we're going to cancel out the bottom bits of g. No more
         * than i can be cancelled out (as we'd be done before that point), and no more than eta+1
         * can be done as its sign will flip once that happens. */
        limit = ((int)eta + 1) > i ? i : ((int)eta + 1);
        /* m is a mask for the bottom min(limit, 8) bits (our table only supports 8 bits). */
        m = (UINT32_MAX >> (32 - limit)) & 255U;
        /* Find what multiple of f must be added to g to cancel its bottom min(limit, 8) bits. */
        w = (g * inv256[(f >> 1) & 127]) & m;
        /* Do so. */
        g += f * w;
        q += u * w;
        r += v * w;
    }
    /* Return data in t and return value. */
    t->u = (int32_t)u;
    t->v = (int32_t)v;
    t->q = (int32_t)q;
    t->r = (int32_t)r;
    return eta;
}

/* Compute (t/2^30) * [d, e] mod modulus, where t is a transition matrix for 30 divsteps.
 *
 * On input and output, d and e are in range (-2*modulus,modulus). All output limbs will be in range
 * (-2^30,2^30).
 *
 * This implements the update_de function from the explanation.
 */
__device__ __inline__ void secp256k1_modinv32_update_de_30(secp256k1_modinv32_signed30* __restrict__  d, secp256k1_modinv32_signed30* __restrict__ e, const secp256k1_modinv32_trans2x2* __restrict__ t, const secp256k1_modinv32_modinfo* __restrict__ modinfo) {

    const int32_t u = t->u, v = t->v, q = t->q, r = t->r;
    int32_t di, ei, md, me, sd, se;
    int64_t cd, ce;
    int i;

    /* [md,me] start as zero; plus [u,q] if d is negative; plus [v,r] if e is negative. */
    sd = d->v[8] >> 31;
    se = e->v[8] >> 31;
    md = (u & sd) + (v & se);
    me = (q & sd) + (r & se);
    /* Begin computing t*[d,e]. */
    di = d->v[0];
    ei = e->v[0];
    cd = (int64_t)u * di + (int64_t)v * ei;
    ce = (int64_t)q * di + (int64_t)r * ei;
    /* Correct md,me so that t*[d,e]+modulus*[md,me] has 30 zero bottom bits. */
    md -= (modinfo->modulus_inv30 * (uint32_t)cd + md) & M30;
    me -= (modinfo->modulus_inv30 * (uint32_t)ce + me) & M30;
    /* Update the beginning of computation for t*[d,e]+modulus*[md,me] now md,me are known. */
    cd += (int64_t)modinfo->modulus.v[0] * md;
    ce += (int64_t)modinfo->modulus.v[0] * me;
    /* Verify that the low 30 bits of the computation are indeed zero, and then throw them away. */
    cd >>= 30;
    ce >>= 30;
    /* Now iteratively compute limb i=1..8 of t*[d,e]+modulus*[md,me], and store them in output
     * limb i-1 (shifting down by 30 bits). */
    for (i = 1; i < 9; ++i) {
        di = d->v[i];
        ei = e->v[i];
        cd += (int64_t)u * di + (int64_t)v * ei;
        ce += (int64_t)q * di + (int64_t)r * ei;
        cd += (int64_t)modinfo->modulus.v[i] * md;
        ce += (int64_t)modinfo->modulus.v[i] * me;
        d->v[i - 1] = (int32_t)cd & M30; cd >>= 30;
        e->v[i - 1] = (int32_t)ce & M30; ce >>= 30;
    }
    /* What remains is limb 9 of t*[d,e]+modulus*[md,me]; store it as output limb 8. */
    d->v[8] = (int32_t)cd;
    e->v[8] = (int32_t)ce;

}

/* Compute (t/2^30) * [f, g], where t is a transition matrix for 30 divsteps.
 *
 * This implements the update_fg function from the explanation.
 */
__device__ __inline__ void secp256k1_modinv32_update_fg_30(secp256k1_modinv32_signed30* __restrict__  f, secp256k1_modinv32_signed30* __restrict__  g, const secp256k1_modinv32_trans2x2* __restrict__  t) {

    const int32_t u = t->u, v = t->v, q = t->q, r = t->r;
    int32_t fi, gi;
    int64_t cf, cg;
    int i;
    /* Start computing t*[f,g]. */
    fi = f->v[0];
    gi = g->v[0];
    cf = (int64_t)u * fi + (int64_t)v * gi;
    cg = (int64_t)q * fi + (int64_t)r * gi;
    /* Verify that the bottom 30 bits of the result are zero, and then throw them away. */
    cf >>= 30;
    cg >>= 30;
    /* Now iteratively compute limb i=1..8 of t*[f,g], and store them in output limb i-1 (shifting
     * down by 30 bits). */
    for (i = 1; i < 9; ++i) {
        fi = f->v[i];
        gi = g->v[i];
        cf += (int64_t)u * fi + (int64_t)v * gi;
        cg += (int64_t)q * fi + (int64_t)r * gi;
        f->v[i - 1] = (int32_t)cf & M30; cf >>= 30;
        g->v[i - 1] = (int32_t)cg & M30; cg >>= 30;
    }
    /* What remains is limb 9 of t*[f,g]; store it as output limb 8. */
    f->v[8] = (int32_t)cf;
    g->v[8] = (int32_t)cg;
}

/* Compute (t/2^30) * [f, g], where t is a transition matrix for 30 divsteps.
 *
 * Version that operates on a variable number of limbs in f and g.
 *
 * This implements the update_fg function from the explanation in modinv64_impl.h.
 */
__device__ __inline__ void secp256k1_modinv32_update_fg_30_var(int len, secp256k1_modinv32_signed30* __restrict__  f, secp256k1_modinv32_signed30* __restrict__  g, const secp256k1_modinv32_trans2x2* __restrict__  t) {

    const int32_t u = t->u, v = t->v, q = t->q, r = t->r;
    int32_t fi, gi;
    int64_t cf, cg;
    int i;
    /* Start computing t*[f,g]. */
    fi = f->v[0];
    gi = g->v[0];
    cf = (int64_t)u * fi + (int64_t)v * gi;
    cg = (int64_t)q * fi + (int64_t)r * gi;
    /* Verify that the bottom 62 bits of the result are zero, and then throw them away. */
    cf >>= 30;
    cg >>= 30;
    /* Now iteratively compute limb i=1..len of t*[f,g], and store them in output limb i-1 (shifting
     * down by 30 bits). */
    for (i = 1; i < len; ++i) {
        fi = f->v[i];
        gi = g->v[i];
        cf += (int64_t)u * fi + (int64_t)v * gi;
        cg += (int64_t)q * fi + (int64_t)r * gi;
        f->v[i - 1] = (int32_t)cf & M30; cf >>= 30;
        g->v[i - 1] = (int32_t)cg & M30; cg >>= 30;
    }
    /* What remains is limb (len) of t*[f,g]; store it as output limb (len-1). */
    f->v[len - 1] = (int32_t)cf;
    g->v[len - 1] = (int32_t)cg;
}

/* Compute the inverse of x modulo modinfo->modulus, and replace x with it (constant time in x). */
__device__ __inline__ void secp256k1_modinv32(secp256k1_modinv32_signed30* __restrict__ x, const secp256k1_modinv32_modinfo* __restrict__  modinfo) {
    /* Start with d=0, e=1, f=modulus, g=x, zeta=-1. */
    secp256k1_modinv32_signed30 d = { {0} };
    secp256k1_modinv32_signed30 e = { {1} };
    secp256k1_modinv32_signed30 f = modinfo->modulus;
    secp256k1_modinv32_signed30 g = *x;
    int i;
    int32_t zeta = -1; /* zeta = -(delta+1/2); delta is initially 1/2. */

    /* Do 20 iterations of 30 divsteps each = 600 divsteps. 590 suffices for 256-bit inputs. */
#pragma unroll
    for (i = 0; i < 20; ++i) {
        /* Compute transition matrix and new zeta after 30 divsteps. */
        secp256k1_modinv32_trans2x2 t;
        zeta = secp256k1_modinv32_divsteps_30(zeta, f.v[0], g.v[0], &t);
        /* Update d,e using that transition matrix. */
        secp256k1_modinv32_update_de_30(&d, &e, &t, modinfo);
        /* Update f,g using that transition matrix. */

        secp256k1_modinv32_update_fg_30(&f, &g, &t);

    }

    /* At this point sufficient iterations have been performed that g must have reached 0
     * and (if g was not originally 0) f must now equal +/- GCD of the initial f, g
     * values i.e. +/- 1, and d now contains +/- the modular inverse. */


    /* Optionally negate d, normalize to [0,modulus), and return it. */
    secp256k1_modinv32_normalize_30(&d, f.v[8], modinfo);
    *x = d;
}

/* Compute the inverse of x modulo modinfo->modulus, and replace x with it (variable time). */
__device__ __inline__ void secp256k1_modinv32_var(secp256k1_modinv32_signed30* __restrict__  x, const secp256k1_modinv32_modinfo* __restrict__  modinfo) {
    /* Start with d=0, e=1, f=modulus, g=x, eta=-1. */
    secp256k1_modinv32_signed30 d = { {0, 0, 0, 0, 0, 0, 0, 0, 0} };
    secp256k1_modinv32_signed30 e = { {1, 0, 0, 0, 0, 0, 0, 0, 0} };
    secp256k1_modinv32_signed30 f = modinfo->modulus;
    secp256k1_modinv32_signed30 g = *x;
#ifdef VERIFY
    int i = 0;
#endif
    int j, len = 9;
    int32_t eta = -1; /* eta = -delta; delta is initially 1 (faster for the variable-time code) */
    int32_t cond, fn, gn;

    /* Do iterations of 30 divsteps each until g=0. */
    while (1) {
        /* Compute transition matrix and new eta after 30 divsteps. */
        secp256k1_modinv32_trans2x2 t;
        eta = secp256k1_modinv32_divsteps_30_var(eta, f.v[0], g.v[0], &t);
        /* Update d,e using that transition matrix. */
        secp256k1_modinv32_update_de_30(&d, &e, &t, modinfo);
        /* Update f,g using that transition matrix. */

        secp256k1_modinv32_update_fg_30_var(len, &f, &g, &t);
        /* If the bottom limb of g is 0, there is a chance g=0. */
        if (g.v[0] == 0) {
            cond = 0;
            /* Check if all other limbs are also 0. */
            for (j = 1; j < len; ++j) {
                cond |= g.v[j];
            }
            /* If so, we're done. */
            if (cond == 0) break;
        }

        /* Determine if len>1 and limb (len-1) of both f and g is 0 or -1. */
        fn = f.v[len - 1];
        gn = g.v[len - 1];
        cond = ((int32_t)len - 2) >> 31;
        cond |= fn ^ (fn >> 31);
        cond |= gn ^ (gn >> 31);
        /* If so, reduce length, propagating the sign of f and g's top limb into the one below. */
        if (cond == 0) {
            f.v[len - 2] |= (uint32_t)fn << 30;
            g.v[len - 2] |= (uint32_t)gn << 30;
            --len;
        }

    }

    /* Optionally negate d, normalize to [0,modulus), and return it. */
    secp256k1_modinv32_normalize_30(&d, f.v[len - 1], modinfo);
    *x = d;
}

#endif /* SECP256K1_MODINV32_IMPL_H */