#include <stdlib.h>
#include <stdio.h>

#include "Kernel.cuh"

#include "lib/hash/GPUHash.cuh"
//#include "lib/hash/sha3_ver3.cuh"

#ifdef  ECMULT_BIG_TABLE
#include "lib/secp256k1/secp256k1_precomp_big_custom.cuh"
#else
#if ECMULT_GEN_PREC_BITS==2
#include "lib/secp256k1/secp256k1_prec.cuh"
#elif ECMULT_GEN_PREC_BITS == 4
#include "lib/secp256k1/secp256k1_prec4.cuh"
#elif ECMULT_GEN_PREC_BITS == 8
#include "lib/secp256k1/secp256k1_prec8.cuh"
#else
#include "lib/secp256k1/secp256k1_prec_custom.cuh"
#endif
#endif //  ECMULT_BIG_TABLE
#include "lib/secp256k1/secp256k1_group.cuh"
#if ARCH==32
#include "lib/secp256k1/secp256k1_scalar.cuh"
#else
#include "lib/secp256k1/secp256k1_scalar64.cuh"
#endif
#include "lib/secp256k1/secp256k1.cuh"
#include "lib/secp256k1/secp256k1_batch_impl.cuh"

__device__ __constant__ char PREFIX[PREFIX_MAX_LEN];

__device__ __constant__ uint32_t _NUM_TARGET_HASHES[1];
__device__ __constant__ uint32_t* _BLOOM_FILTER[1];
__device__ __constant__ uint32_t _BLOOM_FILTER_MASK[1];
__device__ __constant__ uint64_t _BLOOM_FILTER_MASK64[1];
__device__ __constant__ uint32_t _USE_BLOOM_FILTER[1];

__device__ __constant__ int ITERATION[1];
__device__ __constant__ int SHA_LEVEL[1];
__device__ __constant__ uint32_t SHA_PRE[8];

__global__ void worker(bool* isResult, bool* buffResult, const char* __restrict__ lines, const uint32_t* __restrict__ indexes, const uint32_t indexes_size, const secp256k1_ge_storage* __restrict__ precPtr, const size_t precPitch) {
	volatile int64_t tIx = threadIdx.x + blockIdx.x * blockDim.x;
	if (tIx >= indexes_size) {
		return;
	}
	const int _ITER = ITERATION[0];
	char* stringStart;
	uint64_t len = 0;
	uint32_t d_hash[8];
	unsigned char pubKeys[THREAD_STEPS * 65];
	unsigned char prvKeys[THREAD_STEPS * 32];
	uint32_t hash160[5];
	uint64_t starter = THREAD_STEPS * tIx;

	int privKeyIx = 0;
	while (privKeyIx < THREAD_STEPS && starter < indexes_size) {
		if (privKeyIx == 0) {
			if (starter == 0) {
				stringStart = (char*)lines;
				len = indexes[0];
			}
			else {
				stringStart = ((char*)lines + indexes[starter - 1]);
				len = indexes[starter] - indexes[starter - 1];
			}
		}
		else {
			len = indexes[starter] - indexes[starter - 1];
		}
		char toHash[63] = { 0 };
		size_t ix = 0;
		for (; ix < len; stringStart++) {
			toHash[ix++] = *stringStart;
			if (ix == 63) {
				break;
			}
		}
		toHash[len] = 0x80;
		_GetHashSha256(toHash, len, d_hash, _ITER);
		int pk4 = privKeyIx << 5;
		prvKeys[pk4] = (unsigned char)__byte_perm(d_hash[0], 0, 0x3);
		prvKeys[pk4 + 1] = (unsigned char)__byte_perm(d_hash[0], 0, 0x2);
		prvKeys[pk4 + 2] = (unsigned char)__byte_perm(d_hash[0], 0, 0x1);
		prvKeys[pk4 + 3] = (unsigned char)__byte_perm(d_hash[0], 0, 0x0);
		prvKeys[pk4 + 4] = (unsigned char)__byte_perm(d_hash[1], 0, 0x3);
		prvKeys[pk4 + 5] = (unsigned char)__byte_perm(d_hash[1], 0, 0x2);
		prvKeys[pk4 + 6] = (unsigned char)__byte_perm(d_hash[1], 0, 0x1);
		prvKeys[pk4 + 7] = (unsigned char)__byte_perm(d_hash[1], 0, 0x0);
		prvKeys[pk4 + 8] = (unsigned char)__byte_perm(d_hash[2], 0, 0x3);
		prvKeys[pk4 + 9] = (unsigned char)__byte_perm(d_hash[2], 0, 0x2);
		prvKeys[pk4 + 10] = (unsigned char)__byte_perm(d_hash[2], 0, 0x1);
		prvKeys[pk4 + 11] = (unsigned char)__byte_perm(d_hash[2], 0, 0x0);
		prvKeys[pk4 + 12] = (unsigned char)__byte_perm(d_hash[3], 0, 0x3);
		prvKeys[pk4 + 13] = (unsigned char)__byte_perm(d_hash[3], 0, 0x2);
		prvKeys[pk4 + 14] = (unsigned char)__byte_perm(d_hash[3], 0, 0x1);
		prvKeys[pk4 + 15] = (unsigned char)__byte_perm(d_hash[3], 0, 0x0);
		prvKeys[pk4 + 16] = (unsigned char)__byte_perm(d_hash[4], 0, 0x3);
		prvKeys[pk4 + 17] = (unsigned char)__byte_perm(d_hash[4], 0, 0x2);
		prvKeys[pk4 + 18] = (unsigned char)__byte_perm(d_hash[4], 0, 0x1);
		prvKeys[pk4 + 19] = (unsigned char)__byte_perm(d_hash[4], 0, 0x0);
		prvKeys[pk4 + 20] = (unsigned char)__byte_perm(d_hash[5], 0, 0x3);
		prvKeys[pk4 + 21] = (unsigned char)__byte_perm(d_hash[5], 0, 0x2);
		prvKeys[pk4 + 22] = (unsigned char)__byte_perm(d_hash[5], 0, 0x1);
		prvKeys[pk4 + 23] = (unsigned char)__byte_perm(d_hash[5], 0, 0x0);
		prvKeys[pk4 + 24] = (unsigned char)__byte_perm(d_hash[6], 0, 0x3);
		prvKeys[pk4 + 25] = (unsigned char)__byte_perm(d_hash[6], 0, 0x2);
		prvKeys[pk4 + 26] = (unsigned char)__byte_perm(d_hash[6], 0, 0x1);
		prvKeys[pk4 + 27] = (unsigned char)__byte_perm(d_hash[6], 0, 0x0);
		prvKeys[pk4 + 28] = (unsigned char)__byte_perm(d_hash[7], 0, 0x3);
		prvKeys[pk4 + 29] = (unsigned char)__byte_perm(d_hash[7], 0, 0x2);
		prvKeys[pk4 + 30] = (unsigned char)__byte_perm(d_hash[7], 0, 0x1);
		prvKeys[pk4 + 31] = (unsigned char)__byte_perm(d_hash[7], 0, 0x0);
		privKeyIx++;
		starter++;
	}


	secp256k1_ec_pubkey_create_serialized_batch_myunsafe(pubKeys, prvKeys, THREAD_STEPS, precPtr, precPitch);
	//printf("   pubkey %d\n", pubKeys[1]);

	int keyLenSkip;
	int pkField = 0;
	starter = THREAD_STEPS * tIx;
	for (pkField = 0; pkField < THREAD_STEPS; pkField++, starter++) {
		keyLenSkip = 65 * pkField;
		_GetHash160(pubKeys, keyLenSkip, (uint8_t*)hash160);
		if (checkHash(hash160)) {
			buffResult[starter] = true;
			isResult[0] = true;
			return;
		}
		keyLenSkip = 65 * pkField;
		_GetHash160Comp(pubKeys, keyLenSkip, (uint8_t*)hash160);
		if (checkHash(hash160)) {
			buffResult[starter] = true;
			isResult[0] = true;
			return;
		}
		keyLenSkip = 65 * pkField;
		_GetHash160P2SHCompFromHash(hash160, hash160);
		if (checkHash(hash160)) {
			buffResult[starter] = true;
			isResult[0] = true;
			return;
		}
	}
}
__global__ void workerEth(bool* isResult, bool* buffResult, const char* __restrict__ lines, const uint32_t* __restrict__ indexes, const uint32_t indexes_size, const short workType, const bool sha256hash, const secp256k1_ge_storage* __restrict__ precPtr, const size_t precPitch) {
	volatile int64_t tIx = threadIdx.x + blockIdx.x * blockDim.x;
	if (tIx >= indexes_size) {
		return;
	}
	const int _ITER = ITERATION[0];
	char* stringStart;
	int len = 0;
	unsigned char d_hash[32];	
	unsigned char pubKeys[THREAD_STEPS * 65];
	unsigned char prvKeys[THREAD_STEPS * 32];
	char toHashWords[THREAD_STEPS][256];
	int toHashWordsLen[THREAD_STEPS];
	
	//read dictionary
	int privKeyIx = 0;
	uint64_t starter = THREAD_STEPS * tIx;
	int wordsRead = 0;
	while (privKeyIx < THREAD_STEPS && starter < indexes_size) {
		if (privKeyIx == 0) {
			if (starter == 0) {
				stringStart = (char*)lines;
				len = indexes[0];
			}
			else {
				stringStart = ((char*)lines + indexes[starter - 1]);
				len = indexes[starter] - indexes[starter - 1];
			}
		}
		else {
			len = indexes[starter] - indexes[starter - 1];
		}
		size_t ix = 0;
		for (; ix < len; stringStart++) {
			toHashWords[privKeyIx][ix++] = *stringStart;
		}
		toHashWordsLen[privKeyIx] = len;
		wordsRead++;
		starter++;
		privKeyIx++;
	}

	uint32_t d_hashsha[8];

	for (short turn = (workType%2==0)?2:1; turn <= workType; turn*=2) {				
		privKeyIx = 0;
		while (privKeyIx < wordsRead) {
			if (sha256hash) {
				toHashWords[privKeyIx][toHashWordsLen[privKeyIx]] = 0x80;
				_GetHashSha256(toHashWords[privKeyIx], toHashWordsLen[privKeyIx], d_hashsha, _ITER);
				int pk4 = privKeyIx << 5;
				prvKeys[pk4] = (unsigned char)__byte_perm(d_hashsha[0], 0, 0x3);
				prvKeys[pk4 + 1] = (unsigned char)__byte_perm(d_hashsha[0], 0, 0x2);
				prvKeys[pk4 + 2] = (unsigned char)__byte_perm(d_hashsha[0], 0, 0x1);
				prvKeys[pk4 + 3] = (unsigned char)__byte_perm(d_hashsha[0], 0, 0x0);
				prvKeys[pk4 + 4] = (unsigned char)__byte_perm(d_hashsha[1], 0, 0x3);
				prvKeys[pk4 + 5] = (unsigned char)__byte_perm(d_hashsha[1], 0, 0x2);
				prvKeys[pk4 + 6] = (unsigned char)__byte_perm(d_hashsha[1], 0, 0x1);
				prvKeys[pk4 + 7] = (unsigned char)__byte_perm(d_hashsha[1], 0, 0x0);
				prvKeys[pk4 + 8] = (unsigned char)__byte_perm(d_hashsha[2], 0, 0x3);
				prvKeys[pk4 + 9] = (unsigned char)__byte_perm(d_hashsha[2], 0, 0x2);
				prvKeys[pk4 + 10] = (unsigned char)__byte_perm(d_hashsha[2], 0, 0x1);
				prvKeys[pk4 + 11] = (unsigned char)__byte_perm(d_hashsha[2], 0, 0x0);
				prvKeys[pk4 + 12] = (unsigned char)__byte_perm(d_hashsha[3], 0, 0x3);
				prvKeys[pk4 + 13] = (unsigned char)__byte_perm(d_hashsha[3], 0, 0x2);
				prvKeys[pk4 + 14] = (unsigned char)__byte_perm(d_hashsha[3], 0, 0x1);
				prvKeys[pk4 + 15] = (unsigned char)__byte_perm(d_hashsha[3], 0, 0x0);
				prvKeys[pk4 + 16] = (unsigned char)__byte_perm(d_hashsha[4], 0, 0x3);
				prvKeys[pk4 + 17] = (unsigned char)__byte_perm(d_hashsha[4], 0, 0x2);
				prvKeys[pk4 + 18] = (unsigned char)__byte_perm(d_hashsha[4], 0, 0x1);
				prvKeys[pk4 + 19] = (unsigned char)__byte_perm(d_hashsha[4], 0, 0x0);
				prvKeys[pk4 + 20] = (unsigned char)__byte_perm(d_hashsha[5], 0, 0x3);
				prvKeys[pk4 + 21] = (unsigned char)__byte_perm(d_hashsha[5], 0, 0x2);
				prvKeys[pk4 + 22] = (unsigned char)__byte_perm(d_hashsha[5], 0, 0x1);
				prvKeys[pk4 + 23] = (unsigned char)__byte_perm(d_hashsha[5], 0, 0x0);
				prvKeys[pk4 + 24] = (unsigned char)__byte_perm(d_hashsha[6], 0, 0x3);
				prvKeys[pk4 + 25] = (unsigned char)__byte_perm(d_hashsha[6], 0, 0x2);
				prvKeys[pk4 + 26] = (unsigned char)__byte_perm(d_hashsha[6], 0, 0x1);
				prvKeys[pk4 + 27] = (unsigned char)__byte_perm(d_hashsha[6], 0, 0x0);
				prvKeys[pk4 + 28] = (unsigned char)__byte_perm(d_hashsha[7], 0, 0x3);
				prvKeys[pk4 + 29] = (unsigned char)__byte_perm(d_hashsha[7], 0, 0x2);
				prvKeys[pk4 + 30] = (unsigned char)__byte_perm(d_hashsha[7], 0, 0x1);
				prvKeys[pk4 + 31] = (unsigned char)__byte_perm(d_hashsha[7], 0, 0x0);
			}
			else {
				keccak(toHashWords[privKeyIx], toHashWordsLen[privKeyIx], d_hash, 32);
				if (turn == 2) {
#pragma unroll 2030
					for (int camp2Turns = 0; camp2Turns < 2030; camp2Turns++) {
#pragma unroll 32
						for (int xxx = 0; xxx < 32; xxx++) {
							toHashWords[privKeyIx][xxx] = d_hash[xxx];
						}
						keccak(toHashWords[privKeyIx], 32, d_hash, 32);
					}
				}
#pragma unroll 32
				for (int pkField = 32 * privKeyIx, len = 0; len < 32; len++) {
					prvKeys[pkField++] = d_hash[len];
				}
			}			
			privKeyIx++;
		}

		secp256k1_ec_pubkey_create_serialized_batch_myunsafe(pubKeys, prvKeys, wordsRead, precPtr, precPitch);
		
		int keyLenSkip;
		int pkField = 0;
		starter = THREAD_STEPS * tIx;
		for (pkField = 0; pkField < wordsRead; pkField++, starter++) {
			keyLenSkip = 65 * pkField + 1;
			char toHash[64];
			for (int kk = keyLenSkip, len = 0; len < 64; len++) {
				toHash[len] = pubKeys[kk++];
			}
			keccak(toHash, 64, d_hash, 32);
			if (checkHashEth(d_hash)) {
				buffResult[starter] = true;
				isResult[0] = true;
			}
		}
	}
}

__global__ void workerBruteSuffix(bool* isResult, uint16_t* buffResult, size_t const prefixLen, uint64_t const startSuffix, const secp256k1_ge_storage* __restrict__ precPtr, const size_t precPitch) {
	volatile int64_t tIx = threadIdx.x + blockIdx.x * blockDim.x;
	uint32_t d_hash[8];
	const int _ITER = ITERATION[0];
	char toHash[48];
	int pkField = 0;
	unsigned char pubKeys[THREAD_STEPS_BRUTE * 65];
	unsigned char prvKeys[THREAD_STEPS_BRUTE * 32];
	uint64_t suffixNb = startSuffix;
	suffixNb += tIx * THREAD_STEPS_BRUTE;
	__shared__ __device__ char _DIC[ALPHABET_LEN];
	for (pkField = threadIdx.x; pkField < ALPHABET_LEN; pkField += blockDim.x) {
		_DIC[pkField] = _ALPHABET[pkField];
	}
	__shared__ __device__ char _PREFIX[PREFIX_MAX_LEN];
	for (pkField = threadIdx.x; pkField < PREFIX_MAX_LEN; pkField += blockDim.x) {
		_PREFIX[pkField] = PREFIX[pkField];
	}
	__syncthreads();
	pkField = 0;
	char buf[64] = { 0 };
	while (pkField < THREAD_STEPS_BRUTE) {
		int toHashIx = 0;		
		uint64_t s = suffixNb;
		int ix = lltoa(&s, buf, _DIC);
		while (ix < 64) {
			toHash[toHashIx++] = buf[ix++];
		}
		for (int ix = 0; ix < prefixLen; ix++) {
			toHash[toHashIx++] = _PREFIX[ix];
		}
		toHash[toHashIx] = 0x80;
		_GetHashSha256(toHash, toHashIx, d_hash, _ITER);
		int pk4 = (pkField << 5);
		prvKeys[pk4] = (unsigned char)__byte_perm(d_hash[0], 0, 0x3);
		prvKeys[pk4 + 1] = (unsigned char)__byte_perm(d_hash[0], 0, 0x2);
		prvKeys[pk4 + 2] = (unsigned char)__byte_perm(d_hash[0], 0, 0x1);
		prvKeys[pk4 + 3] = (unsigned char)__byte_perm(d_hash[0], 0, 0x0);
		prvKeys[pk4 + 4] = (unsigned char)__byte_perm(d_hash[1], 0, 0x3);
		prvKeys[pk4 + 5] = (unsigned char)__byte_perm(d_hash[1], 0, 0x2);
		prvKeys[pk4 + 6] = (unsigned char)__byte_perm(d_hash[1], 0, 0x1);
		prvKeys[pk4 + 7] = (unsigned char)__byte_perm(d_hash[1], 0, 0x0);
		prvKeys[pk4 + 8] = (unsigned char)__byte_perm(d_hash[2], 0, 0x3);
		prvKeys[pk4 + 9] = (unsigned char)__byte_perm(d_hash[2], 0, 0x2);
		prvKeys[pk4 + 10] = (unsigned char)__byte_perm(d_hash[2], 0, 0x1);
		prvKeys[pk4 + 11] = (unsigned char)__byte_perm(d_hash[2], 0, 0x0);
		prvKeys[pk4 + 12] = (unsigned char)__byte_perm(d_hash[3], 0, 0x3);
		prvKeys[pk4 + 13] = (unsigned char)__byte_perm(d_hash[3], 0, 0x2);
		prvKeys[pk4 + 14] = (unsigned char)__byte_perm(d_hash[3], 0, 0x1);
		prvKeys[pk4 + 15] = (unsigned char)__byte_perm(d_hash[3], 0, 0x0);
		prvKeys[pk4 + 16] = (unsigned char)__byte_perm(d_hash[4], 0, 0x3);
		prvKeys[pk4 + 17] = (unsigned char)__byte_perm(d_hash[4], 0, 0x2);
		prvKeys[pk4 + 18] = (unsigned char)__byte_perm(d_hash[4], 0, 0x1);
		prvKeys[pk4 + 19] = (unsigned char)__byte_perm(d_hash[4], 0, 0x0);
		prvKeys[pk4 + 20] = (unsigned char)__byte_perm(d_hash[5], 0, 0x3);
		prvKeys[pk4 + 21] = (unsigned char)__byte_perm(d_hash[5], 0, 0x2);
		prvKeys[pk4 + 22] = (unsigned char)__byte_perm(d_hash[5], 0, 0x1);
		prvKeys[pk4 + 23] = (unsigned char)__byte_perm(d_hash[5], 0, 0x0);
		prvKeys[pk4 + 24] = (unsigned char)__byte_perm(d_hash[6], 0, 0x3);
		prvKeys[pk4 + 25] = (unsigned char)__byte_perm(d_hash[6], 0, 0x2);
		prvKeys[pk4 + 26] = (unsigned char)__byte_perm(d_hash[6], 0, 0x1);
		prvKeys[pk4 + 27] = (unsigned char)__byte_perm(d_hash[6], 0, 0x0);
		prvKeys[pk4 + 28] = (unsigned char)__byte_perm(d_hash[7], 0, 0x3);
		prvKeys[pk4 + 29] = (unsigned char)__byte_perm(d_hash[7], 0, 0x2);
		prvKeys[pk4 + 30] = (unsigned char)__byte_perm(d_hash[7], 0, 0x1);
		prvKeys[pk4 + 31] = (unsigned char)__byte_perm(d_hash[7], 0, 0x0);
		suffixNb++;
		pkField++;
	}
	__syncwarp();
	secp256k1_ec_pubkey_create_serialized_batch_myunsafe_brute(pubKeys, prvKeys, THREAD_STEPS_BRUTE, precPtr, precPitch);
	__syncwarp();
	int keyLenSkip = 0;
	for (pkField = 0; pkField < THREAD_STEPS_BRUTE; pkField++) {
		keyLenSkip = 65 * pkField;
		_GetHash160(pubKeys, keyLenSkip, (uint8_t*)d_hash);
		if (checkHash(d_hash)) {
			buffResult[tIx] = pkField;
			isResult[0] = true;
			return;
		}
		keyLenSkip = 65 * pkField;
		_GetHash160Comp(pubKeys, keyLenSkip, (uint8_t*)d_hash);
		if (checkHash(d_hash)) {
			buffResult[tIx] = pkField;
			isResult[0] = true;
			return;
		}
		keyLenSkip = 65 * pkField;
		_GetHash160P2SHCompFromHash(d_hash, d_hash);
		if (checkHash(d_hash)) {
			buffResult[tIx] = pkField;
			isResult[0] = true;
			return;
		}
	}
}
__global__ void workerBrute(bool* isResult, uint16_t* buffResult, size_t const prefixLen, uint64_t const startSuffix, const secp256k1_ge_storage* __restrict__ precPtr, const size_t precPitch) {
	volatile int64_t tIx = threadIdx.x + blockIdx.x * blockDim.x;
	uint32_t d_hash[8];
	const int _ITER = ITERATION[0];
	char toHash[48];
	int pkField = 0;
	unsigned char pubKeys[THREAD_STEPS_BRUTE * 65];
	unsigned char prvKeys[THREAD_STEPS_BRUTE * 32];
	uint64_t suffixNb = startSuffix;
	suffixNb += tIx * THREAD_STEPS_BRUTE;
	for (int ix = 0; ix < prefixLen; ix++) {
		toHash[ix] = PREFIX[ix];
	}

	__shared__ __device__ uint32_t _SHA_PRE[8];
	for (int i = threadIdx.x; i < 8; i+= blockDim.x) {
		_SHA_PRE[i] = SHA_PRE[i];
	}
	__shared__ __device__ char _DIC[ALPHABET_LEN];
	for (pkField = threadIdx.x; pkField < ALPHABET_LEN; pkField += blockDim.x) {
		_DIC[pkField] = _ALPHABET[pkField];
	}
	__shared__ __device__ int _SHA_LEVE_SHARED[1];
	if (threadIdx.x == 0) {
		_SHA_LEVE_SHARED[0] = SHA_LEVEL[0];
	}

	__syncthreads();
	pkField = 0;
	while (pkField < THREAD_STEPS_BRUTE) {
		_GetHashSha256FromPre(toHash, _SHA_LEVE_SHARED[0], suffixNb, prefixLen, d_hash, _SHA_PRE, _DIC, _ITER);
		int pk4 = (pkField << 5);
		prvKeys[pk4] = (unsigned char)__byte_perm(d_hash[0], 0, 0x3);
		prvKeys[pk4 + 1] = (unsigned char)__byte_perm(d_hash[0], 0, 0x2);
		prvKeys[pk4 + 2] = (unsigned char)__byte_perm(d_hash[0], 0, 0x1);
		prvKeys[pk4 + 3] = (unsigned char)__byte_perm(d_hash[0], 0, 0x0);
		prvKeys[pk4 + 4] = (unsigned char)__byte_perm(d_hash[1], 0, 0x3);
		prvKeys[pk4 + 5] = (unsigned char)__byte_perm(d_hash[1], 0, 0x2);
		prvKeys[pk4 + 6] = (unsigned char)__byte_perm(d_hash[1], 0, 0x1);
		prvKeys[pk4 + 7] = (unsigned char)__byte_perm(d_hash[1], 0, 0x0);
		prvKeys[pk4 + 8] = (unsigned char)__byte_perm(d_hash[2], 0, 0x3);
		prvKeys[pk4 + 9] = (unsigned char)__byte_perm(d_hash[2], 0, 0x2);
		prvKeys[pk4 + 10] = (unsigned char)__byte_perm(d_hash[2], 0, 0x1);
		prvKeys[pk4 + 11] = (unsigned char)__byte_perm(d_hash[2], 0, 0x0);
		prvKeys[pk4 + 12] = (unsigned char)__byte_perm(d_hash[3], 0, 0x3);
		prvKeys[pk4 + 13] = (unsigned char)__byte_perm(d_hash[3], 0, 0x2);
		prvKeys[pk4 + 14] = (unsigned char)__byte_perm(d_hash[3], 0, 0x1);
		prvKeys[pk4 + 15] = (unsigned char)__byte_perm(d_hash[3], 0, 0x0);
		prvKeys[pk4 + 16] = (unsigned char)__byte_perm(d_hash[4], 0, 0x3);
		prvKeys[pk4 + 17] = (unsigned char)__byte_perm(d_hash[4], 0, 0x2);
		prvKeys[pk4 + 18] = (unsigned char)__byte_perm(d_hash[4], 0, 0x1);
		prvKeys[pk4 + 19] = (unsigned char)__byte_perm(d_hash[4], 0, 0x0);
		prvKeys[pk4 + 20] = (unsigned char)__byte_perm(d_hash[5], 0, 0x3);
		prvKeys[pk4 + 21] = (unsigned char)__byte_perm(d_hash[5], 0, 0x2);
		prvKeys[pk4 + 22] = (unsigned char)__byte_perm(d_hash[5], 0, 0x1);
		prvKeys[pk4 + 23] = (unsigned char)__byte_perm(d_hash[5], 0, 0x0);
		prvKeys[pk4 + 24] = (unsigned char)__byte_perm(d_hash[6], 0, 0x3);
		prvKeys[pk4 + 25] = (unsigned char)__byte_perm(d_hash[6], 0, 0x2);
		prvKeys[pk4 + 26] = (unsigned char)__byte_perm(d_hash[6], 0, 0x1);
		prvKeys[pk4 + 27] = (unsigned char)__byte_perm(d_hash[6], 0, 0x0);
		prvKeys[pk4 + 28] = (unsigned char)__byte_perm(d_hash[7], 0, 0x3);
		prvKeys[pk4 + 29] = (unsigned char)__byte_perm(d_hash[7], 0, 0x2);
		prvKeys[pk4 + 30] = (unsigned char)__byte_perm(d_hash[7], 0, 0x1);
		prvKeys[pk4 + 31] = (unsigned char)__byte_perm(d_hash[7], 0, 0x0);
		//printf("  privkey %d\n", prvKeys[0]);
		suffixNb++;
		pkField++;
	}
	__syncwarp();
	secp256k1_ec_pubkey_create_serialized_batch_myunsafe_brute(pubKeys, prvKeys, THREAD_STEPS_BRUTE, precPtr, precPitch);
	__syncwarp();
	int keyLenSkip = 0;
	for (pkField = 0; pkField < THREAD_STEPS_BRUTE; pkField++) {
		keyLenSkip = 65 * pkField;
		_GetHash160(pubKeys, keyLenSkip, (uint8_t*)d_hash);
		if (checkHash(d_hash)) {
			buffResult[tIx] = pkField;
			isResult[0] = true;
			return;
		}
		keyLenSkip = 65 * pkField;
		_GetHash160Comp(pubKeys, keyLenSkip, (uint8_t*)d_hash);
		if (checkHash(d_hash)) {
			buffResult[tIx] = pkField;
			isResult[0] = true;
			return;
		}
		keyLenSkip = 65 * pkField;
		_GetHash160P2SHCompFromHash(d_hash, d_hash);
		if (checkHash(d_hash)) {
			buffResult[tIx] = pkField;
			isResult[0] = true;
			return;
		}
	}
}


__global__ void workerBruteEthSuffix(bool* isResult, uint16_t* buffResult, size_t const prefixLen, uint64_t const startSuffix, const short workType, const bool sha256hash, const secp256k1_ge_storage* __restrict__ precPtr, const size_t precPitch) {
	volatile int64_t tIx = threadIdx.x + blockIdx.x * blockDim.x;
	unsigned char d_hash[32];
	const int _ITER = ITERATION[0];
	char toHash[64];
	int pkField = 0;
	unsigned char pubKeys[THREAD_STEPS_BRUTE * 65];
	unsigned char prvKeys[THREAD_STEPS_BRUTE * 32];
	uint64_t suffixNb = startSuffix;
	suffixNb += tIx * THREAD_STEPS_BRUTE;

	__shared__ __device__ char _DIC[ALPHABET_LEN];
	for (pkField = threadIdx.x; pkField < ALPHABET_LEN; pkField += blockDim.x) {
		_DIC[pkField] = _ALPHABET[pkField];
	}
	__shared__ __device__ char _PREFIX[PREFIX_MAX_LEN];
	for (pkField = threadIdx.x; pkField < PREFIX_MAX_LEN; pkField += blockDim.x) {
		_PREFIX[pkField] = PREFIX[pkField];
	}
	uint32_t d_hashsha[8];
	pkField = 0;
	__syncthreads();
	char buf[64] = { 0 };
	for (short turn = (workType % 2 == 0) ? 2 : 1; turn <= workType; turn *= 2) {
		while (pkField < THREAD_STEPS_BRUTE) {
			uint64_t s = suffixNb;
			int ix = lltoa(&s, buf, _DIC);
			int toHashIx = 0;
			while (ix < 64) {
				toHash[toHashIx++] = buf[ix++];
			}
			for (int ix = 0; ix < prefixLen; ix++) {
				toHash[toHashIx++] = _PREFIX[ix];
			}
			if (sha256hash) {											
				toHash[toHashIx] = 0x80;
				_GetHashSha256(toHash, toHashIx, d_hashsha, _ITER);
				int pk4 = pkField << 5;
				prvKeys[pk4] = (unsigned char)__byte_perm(d_hashsha[0], 0, 0x3);
				prvKeys[pk4 + 1] = (unsigned char)__byte_perm(d_hashsha[0], 0, 0x2);
				prvKeys[pk4 + 2] = (unsigned char)__byte_perm(d_hashsha[0], 0, 0x1);
				prvKeys[pk4 + 3] = (unsigned char)__byte_perm(d_hashsha[0], 0, 0x0);
				prvKeys[pk4 + 4] = (unsigned char)__byte_perm(d_hashsha[1], 0, 0x3);
				prvKeys[pk4 + 5] = (unsigned char)__byte_perm(d_hashsha[1], 0, 0x2);
				prvKeys[pk4 + 6] = (unsigned char)__byte_perm(d_hashsha[1], 0, 0x1);
				prvKeys[pk4 + 7] = (unsigned char)__byte_perm(d_hashsha[1], 0, 0x0);
				prvKeys[pk4 + 8] = (unsigned char)__byte_perm(d_hashsha[2], 0, 0x3);
				prvKeys[pk4 + 9] = (unsigned char)__byte_perm(d_hashsha[2], 0, 0x2);
				prvKeys[pk4 + 10] = (unsigned char)__byte_perm(d_hashsha[2], 0, 0x1);
				prvKeys[pk4 + 11] = (unsigned char)__byte_perm(d_hashsha[2], 0, 0x0);
				prvKeys[pk4 + 12] = (unsigned char)__byte_perm(d_hashsha[3], 0, 0x3);
				prvKeys[pk4 + 13] = (unsigned char)__byte_perm(d_hashsha[3], 0, 0x2);
				prvKeys[pk4 + 14] = (unsigned char)__byte_perm(d_hashsha[3], 0, 0x1);
				prvKeys[pk4 + 15] = (unsigned char)__byte_perm(d_hashsha[3], 0, 0x0);
				prvKeys[pk4 + 16] = (unsigned char)__byte_perm(d_hashsha[4], 0, 0x3);
				prvKeys[pk4 + 17] = (unsigned char)__byte_perm(d_hashsha[4], 0, 0x2);
				prvKeys[pk4 + 18] = (unsigned char)__byte_perm(d_hashsha[4], 0, 0x1);
				prvKeys[pk4 + 19] = (unsigned char)__byte_perm(d_hashsha[4], 0, 0x0);
				prvKeys[pk4 + 20] = (unsigned char)__byte_perm(d_hashsha[5], 0, 0x3);
				prvKeys[pk4 + 21] = (unsigned char)__byte_perm(d_hashsha[5], 0, 0x2);
				prvKeys[pk4 + 22] = (unsigned char)__byte_perm(d_hashsha[5], 0, 0x1);
				prvKeys[pk4 + 23] = (unsigned char)__byte_perm(d_hashsha[5], 0, 0x0);
				prvKeys[pk4 + 24] = (unsigned char)__byte_perm(d_hashsha[6], 0, 0x3);
				prvKeys[pk4 + 25] = (unsigned char)__byte_perm(d_hashsha[6], 0, 0x2);
				prvKeys[pk4 + 26] = (unsigned char)__byte_perm(d_hashsha[6], 0, 0x1);
				prvKeys[pk4 + 27] = (unsigned char)__byte_perm(d_hashsha[6], 0, 0x0);
				prvKeys[pk4 + 28] = (unsigned char)__byte_perm(d_hashsha[7], 0, 0x3);
				prvKeys[pk4 + 29] = (unsigned char)__byte_perm(d_hashsha[7], 0, 0x2);
				prvKeys[pk4 + 30] = (unsigned char)__byte_perm(d_hashsha[7], 0, 0x1);
				prvKeys[pk4 + 31] = (unsigned char)__byte_perm(d_hashsha[7], 0, 0x0);
			}
			else {
				keccak(toHash, toHashIx, d_hash, 32);
				if (turn == 2) {
					char toHash2[48];
#pragma unroll 2030
					for (int camp2Turns = 0; camp2Turns < 2030; camp2Turns++) {
#pragma unroll 32
						for (int xxx = 0; xxx < 32; xxx++) {
							toHash2[xxx] = d_hash[xxx];
						}
						keccak(toHash2, 32, d_hash, 32);
					}
				}
#pragma unroll 32
				for (int ixTemp = 32 * pkField, len = 0; len < 32; len++) {
					prvKeys[ixTemp++] = d_hash[len];
				}
			}
			pkField++;
			suffixNb++;
		}

		__syncwarp();
		secp256k1_ec_pubkey_create_serialized_batch_myunsafe_brute(pubKeys, prvKeys, THREAD_STEPS_BRUTE, precPtr, precPitch);
		__syncwarp();
		for (int pk = 0; pk < THREAD_STEPS_BRUTE; pk++) {
			for (int kk = 65 * pk + 1, len = 0; len < 64; len++) {
				toHash[len] = pubKeys[kk++];
			}
			keccak(toHash, 64, d_hash, 32);
			if (checkHashEth(d_hash)) {
				buffResult[tIx] = pk;
				isResult[0] = true;
				return;
			}
		}
	}
}

__global__ void workerBruteEth(bool* isResult, uint16_t* buffResult, size_t const prefixLen, uint64_t const startSuffix, const short workType, const bool sha256hash, const secp256k1_ge_storage* __restrict__ precPtr, const size_t precPitch) {
	volatile int64_t tIx = threadIdx.x + blockIdx.x * blockDim.x;
	unsigned char d_hash[32];
	const int _ITER = ITERATION[0];
	char toHash[64];
	int pkField = 0;
	unsigned char pubKeys[THREAD_STEPS_BRUTE * 65];
	unsigned char prvKeys[THREAD_STEPS_BRUTE * 32];
	uint64_t suffixNb = startSuffix;
	suffixNb += tIx * THREAD_STEPS_BRUTE;
	for (int ix = 0; ix < prefixLen; ix++) {
		toHash[ix] = PREFIX[ix];
	}

	__shared__ __device__ uint32_t _SHA_PRE[8];
	for (int i = threadIdx.x; i < 8; i += blockDim.x) {
		_SHA_PRE[i] = SHA_PRE[i];
	}
	__shared__ __device__ char _DIC[ALPHABET_LEN];
	for (pkField = threadIdx.x; pkField < ALPHABET_LEN; pkField += blockDim.x) {
		_DIC[pkField] = _ALPHABET[pkField];
	}
	__shared__ __device__ int _SHA_LEVE_SHARED[1];
	if (threadIdx.x == 0) {
		_SHA_LEVE_SHARED[0] = SHA_LEVEL[0];
	}
	uint32_t d_hashsha[8];
	pkField = 0;
	__syncthreads();
		
	for (short turn = (workType % 2 == 0) ? 2 : 1; turn <= workType; turn *= 2) {
		while (pkField < THREAD_STEPS_BRUTE) {
			if (sha256hash) {
				_GetHashSha256FromPre(toHash, _SHA_LEVE_SHARED[0], suffixNb, prefixLen, d_hashsha, _SHA_PRE, _DIC, _ITER);
				int pk4 = pkField << 5;
				prvKeys[pk4] = (unsigned char)__byte_perm(d_hashsha[0], 0, 0x3);
				prvKeys[pk4 + 1] = (unsigned char)__byte_perm(d_hashsha[0], 0, 0x2);
				prvKeys[pk4 + 2] = (unsigned char)__byte_perm(d_hashsha[0], 0, 0x1);
				prvKeys[pk4 + 3] = (unsigned char)__byte_perm(d_hashsha[0], 0, 0x0);
				prvKeys[pk4 + 4] = (unsigned char)__byte_perm(d_hashsha[1], 0, 0x3);
				prvKeys[pk4 + 5] = (unsigned char)__byte_perm(d_hashsha[1], 0, 0x2);
				prvKeys[pk4 + 6] = (unsigned char)__byte_perm(d_hashsha[1], 0, 0x1);
				prvKeys[pk4 + 7] = (unsigned char)__byte_perm(d_hashsha[1], 0, 0x0);
				prvKeys[pk4 + 8] = (unsigned char)__byte_perm(d_hashsha[2], 0, 0x3);
				prvKeys[pk4 + 9] = (unsigned char)__byte_perm(d_hashsha[2], 0, 0x2);
				prvKeys[pk4 + 10] = (unsigned char)__byte_perm(d_hashsha[2], 0, 0x1);
				prvKeys[pk4 + 11] = (unsigned char)__byte_perm(d_hashsha[2], 0, 0x0);
				prvKeys[pk4 + 12] = (unsigned char)__byte_perm(d_hashsha[3], 0, 0x3);
				prvKeys[pk4 + 13] = (unsigned char)__byte_perm(d_hashsha[3], 0, 0x2);
				prvKeys[pk4 + 14] = (unsigned char)__byte_perm(d_hashsha[3], 0, 0x1);
				prvKeys[pk4 + 15] = (unsigned char)__byte_perm(d_hashsha[3], 0, 0x0);
				prvKeys[pk4 + 16] = (unsigned char)__byte_perm(d_hashsha[4], 0, 0x3);
				prvKeys[pk4 + 17] = (unsigned char)__byte_perm(d_hashsha[4], 0, 0x2);
				prvKeys[pk4 + 18] = (unsigned char)__byte_perm(d_hashsha[4], 0, 0x1);
				prvKeys[pk4 + 19] = (unsigned char)__byte_perm(d_hashsha[4], 0, 0x0);
				prvKeys[pk4 + 20] = (unsigned char)__byte_perm(d_hashsha[5], 0, 0x3);
				prvKeys[pk4 + 21] = (unsigned char)__byte_perm(d_hashsha[5], 0, 0x2);
				prvKeys[pk4 + 22] = (unsigned char)__byte_perm(d_hashsha[5], 0, 0x1);
				prvKeys[pk4 + 23] = (unsigned char)__byte_perm(d_hashsha[5], 0, 0x0);
				prvKeys[pk4 + 24] = (unsigned char)__byte_perm(d_hashsha[6], 0, 0x3);
				prvKeys[pk4 + 25] = (unsigned char)__byte_perm(d_hashsha[6], 0, 0x2);
				prvKeys[pk4 + 26] = (unsigned char)__byte_perm(d_hashsha[6], 0, 0x1);
				prvKeys[pk4 + 27] = (unsigned char)__byte_perm(d_hashsha[6], 0, 0x0);
				prvKeys[pk4 + 28] = (unsigned char)__byte_perm(d_hashsha[7], 0, 0x3);
				prvKeys[pk4 + 29] = (unsigned char)__byte_perm(d_hashsha[7], 0, 0x2);
				prvKeys[pk4 + 30] = (unsigned char)__byte_perm(d_hashsha[7], 0, 0x1);
				prvKeys[pk4 + 31] = (unsigned char)__byte_perm(d_hashsha[7], 0, 0x0);
			}
			else {
				int toHashIx = prefixLen;
				char buf[64] = { 0 };
				uint64_t suffixNbPass = suffixNb;
				int ix = lltoa(&suffixNbPass, buf, _DIC);
				while (ix < 64) {
					toHash[toHashIx++] = buf[ix++];
				}
				keccak(toHash, toHashIx, d_hash, 32);
				if (turn == 2) {
					char toHash2[48];
#pragma unroll 2030
					for (int camp2Turns = 0; camp2Turns < 2030; camp2Turns++) {
#pragma unroll 32
						for (int xxx = 0; xxx < 32; xxx++) {
							toHash2[xxx] = d_hash[xxx];
						}
						keccak(toHash2, 32, d_hash, 32);
					}
				}
#pragma unroll 32
				for (int ixTemp = 32 * pkField, len = 0; len < 32; len++) {
					prvKeys[ixTemp++] = d_hash[len];
				}
			}
			pkField++;
			suffixNb++;
		}

		__syncwarp();
		secp256k1_ec_pubkey_create_serialized_batch_myunsafe_brute(pubKeys, prvKeys, THREAD_STEPS_BRUTE, precPtr, precPitch);
		__syncwarp();
		for (int pk = 0; pk < THREAD_STEPS_BRUTE; pk++) {
			for (int kk = 65 * pk + 1, len = 0; len < 64; len++) {
				toHash[len] = pubKeys[kk++];
			}
			keccak(toHash, 64, d_hash, 32);
			if (checkHashEth(d_hash)) {
				buffResult[tIx] = pk;
				isResult[0] = true;
				return;
			}
		}
	}
}

__global__ void ecmult_big_create(secp256k1_gej * gej_temp, secp256k1_fe * z_ratio, secp256k1_ge_storage* precPtr, size_t precPitch, unsigned int bits) {
	int64_t tIx = threadIdx.x + blockIdx.x * blockDim.x;
	if (tIx != 0) {
		return;
	}
	unsigned int windows;
	size_t window_size;
	size_t i, row;
	secp256k1_fe  fe_zinv;
	secp256k1_ge  ge_temp;
	secp256k1_ge  ge_window_one = secp256k1_ge_const_g;
	secp256k1_gej gej_window_base;

	/* We +1 to account for a possible high 1 bit after converting the privkey to signed digit form.    */
	/* This means our table reaches to 257 bits even though the privkey scalar is at most 256 bits.     */
	//unsigned int bits = (unsigned int)ECMULT_WINDOW_SIZE;
	windows = (256 / bits) + 1;
	window_size = (1 << (bits - 1));
	WINDOWS = windows;
	WINDOW_SIZE = WINDOW_SIZE;
	ECMULT_WINDOW_SIZE = bits;
	//size_t total_size = (256 / bits) * window_size + (1 << (256 % bits));

	//windows = WINDOWS;
	//window_size = WINDOW_SIZE;
	//bits = ECMULT_WINDOW_SIZE;

	/* Total number of required point storage elements.                                 */
	/* This differs from the (windows * window_size) because the last row can be shrunk */
	/*   as it only needs to extend enough to include a possible 1 in the 257th bit.    */
	//total_size = (256 / bits) * window_size + (1 << (256 % bits));

	//rtn->gej_temp = (secp256k1_gej*)checked_malloc(&ctx->error_callback, sizeof(secp256k1_gej) * window_size);
	//rtn->z_ratio = (secp256k1_fe*)checked_malloc(&ctx->error_callback, sizeof(secp256k1_fe) * window_size);

	/************ Precomputed Table Initialization ************/
	secp256k1_gej_set_ge(&gej_window_base, &ge_window_one);

	/* This is the same for all windows.    */
	secp256k1_fe_set_int(&(z_ratio[0]), 0);


	for (row = 0; row < windows; row++) {
		/* The last row is a bit smaller, only extending to include the 257th bit. */
		window_size = (row == windows - 1 ? (1 << (256 % bits)) : (1 << (bits - 1)));

		/* The base element of each row is 2^bits times the previous row's base. */
		if (row > 0) {
			for (i = 0; i < bits; i++) {
				secp256k1_gej_double_var(&gej_window_base, &gej_window_base, NULL); 
			}
		}
		gej_temp[0] = gej_window_base;

		/* The base element is also our "one" value for this row.   */
		/* If we are at offset 2^X, adding "one" should add 2^X.    */
		secp256k1_ge_set_gej(&ge_window_one, &gej_window_base);


		/* Repeated + 1s to fill the rest of the row.   */

		/* We capture the Z ratios between consecutive points for quick Z inversion.    */
		/*   gej_temp[i-1].z * z_ratio[i] => gej_temp[i].z                              */
		/* This means that z_ratio[i] = (gej_temp[i-1].z)^-1 * gej_temp[i].z            */
		/* If we know gej_temp[i].z^-1, we can get gej_temp[i-1].z^1 using z_ratio[i]   */
		/* Visually:                                    */
		/* i            0           1           2       */
		/* gej_temp     a           b           c       */
		/* z_ratio     NaN      (a^-1)*b    (b^-1)*c    */
		for (i = 1; i < window_size; i++) {
			secp256k1_gej_add_ge_var(&(gej_temp[i]), &(gej_temp[i - 1]), &ge_window_one, &(z_ratio[i]));
		}


		/* An unpacked version of secp256k1_ge_set_table_gej_var() that works   */
		/*   element by element instead of requiring a secp256k1_ge *buffer.    */	

		/* Invert the last Z coordinate manually.   */
		i = window_size - 1;
		secp256k1_fe_inv(&fe_zinv, &(gej_temp[i].z));
		secp256k1_ge_set_gej_zinv(&ge_temp, &(gej_temp[i]), &fe_zinv);
		//secp256k1_ge_to_storage(&(prec[row][i]), &ge_temp);
		secp256k1_ge_storage* ROW_PREC = (secp256k1_ge_storage*)((char*)precPtr + row * precPitch) + i;
		secp256k1_ge_to_storage(ROW_PREC, &ge_temp);

		/* Use the last element's known Z inverse to determine the previous' Z inverse. */
		for (; i > 0; i--) {
			/* fe_zinv = (gej_temp[i].z)^-1                 */
			/* (gej_temp[i-1].z)^-1 = z_ratio[i] * fe_zinv  */
			secp256k1_fe_mul(&fe_zinv, &fe_zinv, &(z_ratio[i]));
			/* fe_zinv = (gej_temp[i-1].z)^-1               */

			secp256k1_ge_set_gej_zinv(&ge_temp, &(gej_temp[i - 1]), &fe_zinv);
			//secp256k1_ge_to_storage(&(prec[row][i - 1]), &ge_temp);

			secp256k1_ge_storage* ROW_PRECi_1 = (secp256k1_ge_storage*)((char*)precPtr + row * precPitch) + (i-1);
			secp256k1_ge_to_storage(ROW_PRECi_1, &ge_temp);
		}
	}
}

__global__ void shaPre(size_t const prefixLen, uint32_t* hashPre) {
	char toHash[48];
	uint32_t d_hash[8];
	for (int i = 0; i < prefixLen; i++) {
		toHash[i] = PREFIX[i];
	}
	_GetHashSha256Pre(toHash, SHA_LEVEL[0], prefixLen, d_hash, hashPre);
}

__device__ bool checkHashEth(unsigned char d_hash[32]) {
	uint32_t hash[5];
	for (int h = 12, i = 0; i < 5; i++) {
		hash[i] = (d_hash[h++]) | ((d_hash[h++] << 8) & 0x0000ff00) | ((d_hash[h++] << 16) & 0x00ff0000) | ((d_hash[h++] << 24) & 0xff000000) ;
	}
	return checkHash(hash);
}

__device__ bool checkHash(const uint32_t hash[5]) {
	uint64_t mask = _BLOOM_FILTER_MASK64[0];
	uint32_t* bloomFilter = _BLOOM_FILTER[0];
	uint64_t idx[5];

	idx[0] = ((uint64_t)hash[0] << 32 | hash[1]) & mask;
	idx[1] = ((uint64_t)hash[2] << 32 | hash[3]) & mask;
	idx[2] = ((uint64_t)(hash[0] ^ hash[1]) << 32 | (hash[1] ^ hash[2])) & mask;
	idx[3] = ((uint64_t)(hash[2] ^ hash[3]) << 32 | (hash[3] ^ hash[4])) & mask;
	idx[4] = ((uint64_t)(hash[0] ^ hash[3]) << 32 | (hash[1] ^ hash[3])) & mask;

	for (int i = 0; i < 5; i++) {
		uint32_t f = bloomFilter[idx[i] / 32];

		if ((f & (0x01 << (idx[i] % 32))) == 0) {
			return false;
		}
	}

	return true;
}

cudaError_t loadPrefix(const char* _prefix, size_t const prefixLen) {
	return cudaMemcpyToSymbol(PREFIX, _prefix, prefixLen * sizeof(char));
}
cudaError_t loadIteration(int _i) {
	int _l[1];
	_l[0] = _i;
	return cudaMemcpyToSymbol(ITERATION, _l, 1 * sizeof(int));
}
cudaError_t loadLevel(int _level) {
	int _l[1];
	_l[0] = _level;
	return cudaMemcpyToSymbol(SHA_LEVEL, _l, 1 * sizeof(int));
}
cudaError_t loadShaPre(uint32_t* _pre) {
	return cudaMemcpyToSymbol(SHA_PRE, _pre, 8 * sizeof(uint32_t));
}
cudaError_t cudaMemcpyToSymbol_BLOOM_FILTER(uint32_t* _bloomFilterPtr) {
	return cudaMemcpyToSymbol(_BLOOM_FILTER, &_bloomFilterPtr, sizeof(uint32_t*));
}
cudaError_t cudaMemcpyToSymbol_BLOOM_FILTER_MASK64(uint64_t* bloomFilterMask) {
	return cudaMemcpyToSymbol(_BLOOM_FILTER_MASK64, bloomFilterMask, sizeof(uint64_t*));
}
cudaError_t cudaMemcpyToSymbol_USE_BLOOM_FILTER(uint32_t* useBloomFilter) {
	return cudaMemcpyToSymbol(_USE_BLOOM_FILTER, useBloomFilter, sizeof(uint32_t));
}

cudaError_t loadWindow(unsigned int windowSize, unsigned int windows) {
	int _l[1];
	_l[0] = windows;
	cudaMemcpyToSymbol(WINDOWS_SIZE_CONST, _l, 1 * sizeof(unsigned int));
	_l[0] = windowSize;
	return cudaMemcpyToSymbol(ECMULT_WINDOW_SIZE_CONST, _l, 1 * sizeof(unsigned int));
}

__device__ __inline__ void _GetHashSha256(const char* __restrict__ toHash, size_t const prefixLen, uint32_t* d_hash, const int count) {
	_GetHashSha256(toHash, prefixLen, d_hash);
	for (int c = 1; c < count; c++) {
		uint32_t data[16] = { 0 };
#pragma unroll 8
		for (int i = 0; i < 8; i++) {
			data[i] = d_hash[i];
		}
		data[8] = 0x80000000;
		data[15] = 0x100;
		SHA256Initialize(d_hash);
		SHA256Transform(d_hash, data);
	}
}

__device__ __inline__ void _GetHashSha256(const char* __restrict__ toHash, size_t const prefixLen, uint32_t* d_hash) {
	int ix = 0;
	uint32_t data[16];
	data[0] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
	data[1] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
	data[2] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
	data[3] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
	data[4] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
	data[5] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
	data[6] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
	data[7] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
	data[8] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
	data[9] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
	data[10] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
	data[11] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
	data[12] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
	data[13] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
	data[14] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
	data[15] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | (prefixLen * 8);
	SHA256Initialize(d_hash);
	SHA256Transform(d_hash, data);
}
__device__ __inline__ void _GetHashSha256Pre(char* toHash, int level, size_t const prefixLen, uint32_t* d_hash, uint32_t* hashPre) {
	int ix = 0;
	uint32_t data[16];
	//int level = 0;
	if (level > 0) {
		if (prefixLen >= 4) {
			data[0] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
			if (prefixLen >= 8) {
				data[1] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
				if (prefixLen >= 12) {
					data[2] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
					if (prefixLen >= 16) {
						data[3] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
						if (prefixLen >= 20) {
							data[4] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
						}
					}
				}
			}
		}
		SHA256Initialize(d_hash);
		SHA256TransformPre(level, d_hash, data, hashPre);
	}
}

__device__ __noinline__ void _GetHashSha256FromPre(char* toHash, int& level, uint64_t suffix, size_t const prefixLen, uint32_t* d_hash, uint32_t* hashPre, const char* dict, const int count) {
	_GetHashSha256FromPre(toHash, level, suffix, prefixLen, d_hash, hashPre, dict);
	for (int c = 1; c < count; c++) {
		uint32_t data[16] = { 0 };
#pragma unroll 8
		for (int i = 0; i < 8; i++) {
			data[i] = d_hash[i];
		}
		data[8] = 0x80000000;
		data[15] = 0x100;
		SHA256Initialize(d_hash);
		SHA256Transform(d_hash, data);
	}
}

__device__ __noinline__ void _GetHashSha256FromPre(char* toHash, int& level, uint64_t suffix, size_t const prefixLen, uint32_t* d_hash, uint32_t* hashPre, const char * dict) {
	int toHashIx = prefixLen;
	char buf[64] = { 0 };
	int ix = lltoa(&suffix, buf, dict);
	while (ix < 64) {
		toHash[toHashIx++] = buf[ix++];
	}
	toHash[toHashIx] = 0x80;
	ix = 0;
	uint32_t data[16];

#pragma unroll 12
	for (int z = 0; z < 12; z++) {
		data[z] = ((toHash[ix++] << 24) & 0xff000000) | ((toHash[ix++] << 16) & 0x00ff0000) | ((toHash[ix++] << 8) & 0x0000ff00) | ((toHash[ix++]) & 0x000000ff);
	}

	data[12] = 0x00000000;
	data[13] = 0x00000000;
	data[14] = 0x00000000;
	data[15] = toHashIx * 8;
	SHA256Initialize(d_hash);
	SHA256TransformFromPre(level, d_hash, data, hashPre);
}

__device__ __noinline__ int lltoa(uint64_t * __restrict__ val, char* __restrict__ buf, const char * __restrict__ dict) {
	int i = 63;
	for (; *val && i; --i, *val /= ALPHABET_LEN) {
		buf[i] = dict[*val % ALPHABET_LEN];
	}
	return i + 1;
}