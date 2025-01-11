#ifndef _HASH_LOOKUP_HOST_H
#define _HASH_LOOKUP_HOST_H

typedef struct hash160 {
	uint32_t h[5];

	hash160(const uint32_t hash[5])
	{
		memcpy(h, hash, sizeof(uint32_t) * 5);
	}
}hash160;

class CudaHashLookup {

private:
	uint32_t *_bloomFilterPtr;

	cudaError_t setTargetBloomFilter(const std::vector<struct hash160> &targets, uint32_t bbits);
	
	uint32_t getOptimalBloomFilterBits(double p, size_t n);

	void cleanup();
	
	void initializeBloomFilter64(const std::vector<struct hash160> &targets, uint32_t*filter, uint64_t mask);

public:

	CudaHashLookup()
	{
		_bloomFilterPtr = NULL;
	}

	~CudaHashLookup()
	{
		cleanup();
	}

	cudaError_t setTargets(const std::vector<struct hash160> &targets, uint32_t bbits);
};

#endif