#include <iostream>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <chrono>
#include <set>

#include "Kernel.cuh"
#include "lib/util.h"
#include "lib/Bech32.h"
#include "lib/V/VBase58.h"
#include "lib/hash/sha256.h"
#include "lib/hash/sha3.cuh"
#include "lib/Int.h"
#include "lib/SECP256k1.h"
#include "lib/CudaHashLookup.h"

using namespace std;

cudaError_t prepareCuda();
cudaError_t processCuda2(std::istream& in);
cudaError_t processCuda2Eth(std::istream& stream);

cudaError_t processCudaBrute();
cudaError_t processCudaBruteEth();

bool readArgs(int argc, char** argv);
bool checkDevice();
bool processCandidate(Int& toTest, string text);
bool processCandidateEth(Int& toTest, string text, bool camp2);
bool readFileAddress(const std::string& file_name);
void printSpeed(double speed);
void saveStatus();
void saveStatusBrute();
void checkResult(bool* buffDeviceResult, uint64_t outputSizeB_T, string combined, vector<uint32_t> indexes);
void checkResultEth(bool* buffDeviceResult, uint64_t outputSizeB_T, string combined, vector<uint32_t> indexes);
int _lltoa(uint64_t val, char* buf);

char* __strlwr(char* str);

const uint16_t NO_RESULT = UINT16_MAX - 1;

int DEVICE_NR = 0;
unsigned int BLOCK_THREADS = 0;
unsigned int BLOCK_NUMBER = 0;
uint64_t workSize;
uint64_t counterTotal = 0;
const int fileStatusInterval = 60;
bool IS_VERBOSE = false;

string fileResult = "result.txt";
string fileStatus = "fileStatus.txt";
string fileAddress = "";
string filePhrase = "";
bool PHRASE_IN = false;

bool IS_ETH = false;
bool IS_CAMP2 = false;
bool IS_ETH_SHA256 = false;

bool IS_BRUTE = false;
string ROOT = "";
bool IS_ROOT = false;
string ROOT_SUFFIX = "";
bool IS_ROOT_SUFFIX = false;
size_t ROOT_LEN = 0;
uint64_t SUFFIX = 0;

bool IS_RESULT = false;
string STATUS = "";

bool opt_show_false = false;
uint32_t opt_bloom_bits = 0;

int HASH_ITERATION = 1;

#ifdef _DEBUG
unsigned int PARAM_ECMULT_WINDOW_SIZE = 4;
#else
unsigned int PARAM_ECMULT_WINDOW_SIZE = 18;
#endif

union hash160_20 {
    uint32_t bits5[5];
    uint8_t bits20[20];
};

size_t hashCount = 0;
set<string>addresses;
//set<string>pubkeys;
CudaHashLookup _targetLookup;
std::vector<hash160> _targets;

Secp256K1* secp;

secp256k1_ge_storage* _dev_precomp;
size_t pitch;

int main(int argc, char** argv)
{
    printf("BrainWords Btc+Eth v 0.8\n\n");

    if (argc == 1 || !readArgs(argc, argv)) {
        printf("Usage:\n");
        printf("--bbits BITS\n\tBloom filter optimal bits to override (12-36).\n");
        printf("--bits ECMULT_WINDOW_SIZE\n\tNumber of bits to generate ecmult table (e.g. 4, 8, 16, 20, 24).\n");
        printf("--camp2\n\tUse camp2 calculation.\n");
        printf("--eth\n\tEnable Ethereum mode.\n");
        printf("--ethsha256\n\tUse eth's sha256 calculation.\n");
        printf("--fstatus FILE_STATUS\n\tFile status.\n");
        printf("--inputAddress FILE\n\tFile with input addresses (in hex).\n");
        printf("--inputIn\n\tRead phrases from the standard input.\n");
        printf("--inputPhrase FILE\n\tFile with list of phrases to check.\n");
        printf("--iteration NUM\n\tNumber of hash iterations to perform.\n");
        printf("--root PHRASE\n\tUse the root phrase for the combinations (ROOT + combinations).\n");
        printf("--rootsuffix PHRASE\n\tUse the root phrase for the combinations (combinations + ROOTSUFFIX).\n");
        printf("--suffix RESUME_NO\n\tSpecify number of combinations to resume from.\n");
        printf("-b BLOCK_NUMBER\n\tNumber of block.\n");
        printf("-d DEVICE_NR\n\tDevice number to use.\n");
        printf("-fp\n\tShow false positive results.\n");
        printf("-o FILE\n\tResult file (default: result.txt).\n");
        printf("-os FILE\n\tStatus file (default: fileStatus.txt).\n");
        printf("-t BLOCK_THREADS\n\tNumber of block threads.\n");
        printf("-v\n\tVerbose mode\n");
        printf("\nExamples:\n");
        printf("cat test.txt | ./brainWords --inputAddress addr.txt --inputIn [-v] [--bits 4/8/16..20...24] [-b blocks] [-t threads]\n");
        printf("./brainWords -v --eth --inputPhrase ../testEth.txt --inputAddress ../addrEth.txt\n");
        printf("./brainWords -v --camp2 --inputPhrase ../testEth.txt --inputAddress ../addrEth.txt\n");
        printf("./brainWords -v --eth --camp2 --inputPhrase ../testEth.txt --inputAddress ../addrEth.txt\n");
        printf("./brainWords -v --ethsha256 --inputPhrase ../testEth.txt --inputAddress ../addrEth.txt\n");
        printf("./brainWords -v --eth --root \"Hello W\" --inputAddress ../addrEth.txt\n");
        printf("./brainWords -v --ethsha256 --root \"1234\" --inputAddress ../addrEth.txt\n");
        printf("./brainWords -v --rootsuffix \"123456\" --inputAddress ../addr.txt\n");
        return -1;
    }
    cudaError_t cudaStatus;
    if (!checkDevice()) {
        return -1;
    }

    if (!readFileAddress(fileAddress)) {
        printf("Error: Loading addresses failed!\n");
        return -1;
    }
    else if (hashCount == 0) {
        printf("Error: Targets not loaded!\n");
        return -1;
    }
    printf("Loaded addresses: %zd\n", hashCount);
    printf("Number of blocks: %d\n", BLOCK_NUMBER);
    printf("Number of threads: %d\n", BLOCK_THREADS);
    printf("Number of checks per thread: %d\n", IS_BRUTE? THREAD_STEPS_BRUTE: THREAD_STEPS);
    cudaStatus = prepareCuda();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "Error: Device prepare failed: %s\n", cudaGetErrorString(cudaStatus));
        return 1;
    }

    secp = new Secp256K1();
    secp->Init();

    std::time_t s_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << "Work started at " << std::ctime(&s_time);

    if (IS_BRUTE) {
        char buf[64] = { 0 };
        int bufStart = _lltoa(SUFFIX, buf);
        string source = "";
        if (IS_ROOT) {
            source += ROOT;
            while (bufStart < 64) {
                source += buf[bufStart++];
            }            
        }
        else if (IS_ROOT_SUFFIX) {            
            while (bufStart < 64) {
                source += buf[bufStart++];
            }
            source += ROOT_SUFFIX;
        }
        printf("\nStarting text: %s\n", source.c_str());
        if (IS_ETH || IS_CAMP2) {
            cudaStatus = processCudaBruteEth();
        }
        else {
            cudaStatus = processCudaBrute();
        }
    }
    else {
        if (PHRASE_IN) {
            printf("using phrases from external input\n");
            if (IS_ETH || IS_CAMP2) {
                cudaStatus = processCuda2Eth(cin);
            }
            else {
                cudaStatus = processCuda2(cin);
            }
        }
        else {
            printf("using phrases file: %s\n", filePhrase.c_str());
            std::ifstream stream(filePhrase.c_str(), std::ios::binary);
            std::cout << "Opening file '" << filePhrase << "'" << std::endl;
            if (!stream) {
                std::cout << "Error: Failed to open file '" << filePhrase << "'" << std::endl;
                return -1;
            }
            if (IS_ETH || IS_CAMP2) {
                cudaStatus = processCuda2Eth(stream);
            }
            else {
                cudaStatus = processCuda2(stream);
            }
        }
    }

    s_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << "Work finished at " << std::ctime(&s_time);

    return 0;
}

bool readArgs(int argc, char** argv) {
    int a = 1;
    while (a < argc) {
        if (strcmp(argv[a], "-d") == 0) {
            a++;
            DEVICE_NR = strtol(argv[a], NULL, 10);
        }
        else if (strcmp(argv[a], "-fp") == 0) {
            a++;
            opt_show_false = true;
        }
        else if (strcmp(argv[a], "-v") == 0) {
            IS_VERBOSE = true;
        }
        else if (strcmp(argv[a], "-o") == 0) {
            a++;
            fileResult = string(argv[a]);
        }
        else if (strcmp(argv[a], "-os") == 0) {
            a++;
            fileStatus = string(argv[a]);
        }
        else if (strcmp(argv[a], "-t") == 0) {
            a++;
            BLOCK_THREADS = strtol(argv[a], NULL, 10);
        }
        else if (strcmp(argv[a], "-b") == 0) {
            a++;
            BLOCK_NUMBER = strtol(argv[a], NULL, 10);
        }
        else if (strcmp(argv[a], "--bbits") == 0) {
            a++;
            opt_bloom_bits = strtol(argv[a], NULL, 10);
        }
        else if (strcmp(argv[a], "--bits") == 0) {
            a++;
            PARAM_ECMULT_WINDOW_SIZE = strtol(argv[a], NULL, 10);
        }
        else if (strcmp(argv[a], "--fstatus") == 0) {
            a++;
            fileStatus = string(argv[a]);
        }
        else if (strcmp(argv[a], "--inputAddress") == 0) {
            a++;
            fileAddress = string(argv[a]);
        }
        else if (strcmp(argv[a], "--inputPhrase") == 0) {
            a++;
            filePhrase = string(argv[a]);
        }
        else if (strcmp(argv[a], "--inputIn") == 0) {
            PHRASE_IN = true;
        }
        else if (strcmp(argv[a], "--iteration") == 0) {
            a++;
            HASH_ITERATION = strtol(argv[a], NULL, 10);
        }
        else if (strcmp(argv[a], "--eth") == 0) {
            IS_ETH = true;
        }
        else if (strcmp(argv[a], "--camp2") == 0) {
            IS_CAMP2 = true;
        }
        else if (strcmp(argv[a], "--ethsha256") == 0) {
            IS_ETH_SHA256 = true;
            IS_ETH = true;
        }
        else if (strcmp(argv[a], "--root") == 0) {
            a++;
            ROOT = string(argv[a]);
            ROOT_LEN = ROOT.size();
            IS_ROOT = true;
            IS_BRUTE = true;
        }
        else if (strcmp(argv[a], "--rootsuffix") == 0) {
            a++;
            ROOT_SUFFIX = string(argv[a]);
            ROOT_LEN = ROOT_SUFFIX.size();
            IS_ROOT_SUFFIX = true;
            IS_BRUTE = true;
        }
        else if (strcmp(argv[a], "--suffix") == 0) {
            a++;
            char* end;
            string s = string(argv[a]);
            SUFFIX = strtoull(string(argv[a]).c_str(), &end, 10);
            IS_BRUTE = true;
        }
        a++;
    }
    return true;
}

void readHex(char* buf, const char* txt) {
    char b[3] = "00";
    for (unsigned int i = 0; i < strlen(txt); i += 2) {
        b[0] = *(txt + i);
        b[1] = *(txt + i + 1);
        *(buf + (i >> 1)) = strtoul(b, NULL, 16);
    }
}

char* __strlwr(char* str)
{
    unsigned char* p = (unsigned char*)str;

    while (*p) {
        *p = tolower((unsigned char)*p);
        p++;
    }

    return str;
}

bool readFileAddress(const std::string& file_name) {
    std::ifstream stream(file_name);
    std::cout << "Opening file '" << file_name << "'" << std::endl;
    if (stream.fail() || !stream.good())
    {
        std::cout << "Error: Failed to open file '" << file_name << "'" << std::endl;
        return false;
    }
    std::string buffer;
    int nr = 0;
    while (!stream.eof() && stream.good() && stream.peek() != EOF)
    {
        std::getline(stream, buffer);
        // Convert addresses to plain hex format.
        if (IS_ETH || IS_CAMP2) {
            char* b = (char*) buffer.c_str();
            __strlwr(b);
            if (buffer.at(0) == '0' && buffer.at(1) == 'x') {
                buffer = buffer.substr(2);
            }
        }
        // Add addresses.
        addresses.insert(buffer);
        // Add targets.
        std::vector<unsigned char> r160;
        hash160_20 h;
        if (buffer.length() == 40) {
            char buf[20];
            const char* hex = buffer.c_str();
            readHex(&buf[0], hex);
            for (int i = 0; i < 20 && i < buffer.length(); i++) {
                h.bits20[i] = buf[i];
            }
        }
        else if (buffer.at(0) == '1' || buffer.at(0) == '3') {
            if (DecodeBase58(buffer, r160)) {
                for (int i = 1; i <= 20 && i < r160.size(); i++) {
                    h.bits20[i - 1] = r160.at(i);
                }
            }
        }
        else if (buffer.at(0) == 'b') {
            size_t data_len = 0;
            uint8_t* data = new uint8_t[64];
            if (bech32_decode_my(data, &data_len, buffer.c_str())) {
                for (int i = 0; i < data_len && i < 20; i++) {
                    h.bits20[i] = data[i];
                }
            }
        }
        else {
            fprintf(stderr, "Error: Invalid address format: %s!", buffer.c_str());
            return false;
        }
        // std::cout << "Adding: " << buffer.c_str() << std::endl;
        _targets.push_back(h.bits5);
        nr++;
    }
    stream.close();
    hashCount = std::min(nr, (int)_targets.size());
    return true;
}

bool checkDevice() {
    cudaError_t cudaStatus = cudaSetDevice(DEVICE_NR);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "device %d failed!", DEVICE_NR);
        return false;
    }
    else {
        cudaDeviceProp props;
        cudaStatus = cudaGetDeviceProperties(&props, DEVICE_NR);
        fprintf(stderr, "Using device %d:\n", DEVICE_NR);
        fprintf(stderr, "%s (%2d procs)\n", props.name, props.multiProcessorCount);
        //printf("maxThreadsPerBlock: %2d\n\n", props.maxThreadsPerBlock);
        if (IS_BRUTE) {
            if (BLOCK_NUMBER == 0) {
                BLOCK_NUMBER = props.multiProcessorCount * 4;
            }
            if (BLOCK_THREADS == 0) {
                BLOCK_THREADS = props.maxThreadsPerBlock / 8 * 3;
            }
            workSize = (uint64_t)BLOCK_NUMBER * BLOCK_THREADS * THREAD_STEPS_BRUTE;
        }
        else {
            if (BLOCK_NUMBER == 0) {
                BLOCK_NUMBER = props.multiProcessorCount * 4;
            }
            if (BLOCK_THREADS == 0) {
                BLOCK_THREADS = props.maxThreadsPerBlock / 8 * 2;
            }
            if (IS_CAMP2) {
                BLOCK_NUMBER = props.multiProcessorCount;
                BLOCK_THREADS = props.maxThreadsPerBlock / 8;
            }
            workSize = (uint64_t)BLOCK_NUMBER * BLOCK_THREADS * THREAD_STEPS;
        }
    }
    return true;
}

cudaError_t prepareCuda() {
    unsigned int windows;
    unsigned int bits;
    size_t window_size;
    std::time_t s_time;
    cudaError_t cudaStatus;


#ifdef ECMULT_BIG_TABLE
    bits = PARAM_ECMULT_WINDOW_SIZE;
    printf("prec big gen %d bit, please wait\n", bits);
    s_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << "Gen started at " << std::ctime(&s_time);
    windows = (256 / bits) + 1;
    window_size = (1 << (bits - 1));
    cudaStatus = cudaMallocPitch(&_dev_precomp, &pitch, sizeof(secp256k1_ge_storage) * window_size, windows);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "Cuda prepare failed: kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    secp256k1_gej* _dev_gej_temp;// = new secp256k1_gej[WINDOW_SIZE];
    cudaStatus = cudaMalloc((void**)&_dev_gej_temp, window_size * sizeof(secp256k1_gej));
    secp256k1_fe* _dev_z_ratio;// = new secp256k1_fe[WINDOW_SIZE];
    cudaStatus = cudaMalloc((void**)&_dev_z_ratio, window_size * sizeof(secp256k1_fe));

    cudaStatus = loadWindow(bits, windows);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "Cuda prepare failed: loadWindow(bits) %d!\n", cudaStatus);
        goto Error;
    }

    ecmult_big_create << <1, 1 >> > (_dev_gej_temp, _dev_z_ratio, _dev_precomp, pitch, bits);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "Cuda prepare failed: kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "Cuda prepare failed: cudaDeviceSynchronize returned error code %d after launching kernel 'computeTable'!\n", cudaStatus);
        goto Error;
    }
    printf("prec gen %d bit finished\n", bits);
#else
    int bits = ECMULT_GEN_PREC_BITS;
    if (bits != 8 && bits != 4 && bits != 2) {
        int g = ECMULT_GEN_PREC_G;
        int n = ECMULT_GEN_PREC_N;
        printf("prec gen %d bit, please wait\n", bits);
        secp256k1_ge_storage* _dev_prec_table = new secp256k1_ge_storage[g * n];
        secp256k1_ge* _dev_prec = new secp256k1_ge[g * n];
        secp256k1_gej* _dev_precj = new secp256k1_gej[g * n];

        cudaStatus = cudaMalloc((void**)&_dev_precj, g * n * sizeof(secp256k1_gej));
        cudaStatus = cudaMalloc((void**)&_dev_prec, g * n * sizeof(secp256k1_ge));
        cudaStatus = cudaMalloc((void**)&_dev_prec_table, g * n * sizeof(secp256k1_ge_storage));

        computeTable << <1, 1 >> > (_dev_prec_table, _dev_prec, _dev_precj);
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "Cuda prepare failed: kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
            goto Error;
        }
        cudaStatus = cudaDeviceSynchronize();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "Cuda prepare failed: cudaDeviceSynchronize returned error code %d after launching kernel 'computeTable'!\n", cudaStatus);
            goto Error;
        }
        printf("prec gen %d bit finished\n", bits);
    }
#endif

    fprintf(stderr, "Setting targets [%lld].\n", _targets.size());
    cudaStatus = _targetLookup.setTargets(_targets, opt_bloom_bits);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "\nkernel '_targetLookup.setTarget' launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    cudaStatus = loadIteration(HASH_ITERATION);
    if (cudaSuccess != cudaStatus) {
        fprintf(stderr, "\nCuda prepare failed: 'loadIteration' %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    _targets.clear();
    if (IS_BRUTE) {
        int shaLevel = 0;
        if (IS_ROOT) {
            cudaStatus = loadPrefix(ROOT.c_str(), ROOT_LEN);
            if (cudaSuccess != cudaStatus) {
                fprintf(stderr, "\nCuda prepare failed: 'loadPrefix' %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
            shaLevel = ROOT_LEN / 4;
            cudaStatus = loadLevel(shaLevel);
            if (cudaSuccess != cudaStatus) {
                fprintf(stderr, "\nCuda prepare failed: 'loadLevel' %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
        }
        else if(IS_ROOT_SUFFIX) {
            cudaStatus = loadPrefix(ROOT_SUFFIX.c_str(), ROOT_LEN);
            if (cudaSuccess != cudaStatus) {
                fprintf(stderr, "\nCuda prepare failed: 'loadPrefix' %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
            cudaStatus = loadLevel(0);
            if (cudaSuccess != cudaStatus) {
                fprintf(stderr, "\nCuda prepare failed: 'loadLevel' %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
        }
        if (shaLevel > 0) {
            uint32_t* _shapre = new uint32_t[8];
            cudaStatus = cudaMallocManaged(&_shapre, 8 * sizeof(uint32_t));
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "\ncudaMemcpy 'shaPre' failed: %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
            shaPre << <1, 1 >> > (ROOT_LEN, _shapre);
            cudaStatus = cudaGetLastError();
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "\nkernel launch 'shaPre' failed: %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
            cudaStatus = cudaDeviceSynchronize();
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "\ncudaDeviceSynchronize returned error code %d after launching kernel 'shaPre'!\n", cudaStatus);
                goto Error;
            }
            else {
                loadShaPre(_shapre);
            }
        }
    }
    printf("\nCUDA preparation finished.\n");

#ifdef ECMULT_BIG_TABLE
    Error:
        cudaFree(_dev_gej_temp);
        cudaFree(_dev_z_ratio);
#else
    Error:
        cudaFree(_dev_precj);
        cudaFree(_dev_prec);
        cudaFree(_dev_prec_table);
#endif

    return cudaStatus;
}

cudaError_t processCuda2Eth(std::istream& stream) {
    cudaError_t cudaStatus;

    uint32_t outputSizeB_T = BLOCK_NUMBER * BLOCK_THREADS * THREAD_STEPS;

    bool* buffIsResult = new bool[1];
    cudaStatus = cudaMallocManaged(&buffIsResult, 1 * sizeof(bool));
    buffIsResult[0] = false;

    bool* buffDeviceResult = new bool[outputSizeB_T];
    cudaStatus = cudaMallocManaged(&buffDeviceResult, outputSizeB_T * sizeof(bool));
    for (uint32_t i = 0; i < outputSizeB_T; i++) {
        buffDeviceResult[i] = false;
    }

    std::string buffer;
    std::string combined1;
    std::vector<uint32_t> indexes1;
    std::string combined2;
    std::vector<uint32_t> indexes2;

    std::chrono::steady_clock::time_point beginCountHashrate = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point beginCountStatus = std::chrono::steady_clock::now();
    uint64_t counter = 0;

    bool isData = true;
    bool data1ready = false;
    bool data2ready = false;
    uint32_t* _devIndexes1;
    uint32_t* _devIndexes2;
    char* _devLines1;
    char* _devLines2;
    uint32_t nr1 = 0;
    uint32_t nr2 = 0;
    string nr1Last = "";
    string nr2Last = "";

    cudaMalloc((void**)&_devIndexes1, sizeof(uint32_t));
    cudaMalloc((void**)&_devIndexes2, sizeof(uint32_t));
    cudaMalloc((void**)&_devLines1, 1);
    cudaMalloc((void**)&_devLines2, 1);

    short workType = 0;
    if (IS_ETH) {
        workType += 1;
    }
    if (IS_CAMP2) {
        workType += 2;
    }

    while (isData) {
        //read file 1
        combined1.clear();
        indexes1.clear();
        nr1 = 0;
        data1ready = false;
        std::chrono::steady_clock::time_point beginCountRead = std::chrono::steady_clock::now();
        while (std::getline(stream, buffer))
        {
            if (buffer.length() == 0) {
                continue;
            }
            if (buffer.length() > 63) {
                buffer = buffer.substr(0, 63);
            }
            combined1 += buffer;
            nr1Last = buffer;
            indexes1.emplace_back(combined1.size());
            nr1++;
            if (nr1 < outputSizeB_T) {
                continue;
            }
            data1ready = true;
            break;
        }
        if (!data1ready) {
            isData = false;
        }

        //copy data 1 -> device
        if (nr1 > 0) {
            cudaFree(_devLines1);
            cudaMalloc((void**)&_devLines1, combined1.size());
            cudaStatus = cudaMemcpyAsync(_devLines1, combined1.data(), combined1.size(), cudaMemcpyHostToDevice);
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaMemcpy combined1.data() failed: %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
            cudaFree(_devIndexes1);
            cudaMalloc((void**)&_devIndexes1, indexes1.size() * sizeof(uint32_t));
            cudaStatus = cudaMemcpyAsync(_devIndexes1, indexes1.data(), indexes1.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaMemcpy indexes1.data() failed: %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
        }
        //show status
        long long tHash = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - beginCountHashrate).count();
        if (tHash > 5000) {
            printSpeed((double)((double)counter / tHash) / 1000.0);
            counter = 0;
            beginCountHashrate = std::chrono::steady_clock::now();
            if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - beginCountStatus).count() >= fileStatusInterval) {
                saveStatus();
                beginCountStatus = std::chrono::steady_clock::now();
            }
        }

        //sync kernel 2
        if (nr2 > 0) {
            cudaStatus = cudaDeviceSynchronize();
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching kernel!\n", cudaStatus);
                goto Error;
            }
            //check result 2
            if (buffIsResult[0]) {
                buffIsResult[0] = false;
                checkResultEth(buffDeviceResult, outputSizeB_T, combined2, indexes2);
            }
        }

        //launch kernel 1
        if (nr1 > 0) {
            //worker << <BLOCK_NUMBER, BLOCK_THREADS >> > (buffIsResult, buffDeviceResult, _devLines1, _devIndexes1, nr1, _dev_precomp, pitch);
            workerEth << <BLOCK_NUMBER, BLOCK_THREADS >> > (buffIsResult, buffDeviceResult, _devLines1, _devIndexes1, nr1, workType, IS_ETH_SHA256, _dev_precomp, pitch);
            cudaStatus = cudaGetLastError();
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
        }

        STATUS = nr2Last;
        counter += nr2;
        counterTotal += nr2;

        //read file 1
        combined2.clear();
        indexes2.clear();
        nr2 = 0;
        data2ready = false;
        while (std::getline(stream, buffer))
        {
            if (buffer.length() == 0) {
                continue;
            }
            if (buffer.length() > 63) {
                buffer = buffer.substr(0, 63);
            }
            combined2 += buffer;
            nr2Last = buffer;
            indexes2.emplace_back(combined2.size());
            nr2++;
            if (nr2 < outputSizeB_T) {
                continue;
            }
            data2ready = true;
            break;
        }
        if (!data2ready) {
            isData = false;
        }
        if (nr2 > 0) {
            //copy data 2 -> device
            cudaFree(_devLines2);
            cudaMalloc((void**)&_devLines2, combined2.size());
            cudaStatus = cudaMemcpyAsync(_devLines2, combined2.data(), combined2.size(), cudaMemcpyHostToDevice);
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaMemcpy combined2.data() failed: %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
            cudaFree(_devIndexes2);
            cudaMalloc((void**)&_devIndexes2, indexes2.size() * sizeof(uint32_t));
            cudaStatus = cudaMemcpyAsync(_devIndexes2, indexes2.data(), indexes2.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaMemcpy indexes2.data() failed: %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
        }
        //sync kernel 1
        if (nr1 > 0) {
            cudaStatus = cudaDeviceSynchronize();
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching kernel!\n", cudaStatus);
                goto Error;
            }
            //check result 1
            if (buffIsResult[0]) {
                buffIsResult[0] = false;
                checkResultEth(buffDeviceResult, outputSizeB_T, combined1, indexes1);
            }
        }
        //launch kernel 2
        if (nr2 > 0) {
            //worker << <BLOCK_NUMBER, BLOCK_THREADS >> > (buffIsResult, buffDeviceResult, _devLines2, _devIndexes2, nr2, _dev_precomp, pitch);
            workerEth << <BLOCK_NUMBER, BLOCK_THREADS >> > (buffIsResult, buffDeviceResult, _devLines2, _devIndexes2, nr2, workType, IS_ETH_SHA256, _dev_precomp, pitch);
            cudaStatus = cudaGetLastError();
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
        }
        STATUS = nr1Last;
        counter += nr1;
        counterTotal += nr1;
    }

Error:
    cudaFree(_dev_precomp);
    cudaFree(_devIndexes1);
    cudaFree(_devIndexes2);
    cudaFree(_devLines1);
    cudaFree(_devLines2);
    return cudaStatus;
}

cudaError_t processCuda2(std::istream& stream) {
    cudaError_t cudaStatus;

    uint32_t outputSizeB_T = BLOCK_NUMBER * BLOCK_THREADS * THREAD_STEPS;

    bool* buffIsResult = new bool[1];
    cudaStatus = cudaMallocManaged(&buffIsResult, 1 * sizeof(bool));
    buffIsResult[0] = false;

    bool* buffDeviceResult = new bool[outputSizeB_T];
    cudaStatus = cudaMallocManaged(&buffDeviceResult, outputSizeB_T * sizeof(bool));
    for (uint32_t i = 0; i < outputSizeB_T; i++) {
        buffDeviceResult[i] = false;
    }

    std::string buffer;
    std::string combined1;
    std::vector<uint32_t> indexes1;
    std::string combined2;
    std::vector<uint32_t> indexes2;

    std::chrono::steady_clock::time_point beginCountHashrate = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point beginCountStatus = std::chrono::steady_clock::now();
    uint64_t counter = 0;

    bool isData = true;
    bool data1ready = false;
    bool data2ready = false;
    uint32_t* _devIndexes1;
    uint32_t* _devIndexes2;
    char* _devLines1;
    char* _devLines2;
    uint32_t nr1 = 0;
    uint32_t nr2 = 0;
    string nr1Last = "";
    string nr2Last = "";

    cudaMalloc((void**)&_devIndexes1, sizeof(uint32_t));
    cudaMalloc((void**)&_devIndexes2, sizeof(uint32_t));
    cudaMalloc((void**)&_devLines1, 1);
    cudaMalloc((void**)&_devLines2, 1);

    while (isData) {
        //read file 1
        combined1.clear();
        indexes1.clear();
        nr1 = 0;
        data1ready = false;
        std::chrono::steady_clock::time_point beginCountRead = std::chrono::steady_clock::now();
        while (std::getline(stream, buffer))
        {
            if (buffer.length() == 0) {
                continue;
            }
            if (buffer.length() > 63) {
                buffer = buffer.substr(0, 63);
            }
            combined1 += buffer;
            nr1Last = buffer;
            indexes1.emplace_back(combined1.size());
            nr1++;
            if (nr1 < outputSizeB_T) {
                continue;
            }
            data1ready = true;
            break;
        }
        if (!data1ready) {
            isData = false;
        }

        //copy data 1 -> device
        if (nr1 > 0) {
            cudaFree(_devLines1);
            cudaMalloc((void**)&_devLines1, combined1.size());
            cudaStatus = cudaMemcpyAsync(_devLines1, combined1.data(), combined1.size(), cudaMemcpyHostToDevice);
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaMemcpy combined1.data() failed: %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
            cudaFree(_devIndexes1);
            cudaMalloc((void**)&_devIndexes1, indexes1.size() * sizeof(uint32_t));
            cudaStatus = cudaMemcpyAsync(_devIndexes1, indexes1.data(), indexes1.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaMemcpy indexes1.data() failed: %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
        }
        //show status
        long long tHash = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - beginCountHashrate).count();
        if (tHash > 5000) {
            printSpeed((double)((double)counter / tHash) / 1000.0);
            counter = 0;
            beginCountHashrate = std::chrono::steady_clock::now();
            if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - beginCountStatus).count() >= fileStatusInterval) {
                saveStatus();
                beginCountStatus = std::chrono::steady_clock::now();
            }
        }

        //sync kernel 2
        if (nr2 > 0) {
            cudaStatus = cudaDeviceSynchronize();
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching kernel!\n", cudaStatus);
                goto Error;
            }
            //check result 2
            if (buffIsResult[0]) {
                buffIsResult[0] = false;
                checkResult(buffDeviceResult, outputSizeB_T, combined2, indexes2);
            }
        }

        //launch kernel 1
        if (nr1 > 0) {
            worker << <BLOCK_NUMBER, BLOCK_THREADS >> > (buffIsResult, buffDeviceResult, _devLines1, _devIndexes1, nr1, _dev_precomp, pitch);
            cudaStatus = cudaGetLastError();
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
        }

        STATUS = nr2Last;
        counter += nr2;
        counterTotal += nr2;

        //read file 1
        combined2.clear();
        indexes2.clear();
        nr2 = 0;
        data2ready = false;
        while (std::getline(stream, buffer))
        {
            if (buffer.length() == 0) {
                continue;
            }
            if (buffer.length() > 63) {
                buffer = buffer.substr(0, 63);
            }
            combined2 += buffer;
            nr2Last = buffer;
            indexes2.emplace_back(combined2.size());
            nr2++;
            if (nr2 < outputSizeB_T) {
                continue;
            }
            data2ready = true;
            break;
        }
        if (!data2ready) {
            isData = false;
        }
        if (nr2 > 0) {
            //copy data 2 -> device
            cudaFree(_devLines2);
            cudaMalloc((void**)&_devLines2, combined2.size());
            cudaStatus = cudaMemcpyAsync(_devLines2, combined2.data(), combined2.size(), cudaMemcpyHostToDevice);
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaMemcpy combined2.data() failed: %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
            cudaFree(_devIndexes2);
            cudaMalloc((void**)&_devIndexes2, indexes2.size() * sizeof(uint32_t));
            cudaStatus = cudaMemcpyAsync(_devIndexes2, indexes2.data(), indexes2.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaMemcpy indexes2.data() failed: %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
        }
        //sync kernel 1
        if (nr1 > 0) {
            cudaStatus = cudaDeviceSynchronize();
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching kernel!\n", cudaStatus);
                goto Error;
            }
            //check result 1
            if (buffIsResult[0]) {
                buffIsResult[0] = false;
                checkResult(buffDeviceResult, outputSizeB_T, combined1, indexes1);
            }
        }
        //launch kernel 2
        if (nr2 > 0) {
            worker << <BLOCK_NUMBER, BLOCK_THREADS >> > (buffIsResult, buffDeviceResult, _devLines2, _devIndexes2, nr2, _dev_precomp, pitch);
            cudaStatus = cudaGetLastError();
            if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
                goto Error;
            }
        }
        STATUS = nr1Last;
        counter += nr1;
        counterTotal += nr1;
    }

Error:
    cudaFree(_dev_precomp);
    cudaFree(_devIndexes1);
    cudaFree(_devIndexes2);
    cudaFree(_devLines1);
    cudaFree(_devLines1);
    return cudaStatus;
}

int _lltoa(uint64_t val, char* buf) {
    int i = 63;
    for (; val && i; --i, val /= ALPHABET_LEN) {
        buf[i] = ALPHABET[val % ALPHABET_LEN];
    }
    return i + 1;
}

void checkResultEth(bool* buffDeviceResult, uint64_t outputSizeB_T, string combined, vector<uint32_t> indexes) {
    for (uint32_t i = 0; i < indexes.size(); i++) {
        if (buffDeviceResult[i]) {
            uint64_t s = i == 0 ? 0 : indexes.at(i - 1);
            uint64_t len = indexes.at(i) - s;
            string source;
            size_t ix = 0;
            for (; ix < len; ix++) {
                source += combined.at(s++);
                if (ix == 63) {
                    break;
                }
            }
            for (int turn = 1; turn <= 1 + (IS_CAMP2 ? 1 : 0); turn++) {
                Int privKey;
                privKey.SetInt32(0);
                int32_t rewrittenKey[NB32BLOCK];
                if (IS_ETH_SHA256) {
                    sha256((unsigned char*)source.c_str(), source.size(), (unsigned char*)privKey.bits64);
                    for (int hIter = 1; hIter < HASH_ITERATION; hIter++) {
                        unsigned char toHash[32];
                        for (int rk = 0; rk < 8; rk++) {
                            int rkk = 4 * rk;
                            toHash[rkk] = (unsigned char)((privKey.bits[rk] & 0x000000ff));
                            toHash[rkk + 1] = (unsigned char)((privKey.bits[rk] & 0x0000ff00) >> 8);
                            toHash[rkk + 2] = (unsigned char)((privKey.bits[rk] & 0x00ff0000) >> 16);
                            toHash[rkk + 3] = (unsigned char)((privKey.bits[rk] & 0xff000000) >> 24);
                        }
                        sha256(toHash, 32, (unsigned char*)privKey.bits);
                    }
                    for (int rk = 0; rk < 8; rk++) {
                        rewrittenKey[rk] = (privKey.bits[rk] & 0x000000ff) << 24;
                        rewrittenKey[rk] |= (privKey.bits[rk] & 0x0000ff00) << 8;
                        rewrittenKey[rk] |= (privKey.bits[rk] & 0x00ff0000) >> 8;
                        rewrittenKey[rk] |= (privKey.bits[rk] & 0xff000000) >> 24;
                        privKey.bits[rk] = rewrittenKey[rk];
                    }
                }
                else {
                    unsigned char d_hash[32];
                    _keccak(source.c_str(), source.size(), d_hash, 32);
                    if (turn == 2) {
                        char toHash[32];
#pragma unroll 2030
                        for (int camp2Turns = 0; camp2Turns < 2030; camp2Turns++) {
#pragma unroll 32
                            for (int xxx = 0; xxx < 32; xxx++) {
                                toHash[xxx] = d_hash[xxx];
                            }
                            _keccak(toHash, 32, d_hash, 32);
                        }
                    }
                    for (int rk = 0; rk < 8; rk++) {
                        rewrittenKey[rk] = ((int32_t)d_hash[rk * 4]) << 24;
                        rewrittenKey[rk] |= (int32_t)(d_hash[rk * 4 + 1]) << 16;
                        rewrittenKey[rk] |= ((int32_t)d_hash[rk * 4 + 2]) << 8;
                        rewrittenKey[rk] |= ((int32_t)d_hash[rk * 4 + 3]);
                        privKey.bits[rk] = rewrittenKey[rk];
                    }
                }
                privKey.SetInt32(0);
                for (int rk = 0; rk < 8; rk++) {
                    privKey.bits[7 - rk] = rewrittenKey[rk];
                }
                if (processCandidateEth(privKey, source, turn==2)) {
                    printf("found: %s \n", source.c_str());
                }
            }
            buffDeviceResult[i] = false;
        }
    }
}

void checkResult(bool * buffDeviceResult, uint64_t outputSizeB_T, string combined, vector<uint32_t> indexes) {
    for (uint32_t i = 0; i < indexes.size(); i++) {
        if (buffDeviceResult[i]) {
            uint64_t s = i == 0 ? 0 : indexes.at(i - 1);
            uint64_t len = indexes.at(i) - s;
            string source;
            size_t ix = 0;
            for (; ix < len; ix++) {
                source += combined.at(s++);
                if (ix == 63) {
                    break;
                }
            }
            Int privKey;
            privKey.SetInt32(0);
            int32_t rewrittenKey[NB32BLOCK];
            sha256((unsigned char*)source.c_str(), source.size(), (unsigned char*)privKey.bits64);
            for (int hIter = 1; hIter < HASH_ITERATION; hIter++) {
                unsigned char toHash[32];
                for (int rk = 0; rk < 8; rk++) {
                    int rkk = 4 * rk;
                    toHash[rkk] = (unsigned char)((privKey.bits[rk] & 0x000000ff) );
                    toHash[rkk +1] = (unsigned char)((privKey.bits[rk] & 0x0000ff00) >> 8);
                    toHash[rkk +2] = (unsigned char)((privKey.bits[rk] & 0x00ff0000) >> 16);
                    toHash[rkk +3] = (unsigned char)((privKey.bits[rk] & 0xff000000) >> 24);
                }
                sha256(toHash, 32, (unsigned char*)privKey.bits);
            }
            for (int rk = 0; rk < 8; rk++) {
                rewrittenKey[rk] = (privKey.bits[rk] & 0x000000ff) << 24;
                rewrittenKey[rk] |= (privKey.bits[rk] & 0x0000ff00) << 8;
                rewrittenKey[rk] |= (privKey.bits[rk] & 0x00ff0000) >> 8;
                rewrittenKey[rk] |= (privKey.bits[rk] & 0xff000000) >> 24;
                privKey.bits[rk] = rewrittenKey[rk];
            }
            privKey.SetInt32(0);
            for (int rk = 0; rk < 8; rk++) {
                privKey.bits[7 - rk] = rewrittenKey[rk];
            }
            if (processCandidate(privKey, source)) {
                printf("found: %s \n", source.c_str());
            }
            buffDeviceResult[i]=false;
        }
    }
}

bool processCandidateEth(Int& toTest, string text, bool isCamp2) {
    char address[50];
    Point publickey = secp->ComputePublicKey(&toTest);
    char pubKey[64];
    char buf[32];
    readHex(&buf[0], publickey.x.GetBase16().c_str());
    for (int i = 0; i < 32; i++) {
        pubKey[i] = buf[i];
    }
    readHex(&buf[0], publickey.y.GetBase16().c_str());
    for (int i = 0; i < 32; i++) {
        pubKey[32+i] = buf[i];
    }
    unsigned char d_hash[32];
    _keccak(pubKey, 64, d_hash, 32);
    uint32_t hash[5];
    for (int h = 12, i = 0; i < 5; i++) {
        hash[i] = ((d_hash[h++] << 24) & 0xff000000);
        hash[i] |= ((d_hash[h++] << 16) & 0x00ff0000);
        hash[i] |= ((d_hash[h++] << 8) & 0x0000ff00);
        hash[i] |= (d_hash[h++] );
    }
    sprintf(address, "%08x%08x%08x%08x%08x", hash[0], hash[1], hash[2], hash[3], hash[4]);
    string a(address);
    if (addresses.find(a) != addresses.end()) {
        FILE* keys;
        printf("\nfound: %s - %s (%s)\n", address, toTest.GetBase16().c_str(), isCamp2?"c2":"1");
        keys = fopen(fileResult.c_str(), "a+");
        fprintf(keys, "%s (%s)\n", address, isCamp2 ? "c2" : "1");
        fprintf(keys, "%s\n", text.c_str());
        fprintf(keys, "%s\n\n", toTest.GetBase16().c_str());
        fclose(keys);
        return true;
    }
    return false;
}

bool processCandidate(Int& toTest, string text) {
    char rmdhash[21], address[50];
    Point publickey = secp->ComputePublicKey(&toTest);
    bool result = false;
    bool rmdResult = false;
    bool uncompResult = false;
    bool compResult = false;
    bool bech32Result= false;
    bool p2shResult = false;
    {
        secp->GetHash160(P2PKH, false, publickey, (unsigned char*)rmdhash);
        string a = (char*)tohex((char*)rmdhash, 20);
        if (addresses.find(a) != addresses.end()) {
            FILE* keys;
            printf("\nfound: %s - %s (%s)\n", a.c_str(), toTest.GetBase16().c_str(), "H");
            keys = fopen(fileResult.c_str(), "a+");
            fprintf(keys, "%s (%s)\n", a.c_str(), "H");
            fprintf(keys, "%s\n", text.c_str());
            fprintf(keys, "%s\n\n", toTest.GetBase16().c_str());
            fclose(keys);
            rmdResult = true;
            return true;
        }
    }
    {
        secp->GetHash160(P2PKH, false, publickey, (unsigned char*)rmdhash);
        addressToBase58(rmdhash, address, false);
        string a = address;
        if (addresses.find(a) != addresses.end()) {
            FILE* keys;
            printf("\nfound: %s - %s (U)\n", address, toTest.GetBase16().c_str());
            keys = fopen(fileResult.c_str(), "a+");
            fprintf(keys, "%s (U)\n", address);
            fprintf(keys, "%s\n", text.c_str());
            fprintf(keys, "%s\n\n", toTest.GetBase16().c_str());
            fclose(keys);
            uncompResult = true;
            //return true;
        }
    }
    {
        secp->GetHash160(P2PKH, true, publickey, (unsigned char*)rmdhash);
        addressToBase58(rmdhash, address, false);
        string a = address;
        if (addresses.find(a) != addresses.end()) {
            FILE* keys;
            printf("\nfound: %s - %s (C)\n", address, toTest.GetBase16().c_str());
            keys = fopen(fileResult.c_str(), "a+");
            fprintf(keys, "%s (C)\n", address);
            fprintf(keys, "%s\n", text.c_str());
            fprintf(keys, "%s\n\n", toTest.GetBase16().c_str());
            fclose(keys);
            compResult = true;
        }
    }
    {
        char output[128];
        segwit_addr_encode(output, "bc", 0, (unsigned char*)rmdhash, 20);
        string a = output;
        if (addresses.find(a) != addresses.end()) {
            FILE* keys;
            printf("\nfound: %s - %s (B)\n", output, toTest.GetBase16().c_str());
            keys = fopen(fileResult.c_str(), "a+");
            fprintf(keys, "%s (B)\n", output);
            fprintf(keys, "%s\n", text.c_str());
            fprintf(keys, "%s\n\n", toTest.GetBase16().c_str());
            fclose(keys);
            bech32Result = true;
        }
    }
    {
        secp->GetHash160(P2SH, true, publickey, (unsigned char*)rmdhash);
        addressToBase58(rmdhash, address, true);
        string a = address;
        if (addresses.find(a) != addresses.end()) {
            FILE* keys;
            printf("\nfound: %s - %s (P)\n", address, toTest.GetBase16().c_str());
            keys = fopen(fileResult.c_str(), "a+");
            fprintf(keys, "%s (P)\n", address);
            fprintf(keys, "%s\n", text.c_str());
            fprintf(keys, "%s\n\n", toTest.GetBase16().c_str());
            fclose(keys);
            compResult = true;
        }
    }
    result = rmdResult|uncompResult|compResult|bech32Result|p2shResult;
    if (!result && opt_show_false) {
        secp->GetHash160(P2PKH, false, publickey, (unsigned char*)rmdhash);
        string a = (char*)tohex((char*)rmdhash, 20);
        printf("\nfalse: %s - %s (%s) - %s\n", a.c_str(), toTest.GetBase16().c_str(), "F", text.c_str());
    }
    return result;
}

void printSpeed(double speed) {
    std::string speedStr;
    if (speed < 0.01) {
        speedStr = "< 0.01 MKey/s";
    }
    else {
        if (speed < 1000) {
            speedStr = formatDouble("%.3f", speed) + " MKey/s";
        }
        else {
            speed /= 1000;
            if (speed < 1000) {
                speedStr = formatDouble("%.3f", speed) + " GKey/s";
            }
            else {
                speed /= 1000;
                speedStr = formatDouble("%.3f", speed) + " TKey/s";
            }
        }
        }
    if (IS_VERBOSE) {
        if (IS_BRUTE) {
            string source = "";
            char buf[64] = { 0 };
            int bufStart = _lltoa(SUFFIX, buf);
            if (IS_ROOT) {
                source += ROOT;
                while (bufStart < 64) {
                    source += buf[bufStart++];
                }
            }
            else {
                while (bufStart < 64) {
                    source += buf[bufStart++];
                }
                source += ROOT_SUFFIX;
            }
            printf("\r %s, %s|%lld     ", speedStr.c_str(), source.c_str(), SUFFIX);
        }
        else {
            printf("\r%s, %lld|%s   ", speedStr.c_str(), counterTotal, STATUS.c_str());
        }
    }
    else {
        printf("\r%s     ", speedStr.c_str());
    }
    fflush(stdout);
}

void saveStatusBrute() {
    FILE* stat = fopen(fileStatus.c_str(), "w");
    auto time = std::chrono::system_clock::now();
    std::time_t s_time = std::chrono::system_clock::to_time_t(time);
    fprintf(stat, "%s\n\n", std::ctime(&s_time));
    fprintf(stat, "--suffix %lld\n", SUFFIX);
    if (IS_ROOT) {
        fprintf(stat, "--root %s\n", ROOT.c_str());
    }
    else {
        fprintf(stat, "--rootsuffix %s\n", ROOT_SUFFIX.c_str());
    }
    fprintf(stat, "--fstatus %s\n", fileStatus.c_str());
    fprintf(stat, "-d %d\n", DEVICE_NR);
    fprintf(stat, "-t %d\n", BLOCK_THREADS);
    fprintf(stat, "-b %d\n", BLOCK_NUMBER);
    fprintf(stat, "--bits %d\n", PARAM_ECMULT_WINDOW_SIZE);
    fprintf(stat, "--bbits %d\n", opt_bloom_bits);
    fprintf(stat, "--inputAddress %s\n", fileAddress.c_str());
    fclose(stat);
}

void saveStatus() {
    FILE* stat = fopen(fileStatus.c_str(), "w");
    auto time = std::chrono::system_clock::now();
    std::time_t s_time = std::chrono::system_clock::to_time_t(time);
    fprintf(stat, "%s\n\n", std::ctime(&s_time));
    fprintf(stat, "lines processed: %lld\n", counterTotal);
    fprintf(stat, "last line processed: %s\n", STATUS.c_str());
    fprintf(stat, "--fstatus %s\n", fileStatus.c_str());
    fprintf(stat, "-d %d\n", DEVICE_NR);
    fprintf(stat, "-t %d\n", BLOCK_THREADS);
    fprintf(stat, "-b %d\n", BLOCK_NUMBER);
    fprintf(stat, "--bits %d\n", PARAM_ECMULT_WINDOW_SIZE);
    fprintf(stat, "--bbits %d\n", opt_bloom_bits);
    fprintf(stat, "--inputAddress %s\n", fileAddress.c_str());
    fclose(stat);
}

cudaError_t processCudaBrute() {
    cudaError_t cudaStatus;

    bool* buffIsResult = new bool[1];
    cudaStatus = cudaMallocManaged(&buffIsResult, 1 * sizeof(bool));
    buffIsResult[0] = false;

    uint32_t outputSizeB_T = BLOCK_NUMBER * BLOCK_THREADS;
    uint16_t* buffDeviceResult = new uint16_t[outputSizeB_T];
    cudaStatus = cudaMallocManaged(&buffDeviceResult, outputSizeB_T * sizeof(uint16_t));
    for (uint32_t i = 0; i < outputSizeB_T; i++) {
        buffDeviceResult[i] = NO_RESULT;
    }

    uint64_t counter = 0;
    std::chrono::steady_clock::time_point beginCountHashrate = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point beginCountStatus = std::chrono::steady_clock::now();
    while (!IS_RESULT) {
        //std::chrono::steady_clock::time_point beginTkernel = std::chrono::steady_clock::now();
        if (IS_ROOT) {
            workerBrute << <BLOCK_NUMBER, BLOCK_THREADS >> > (buffIsResult, buffDeviceResult, ROOT_LEN, SUFFIX, _dev_precomp, pitch);
        }
        else {
            workerBruteSuffix << <BLOCK_NUMBER, BLOCK_THREADS >> > (buffIsResult, buffDeviceResult, ROOT_LEN, SUFFIX, _dev_precomp, pitch);
        }
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
            goto Error;
        }

        long long tHash = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - beginCountHashrate).count();
        if (tHash > 5000) {
            printSpeed((double)((double)counter / tHash) / 1000.0);
            counter = 0;
            beginCountHashrate = std::chrono::steady_clock::now();
            if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - beginCountStatus).count() >= fileStatusInterval) {
                saveStatusBrute();
                beginCountStatus = std::chrono::steady_clock::now();
            }
        }

        cudaStatus = cudaDeviceSynchronize();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching kernel!\n", cudaStatus);
            goto Error;
        }
        //long long tKernel = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - beginTkernel).count();
        //printf(" tKernel %d\n", tKernel);
        if (buffIsResult[0]) {
            buffIsResult[0] = false;
            for (uint32_t i = 0; i < outputSizeB_T; i++) {
                if (buffDeviceResult[i] != NO_RESULT) {
                    //printf("%d %lld %u\n", c, SUFFIX, buffDeviceResult[i]); continue;
                    uint64_t base = SUFFIX;
                    base += (THREAD_STEPS_BRUTE * i) + buffDeviceResult[i];
                    string source="";
                    char buf[64] = { 0 };
                    int bufStart = _lltoa(base, buf);
                    if (IS_ROOT) {
                        source += ROOT;
                        while (bufStart < 64) {
                            source += buf[bufStart++];
                        }
                    }
                    else {
                        while (bufStart < 64) {
                            source += buf[bufStart++];
                        }
                        source += ROOT_SUFFIX;
                    }
                    Int privKey;
                    privKey.SetInt32(0);
                    sha256((unsigned char*)source.c_str(), source.size(), (unsigned char*)privKey.bits64);
                    int32_t rewrittenKey[NB32BLOCK];
                    for (int rk = 0; rk < 8; rk++) {
                        rewrittenKey[rk] = (privKey.bits[rk] & 0x000000ff) << 24;
                        rewrittenKey[rk] |= (privKey.bits[rk] & 0x0000ff00) << 8;
                        rewrittenKey[rk] |= (privKey.bits[rk] & 0x00ff0000) >> 8;
                        rewrittenKey[rk] |= (privKey.bits[rk] & 0xff000000) >> 24;
                        privKey.bits[rk] = rewrittenKey[rk];
                    }
                    privKey.SetInt32(0);
                    for (int rk = 0; rk < 8; rk++) {
                        privKey.bits[7 - rk] = rewrittenKey[rk];
                    }
                    if (processCandidate(privKey, source)) {
                        printf("found: %s | %llu\n", source.c_str(), base);
                        IS_RESULT = true;
                    }
                    buffDeviceResult[i] = NO_RESULT;
                }
            }
        }
        counter += workSize;
        SUFFIX += workSize;
    }

Error:
    cudaFree(_dev_precomp);
    return cudaStatus;
}

cudaError_t processCudaBruteEth() {
    cudaError_t cudaStatus;

    bool* buffIsResult = new bool[1];
    cudaStatus = cudaMallocManaged(&buffIsResult, 1 * sizeof(bool));
    buffIsResult[0] = false;

    uint32_t outputSizeB_T = BLOCK_NUMBER * BLOCK_THREADS;
    uint16_t* buffDeviceResult = new uint16_t[outputSizeB_T];
    cudaStatus = cudaMallocManaged(&buffDeviceResult, outputSizeB_T * sizeof(uint16_t));
    for (uint32_t i = 0; i < outputSizeB_T; i++) {
        buffDeviceResult[i] = NO_RESULT;
    }

    short workType = 0;
    if (IS_ETH) {
        workType += 1;
    }
    if (IS_CAMP2) {
        workType += 2;
    }

    uint64_t counter = 0;
    std::chrono::steady_clock::time_point beginCountHashrate = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point beginCountStatus = std::chrono::steady_clock::now();
    while (!IS_RESULT) {
        //std::chrono::steady_clock::time_point beginTkernel = std::chrono::steady_clock::now();
        if (IS_ROOT) {
            workerBruteEth << <BLOCK_NUMBER, BLOCK_THREADS >> > (buffIsResult, buffDeviceResult, ROOT_LEN, SUFFIX, workType, IS_ETH_SHA256, _dev_precomp, pitch);
        }
        else {
            workerBruteEthSuffix << <BLOCK_NUMBER, BLOCK_THREADS >> > (buffIsResult, buffDeviceResult, ROOT_LEN, SUFFIX, workType, IS_ETH_SHA256, _dev_precomp, pitch);
        }
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
            goto Error;
        }

        long long tHash = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - beginCountHashrate).count();
        if (tHash > 5000) {
            printSpeed((double)((double)counter / tHash) / 1000.0);
            counter = 0;
            beginCountHashrate = std::chrono::steady_clock::now();
            if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - beginCountStatus).count() >= fileStatusInterval) {
                saveStatusBrute();
                beginCountStatus = std::chrono::steady_clock::now();
            }
        }

        cudaStatus = cudaDeviceSynchronize();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching kernel!\n", cudaStatus);
            goto Error;
        }

        if (buffIsResult[0]) {
            buffIsResult[0] = false;
            for (uint32_t i = 0; i < outputSizeB_T; i++) {
                if (buffDeviceResult[i] != NO_RESULT) {
                    //printf("%d %lld %u\n", c, SUFFIX, buffDeviceResult[i]); continue;
                    uint64_t base = SUFFIX;
                    base += (THREAD_STEPS_BRUTE * i) + buffDeviceResult[i];
                    string source = "";
                    char buf[64] = { 0 };
                    int bufStart = _lltoa(base, buf);
                    if (IS_ROOT) {
                        source += ROOT;
                        while (bufStart < 64) {
                            source += buf[bufStart++];
                        }
                    }
                    else {
                        while (bufStart < 64) {
                            source += buf[bufStart++];
                        }
                        source += ROOT_SUFFIX;
                    }
                    for (int turn = 1; turn <= 1 + (IS_CAMP2 ? 1 : 0); turn++) {
                        Int privKey;
                        privKey.SetInt32(0);
                        int32_t rewrittenKey[NB32BLOCK];
                        if (IS_ETH_SHA256) {
                            sha256((unsigned char*)source.c_str(), source.size(), (unsigned char*)privKey.bits64);
                            for (int rk = 0; rk < 8; rk++) {
                                rewrittenKey[rk] = (privKey.bits[rk] & 0x000000ff) << 24;
                                rewrittenKey[rk] |= (privKey.bits[rk] & 0x0000ff00) << 8;
                                rewrittenKey[rk] |= (privKey.bits[rk] & 0x00ff0000) >> 8;
                                rewrittenKey[rk] |= (privKey.bits[rk] & 0xff000000) >> 24;
                                privKey.bits[rk] = rewrittenKey[rk];
                            }
                        }
                        else {
                            unsigned char d_hash[32];
                            _keccak(source.c_str(), source.size(), d_hash, 32);
                            if (turn == 2) {
                                char toHash[32];
#pragma unroll 2030
                                for (int camp2Turns = 0; camp2Turns < 2030; camp2Turns++) {
#pragma unroll 32
                                    for (int xxx = 0; xxx < 32; xxx++) {
                                        toHash[xxx] = d_hash[xxx];
                                    }
                                    _keccak(toHash, 32, d_hash, 32);
                                }
                            }
                            for (int rk = 0; rk < 8; rk++) {
                                rewrittenKey[rk] = ((int32_t)d_hash[rk * 4]) << 24;
                                rewrittenKey[rk] |= (int32_t)(d_hash[rk * 4 + 1]) << 16;
                                rewrittenKey[rk] |= ((int32_t)d_hash[rk * 4 + 2]) << 8;
                                rewrittenKey[rk] |= ((int32_t)d_hash[rk * 4 + 3]);
                                privKey.bits[rk] = rewrittenKey[rk];
                            }
                        }
                        privKey.SetInt32(0);
                        for (int rk = 0; rk < 8; rk++) {
                            privKey.bits[7 - rk] = rewrittenKey[rk];
                        }
                        if (processCandidateEth(privKey, source, turn == 2)) {
                            printf("found: %s \n", source.c_str());
                            IS_RESULT = true;
                        }
                    }
                    buffDeviceResult[i] = NO_RESULT;
                }
            }
        }
        counter += workSize;
        SUFFIX += workSize;
    }


Error:
    cudaFree(_dev_precomp);
    return cudaStatus;
}

//=================================


/**
Returns the optimal bloom filter size in bits given the probability of false-positives and the
number of hash functions
*/
uint32_t CudaHashLookup::getOptimalBloomFilterBits(double p, size_t n)
{
    double m = 3.6 * ceil((n * log(p)) / log(1 / pow(2, log(2))));

    return (uint32_t)ceil(log(m) / log(2));
}

void CudaHashLookup::initializeBloomFilter64(const std::vector<struct hash160>& targets, uint32_t* filter, uint64_t mask)
{
    for (uint32_t k = 0; k < targets.size(); k++) {

        uint32_t hash[5];

        uint64_t idx[5];

        //undoRMD160FinalRound(targets[k].h, hash);
        for (int i = 0; i < 5; i++) {
            hash[i] = targets[k].h[i];
        }

        idx[0] = ((uint64_t)hash[0] << 32 | hash[1]) & mask;
        idx[1] = ((uint64_t)hash[2] << 32 | hash[3]) & mask;
        idx[2] = ((uint64_t)(hash[0] ^ hash[1]) << 32 | (hash[1] ^ hash[2])) & mask;
        idx[3] = ((uint64_t)(hash[2] ^ hash[3]) << 32 | (hash[3] ^ hash[4])) & mask;
        idx[4] = ((uint64_t)(hash[0] ^ hash[3]) << 32 | (hash[1] ^ hash[3])) & mask;

        for (int i = 0; i < 5; i++) {

            filter[idx[i] / 32] |= (0x01 << (idx[i] % 32));
        }
    }
}

/**
 * Populates the bloom filter with the target hashes.
 */
cudaError_t CudaHashLookup::setTargetBloomFilter(const std::vector<struct hash160>& targets, uint32_t bbits)
{
    uint32_t bloomFilterBits = bbits > 0 ? bbits : getOptimalBloomFilterBits(1.0e-9, targets.size());
    uint64_t bloomFilterSizeWords = (uint64_t)1 << (bloomFilterBits - 5);
    uint64_t bloomFilterBytes = (uint64_t)1 << (bloomFilterBits - 3);
    uint64_t bloomFilterMask = (((uint64_t)1 << bloomFilterBits) - 1);

    uint32_t* filter = NULL;

    fprintf(stderr, "Initializing Bloom filter [bits=%lu,size=%lu].\n",
        (unsigned long)bloomFilterBits, (unsigned long)bloomFilterSizeWords);
    try {
        filter = new uint32_t[bloomFilterSizeWords];
    }
    catch (std::bad_alloc) {
        return cudaErrorMemoryAllocation;
    }

    cudaError_t err = cudaMalloc(&_bloomFilterPtr, bloomFilterBytes);

    if (err) {
        delete[] filter;
        return err;
    }

    memset(filter, 0, sizeof(uint32_t) * bloomFilterSizeWords);

    initializeBloomFilter64(targets, filter, bloomFilterMask);
    /*
    if(bloomFilterBits > 32) {
        initializeBloomFilter64(targets, filter, bloomFilterMask);
    } else {
        initializeBloomFilter(targets, filter, (unsigned int)bloomFilterMask);
    }*/

    // Copy to device
    err = cudaMemcpy(_bloomFilterPtr, filter, sizeof(uint32_t) * bloomFilterSizeWords, cudaMemcpyHostToDevice);
    if (err) {
        printf("cudahashlookup settsrgetbloomfilters: cudamemcpy bloomfilterptr error!");
        cudaFree(_bloomFilterPtr);
        _bloomFilterPtr = NULL;
        delete[] filter;
        return err;
    }

    // Copy device memory pointer to constant memory
    //err = cudaMemcpyToSymbol(_BLOOM_FILTER, &_bloomFilterPtr, sizeof(uint32_t*));
    err = cudaMemcpyToSymbol_BLOOM_FILTER(_bloomFilterPtr);
    if (err) {
        printf("cudahashlookup settsrgetbloomfilters: cudamemcpytosymbol bloom filter error!");
        cudaFree(_bloomFilterPtr);
        _bloomFilterPtr = NULL;
        delete[] filter;
        return err;
    }

    // Copy device memory pointer to constant memory
    /*if (bloomFilterBits <= 32) {
        err = cudaMemcpyToSymbol(_BLOOM_FILTER_MASK, &bloomFilterMask, sizeof(unsigned int *));
        if(err) {
            printf("cudahashlookup settargetbloomfilters: cudamemcopytosymbol bloomfiltermask error!");
            cudaFree(_bloomFilterPtr);
            _bloomFilterPtr = NULL;
            delete[] filter;
            return err;
        }
    } else {*/
    //err = cudaMemcpyToSymbol(_BLOOM_FILTER_MASK64, &bloomFilterMask, sizeof(uint64_t*));
    err = cudaMemcpyToSymbol_BLOOM_FILTER_MASK64(&bloomFilterMask);
    if (err) {
        printf("cudahashlookup settsrgetbloomfilters: cudamemcopytosymbol bloomfiltermask64 error!");
        cudaFree(_bloomFilterPtr);
        _bloomFilterPtr = NULL;
        delete[] filter;
        return err;
    }
    //}

    //unsigned int useBloomFilter = bloomFilterBits <= 32 ? 1 : 2;
    uint32_t useBloomFilter = 2;

    //err = cudaMemcpyToSymbol(_USE_BLOOM_FILTER, &useBloomFilter, sizeof(uint32_t));
    err = cudaMemcpyToSymbol_USE_BLOOM_FILTER(&useBloomFilter);
    if (err) {
        printf("cudahashlookup settsrgetbloomfilters: cudamemcopytosymbol useBloomFilter error!");
        cudaFree(_bloomFilterPtr);
        _bloomFilterPtr = NULL;
        delete[] filter;
        return err;
    }
    delete[] filter;
    return err;
}

/**
*Copies the target hashes to either constant memory, or the bloom filter depending
on how many targets there are
*/
cudaError_t CudaHashLookup::setTargets(const std::vector<struct hash160>& targets, uint32_t bbits)
{
    cleanup();

    return setTargetBloomFilter(targets, bbits);
}

void CudaHashLookup::cleanup()
{
    if (_bloomFilterPtr != NULL) {
        cudaFree(_bloomFilterPtr);
        _bloomFilterPtr = NULL;
    }
}
