// Test: C++ version of vector addition
// Copyright 2024 LavaLab @ University of Virginia. All rights reserved.

#include "libpimsim.h"
#include <iostream>
#include <vector>
#include <getopt.h>
#include <stdint.h>
#include <iomanip>
#include <omp.h>
#include "../util.h"

using namespace std;

// Params ---------------------------------------------------------------------
typedef struct Params
{
    uint64_t vectorLength;
    char* configFile; 
    char* inputFile;
} Params;

void usage()
{
    fprintf(stderr,
            "\nUsage:  ./add [options]"
            "\n"
            "\n    -l    input size (default=8M elements)"
            "\n    -c    dramsim config file"
            "\n    -i    input file containing two vectors (default=generates vector with random numbers)"
            "\n");
}

struct Params getInputParams(int argc, char **argv)
{
    struct Params p;
    p.vectorLength = 65536;
    p.configFile = nullptr;
    p.inputFile = nullptr;

    int opt;
    while ((opt = getopt(argc, argv, "l:c:i:")) >= 0)
    {
        switch (opt)
        {
        case 'h':
            usage();
            exit(0);
            break;
        case 'l':
            p.vectorLength = strtoull(optarg,  NULL, 0);
            break;
        case 'c':
            p.configFile = optarg;
            break;
        case 'i':
            p.inputFile = optarg;
            break;
        default:
            fprintf(stderr, "\nUnrecognized option!\n");
            usage();
            exit(0);
        }
    }
    return p;
}

void vectorAddition(uint64_t vectorLength, std::vector<int>& src1, std::vector<int>& src2, std::vector<int>& dst, char* configFile)
{
  if (configFile == nullptr) {
    unsigned numCores = 48;
    unsigned numRows = 8192;
    unsigned numCols = 8192;
    PimStatus status = pimCreateDevice(PIM_FUNCTIONAL, numCores, numRows, numCols);
    if (status != PIM_OK)
    {
      std::cout << "Abort" << std::endl;
      return ;
    }
  } else {
    PimStatus status = pimCreateDeviceFromConfig(PIM_FUNCTIONAL, configFile);
    if (status != PIM_OK) {
      std::cout << "Abort" << std::endl;
      return ;
    }
  }

  unsigned bitsPerElement = sizeof(int)*8;
  PimObjId srcObj1 = pimAlloc(PIM_ALLOC_V1, vectorLength, bitsPerElement, PIM_INT32);
  if (srcObj1 == -1) {
    std::cout << "Abort" << std::endl;
    return ;
  }
  PimObjId srcObj2 = pimAllocAssociated(PIM_ALLOC_V1, vectorLength, bitsPerElement, srcObj1, PIM_INT32);
  if (srcObj2 == -1) {
    std::cout << "Abort" << std::endl;
    return ;
  }
  PimObjId dstObj = pimAllocAssociated(PIM_ALLOC_V1, vectorLength, bitsPerElement, srcObj1, PIM_INT32);
  if (dstObj == -1) {
    std::cout << "Abort" << std::endl;
    return ;
  }

  PimStatus status = pimCopyHostToDevice(PIM_COPY_V, (void*)src1.data(), srcObj1);
  if (status != PIM_OK) {
    std::cout << "Abort" << std::endl;
    return ;
  }

  status = pimCopyHostToDevice(PIM_COPY_V, (void*)src2.data(), srcObj2);
  if (status != PIM_OK) {
    std::cout << "Abort" << std::endl;
    return ;
  }

  status = pimAdd(srcObj1, srcObj2, dstObj);
  if (status != PIM_OK) {
    std::cout << "Abort" << std::endl;
    return ;
  }

  dst.reserve(vectorLength);
  status = pimCopyDeviceToHost(PIM_COPY_V, dstObj, (void*)dst.data());
  if (status != PIM_OK) {
    std::cout << "Abort" << std::endl;
    return ;
  }

  pimFree(srcObj1);
  pimFree(srcObj2);
  pimFree(dstObj);
  //verify result
  #pragma omp parallel for
  for (unsigned i = 0; i < vectorLength; ++i) {
    int sum = src1[i] + src2[i];
    if (dst[i] != sum) {
      std::cout << "Wrong answer for addition: " << src1[i] << " + " << src2[i] << " = " << dst[i] << " (expected " << sum << ")" << std::endl;
    }
  }
}


int main(int argc, char *argv[])
{
  struct Params params = getInputParams(argc, argv);
  std::cout << "Vector length: " << params.vectorLength << "\n";
  std::vector<int> src1, src2, dst;
  if (params.inputFile == nullptr) {
    getVector(params.vectorLength, src1);
    getVector(params.vectorLength, src2);
  }

  vectorAddition(params.vectorLength, src1, src2, dst, params.configFile);

  pimShowStats();
  

  return 0;
}
