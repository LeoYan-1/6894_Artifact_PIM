#include <iostream>
#include <omp.h>
#include <vector>
using namespace std;

#ifndef _UTIL_H_
#define _UTIL_H_

#ifndef DATA_TYPE
typedef int32_t data_t;
#else
typedef DATA_TYPE data_t;
#endif

void getVector(uint64_t vectorLength, std::vector<int>& srcVector) {
  srand((unsigned)time(NULL));
  srcVector.reserve(vectorLength);
  #pragma omp parallel for
  for (int i = 0; i < vectorLength; ++i)
  {
    srcVector[i] = rand() % (i+1);
  }
}

void getMatrix(int row, int column, int padding, std::vector<std::vector<int>>& inputMatrix) {
    srand((unsigned)time(NULL));
    inputMatrix.resize(row + 2 * padding, std::vector<int>(column + 2 * padding, 0));
    #pragma omp parallel for
    for (int i = padding; i < row + padding; ++i) {
        for (int j = padding; j < column + padding; ++j) {
            inputMatrix[i][j] = rand() % (i * j + 1);
        }
    }
}

#endif
