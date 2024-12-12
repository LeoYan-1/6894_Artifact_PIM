// Triangle Counting Benchmark
// Copyright (c) 2024 University of Virginia
// This file is licensed under the MIT License.
// See the LICENSE file in the root of this repository for more details.

#include <iostream>
#include <vector>
#include <cassert>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <bitset>
#include <unordered_set>
#include <sstream>
#include <getopt.h>
#include <unordered_map>
#if defined(_OPENMP)
#include <omp.h>
#endif

#include "../../util.h"
#include "libpimeval.h"

#define DEBUG 0
#define BITS_PER_INT 32

#define ADJ_MATRIX 0
#define ADJ_LIST 1
// #define INPUT_MODE ADJ_LIST
#define INPUT_MODE ADJ_MATRIX

#define USE_OPT 0

typedef uint32_t UINT32;


using namespace std;


// Params ---------------------------------------------------------------------
typedef struct Params
{
  const char *configFile;
  const char *inputFile;
  bool shouldVerify;
  std::string algorithm;
  UINT32 tileSize = 10;
  UINT32 sliceNum = 1;
  

} Params;

void usage()
{
  fprintf(stderr,
          "\nUsage:  ./tc.out [options]"
          "\n"
          "\n    -c    dramsim config file"
          "\n    -i    input file containing two vectors (default=generates vector with random numbers)"
          "\n    -v    t = verifies PIM output with host output. (default=false)"
          "\n    -a    Algorithms. Choose from \"BASELINE\", \"SERIAL\", \"TILED\", or \"SLICED\""
          "\n    -t    tile size (for TILED algorithm only)"
          "\n    -s    number of slices (for SLICED algorithm only)"
          "\n");
}

struct Params getInputParams(int argc, char **argv)
{
  struct Params p;
  p.configFile = nullptr;
  p.inputFile = "Dataset/v18772_symetric";
  p.shouldVerify = false;
  p.algorithm = "BASELINE";

  int opt;
  while ((opt = getopt(argc, argv, "h:c:i:v:a:t:s:q:")) >= 0)
  {
    switch (opt)
    {
    case 'h':
      usage();
      exit(0);
      break;
    case 'c':
      p.configFile = optarg;
      break;
    case 'i':
      p.inputFile = optarg;
      break;
    case 'v':
      p.shouldVerify = (*optarg == 't') ? true : false;
      break;
    case 'a':
      static const std::unordered_set<std::string> validAlgorithms = {"BASELINE", "SERIAL", "TILED", "SLICED"};
      if (validAlgorithms.find(std::string(optarg)) == validAlgorithms.end()){
        fprintf(stderr, "\nInvalid algorithm specified. Choose from \"BASELINE\", \"SERIAL\", \"TILED\", or \"SLICED\".\n");
      } else {
        p.algorithm = optarg;
      };
      break;
    case 't':
      p.tileSize = atoi(optarg);
      break;
    case 's':
      p.sliceNum = atoi(optarg);
      break;
    case 'q':
      break;
    default:
      fprintf(stderr, "\nUnrecognized option!\n");
      usage();
      exit(0);
    }
  }
  return p;
}

// Function to convert edge list to adjacency matrix
vector<vector<bool>> edgeListToAdjMatrix(const vector<pair<int, int>>& edgeList, int numNodes) {
    vector<vector<bool>> adjMatrix(numNodes+1, vector<bool>(numNodes+1, 0));

    for (const auto& edge : edgeList) {
        int u = edge.first;
        int v = edge.second;
        adjMatrix[u][v] = adjMatrix[v][u] = 1; // assuming undirected graph
    }

    return adjMatrix;
}

// Function to convert standard adjacency matrix to bitwise adjacency matrix
vector<vector<UINT32>> convertToBitwiseAdjMatrix(const vector<vector<bool>>& adjMatrix) {
    int V = adjMatrix.size();
    int numInts = (V + BITS_PER_INT - 1) / BITS_PER_INT; // Number of 32-bit integers needed per row

    vector<vector<UINT32>> bitAdjMatrix(V, vector<UINT32>(numInts, 0));
    int step = V / 10; // Each 10 percent of the total iterations

    for (int i = 0; i < V; ++i) {
        for (int j = 0; j < V; ++j) {
            if (adjMatrix[i][j]) {
                bitAdjMatrix[i][j / BITS_PER_INT] |= (1 << (j % BITS_PER_INT));
            }
        }
        if (i % step == 0) {
            std::cout << "convertToBitwiseAdjMatrix: Progress: " << (i * 100 / V) << "\% completed." << std::endl;
        }
    }

    return bitAdjMatrix;
}

// Function to read edge list from a JSON file
vector<pair<int, int>> readEdgeList(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Unable to open file");
    }

    vector<pair<int, int>> edgeList;

    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        int u, v;
        if (!(iss >> u >> v)) {
            throw runtime_error("Invalid file format");
        }
        edgeList.emplace_back(u, v);
    }
    cout << "Edge list size: " << edgeList.size() << endl;
    return edgeList;
}

int vectorAndPopCntRedSum(uint64_t numElements, std::vector<unsigned int> &src1, std::vector<unsigned int> &src2, std::vector<unsigned int> &dst, std::vector<unsigned int> &popCountSrc) {
    cout << "numElements: " << numElements << endl;

    PimObjId srcObj1 = pimAlloc(PIM_ALLOC_AUTO, numElements, PIM_INT32);
    if (srcObj1 == -1)
    {
        std::cout << "src1: pimAlloc" << std::endl;
        return -1;
    }
    if (DEBUG) cout << "src1 allocated successfully!" << endl;

    PimObjId srcObj2 = pimAllocAssociated(srcObj1, PIM_INT32);
    if (srcObj2 == -1)
    {
        std::cout << "src2: pimAllocAssociated" << std::endl;
        return -1;
    }
    if (DEBUG) cout << "src2 allocated successfully!" << endl;

    PimObjId dstObj = pimAllocAssociated(srcObj1, PIM_INT32);
    if (dstObj == -1)
    {
        std::cout << "dst: pimAllocAssociated" << std::endl;
        return -1;
    }
    if (DEBUG) cout << "dst allocated successfully!" << endl;

    PimObjId popCountSrcObj = pimAllocAssociated(srcObj1, PIM_INT32);
    if (popCountSrcObj == -1)
    {
        std::cout << "popCountSrc: pimAllocAssociated" << std::endl;
        return -1;
    }
    if (DEBUG) cout << "popCountSrc allocated successfully!" << endl;

    PimStatus status = pimCopyHostToDevice((void *)src1.data(), srcObj1);
    if (status != PIM_OK)
    {
        std::cout << "src1: pimCopyHostToDevice Abort" << std::endl;
        return -1;
    }
    if (DEBUG) cout << "src1 copied successfully!" << endl;

    status = pimCopyHostToDevice((void *)src2.data(), srcObj2);
    if (status != PIM_OK)
    {
        std::cout << "src2: pimCopyHostToDevice Abort" << std::endl;
        return -1;
    }
    if (DEBUG) cout << "src2 copied successfully!" << endl;

    status = pimAnd(srcObj1, srcObj2, dstObj);
    if (status != PIM_OK)
    {
        std::cout << "pimAnd Abort" << std::endl;
        return -1;
    }
    if (DEBUG) cout << "pimAnd completed successfully!" << endl;

    status = pimPopCount(dstObj, popCountSrcObj);
    if (status != PIM_OK)
    {
        std::cout << "pimPopCount Abort" << std::endl;
        return -1;
    }
    if (DEBUG) cout << "pimPopCount completed successfully!" << endl;

    int64_t sum = 0;
    status = pimRedSumInt(popCountSrcObj, &sum);
    if (status != PIM_OK)
    {
        std::cout << "pimRedSumInt Abort" << std::endl;
        return -1;
    }
    if (DEBUG) cout << "pimRedSum completed successfully!" << endl;

    pimFree(srcObj1);
    pimFree(srcObj2);
    pimFree(dstObj);
    pimFree(popCountSrcObj);

    return sum;
}

int run_rowmaxusage_opt(const vector<vector<bool>>& adjMatrix, const vector<vector<UINT32>>& bitAdjMatrix, uint64_t words_per_device) {
    uint64_t operandsCount = 4; // src1, src2, dst, popCountSrc
    uint64_t operandMaxNumberOfWords = words_per_device / operandsCount;
    int count = 0;
    int V = bitAdjMatrix.size();
    uint64_t wordsPerMatrixRow = (V + BITS_PER_INT - 1) / BITS_PER_INT; // Number of 32-bit integers needed per row
    cout << "wordsPerMatrixRow: " << wordsPerMatrixRow << endl;
    cout << "words_per_device: " << words_per_device << endl;
    assert(wordsPerMatrixRow <=  operandMaxNumberOfWords && "Number of vertices cannot exceed (words_per_device / 2)");
    int oneCount = 0;
    uint64_t words = 0;
    std::vector<unsigned int> src1;
    std::vector<unsigned int> src2;
    int step = V / 10; // Each 10 percent of the total iterations
    uint16_t iterations = 0;
    double host_time_if = 0.0, host_time_forloop = 0.0;
    for (int i = 0; i < V; ++i) {
        for (int j = 0; j < V; ++j) {
            if (adjMatrix[i][j]) 
            { // If there's an edge between i and j
                ++oneCount;
                auto start = std::chrono::high_resolution_clock::now();
#ifdef ENABLE_PARALLEL
                // Parallelizing the loop with OpenMP
                #pragma omp parallel for reduction(+:host_time_if, words) 
#endif
                for (uint64_t k = 0; k < wordsPerMatrixRow; ++k) {
                    unsigned int op1 = bitAdjMatrix[i][k];
                    unsigned int op2 = bitAdjMatrix[j][k];
                    auto ifstart = std::chrono::high_resolution_clock::now();
                    if (op1 == 0 || op2 == 0)
                    {
                        auto ifend = std::chrono::high_resolution_clock::now();
                        auto ifelapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(ifend - ifstart);
                        host_time_if += ifelapsedTime.count();
                        continue;
                    }
                    auto ifend = std::chrono::high_resolution_clock::now();
                    auto ifelapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(ifend - ifstart);
                    host_time_if += ifelapsedTime.count();
#ifdef ENABLE_PARALLEL
                    #pragma omp atomic
#endif
                    ++words;
#ifdef ENABLE_PARALLEL
                    #pragma omp critical
#endif
                    {
                        src1.push_back(op1);
                        src2.push_back(op2);
                    }
                }
                auto end = std::chrono::high_resolution_clock::now();
                auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
                host_time_forloop += elapsedTime.count();
            }
            if((words + wordsPerMatrixRow > operandMaxNumberOfWords) || ((i == V - 1) && (j == V - 1) && words > 0)){
                cout << "words: " << words << endl;
                cout << "-------------itr[" << iterations << "]-------------" << endl;
                std::vector<unsigned int> dst(words);
                std::vector<unsigned int> popCountSrc(words);
                int sum = vectorAndPopCntRedSum((uint64_t) words, src1, src2, dst, popCountSrc);
           
                if(sum < 0)
                    return -1;
                words = 0;
                iterations++;
                src1.clear();
                src2.clear();
                dst.clear();
                count += sum;
            }

        }
        if (i % step == 0) {
            std::cout << "run_adjmatrix: Progress: " << (i * 100 / V) << "\% rows completed." << std::endl;
        }
    }
    cout << "Host time (for loop): " << std::fixed << std::setprecision(3) << host_time_forloop << " ms." << endl;
    cout << "Host time (if): " << std::fixed << std::setprecision(3) << host_time_if << " ms." << endl;
    cout << "TriangleCount: " << count / 6 << endl;
    // Each triangle is counted 6 times (once at each vertex), so divide the count by 6
    return count / 6;
}

int run_adjmatrix(const vector<vector<bool>>& adjMatrix, const vector<vector<UINT32>>& bitAdjMatrix, uint64_t words_per_device) {
    uint64_t operandsCount = 4; // src1, src2, dst, popCountSrc
    uint64_t operandMaxNumberOfWords = words_per_device / operandsCount;
    int count = 0;
    int V = bitAdjMatrix.size();
    uint64_t wordsPerMatrixRow = (V + BITS_PER_INT - 1) / BITS_PER_INT; // Number of 32-bit integers needed per row
    assert(wordsPerMatrixRow <=  operandMaxNumberOfWords && "Number of vertices cannot exceed (words_per_device / 2)");
    int oneCount = 0;
    uint64_t words = 0;
    std::vector<unsigned int> src1;
    std::vector<unsigned int> src2;
    int step = V / 10; // Each 10 percent of the total iterations
    uint16_t iterations = 0;

    for (int i = 0; i < V; ++i) {
        for (int j = 0; j < V; ++j) {
            if (adjMatrix[i][j]) 
            { // If there's an edge between i and j
                ++oneCount;
                for (uint64_t k = 0; k < wordsPerMatrixRow; ++k) {
                    unsigned int op1 = bitAdjMatrix[i][k];
                    unsigned int op2 = bitAdjMatrix[j][k];
                    ++words;
                    src1.push_back(op1);
                    src2.push_back(op2);
                }
            }
            if((words + wordsPerMatrixRow > operandMaxNumberOfWords) || ((i == V - 1) && (j == V - 1) && words > 0)){
                cout << "-------------itr[" << iterations << "]-------------" << endl;
                cout << "Number of words that are processed in this iteration: " << words << endl;
                std::vector<unsigned int> dst(words);
                std::vector<unsigned int> popCountSrc(words);
                int sum = vectorAndPopCntRedSum((uint64_t) words, src1, src2, dst, popCountSrc);

                if(sum < 0)
                    return -1;
                
                words = 0;
                iterations++;
                src1.clear();
                src2.clear();
                dst.clear();
                count += sum;
            }

        }
        if (i % step == 0) {
            std::cout << "Progress: " << (i * 100 / V) << "\% rows completed." << std::endl;
        }
    }
    // Each triangle is counted 6 times (once at each vertex), so divide the count by 6
    return count / 6;
}

#include <algorithm>
#include <queue>

//record the allocated PIM colum for data reuse
class PIMColListElem{
public:    
    int         ColIdx;
    PimObjId    ObjId;

    PIMColListElem(int colidx, PimObjId objid): ColIdx(colidx), ObjId(objid){}
    ~PIMColListElem(){}
};



int run_adjmatrix_serial(const vector<vector<bool>>& adjMatrix, const vector<vector<UINT32>>& bitAdjMatrix, uint64_t words_per_device) {
    uint64_t operandsCount = 4; // src1, src2, dst, popCountSrc
    uint64_t operandMaxNumberOfWords = words_per_device / operandsCount;
    int count = 0;
    int V = bitAdjMatrix.size();
    uint64_t wordsPerMatrixRow = (V + BITS_PER_INT - 1) / BITS_PER_INT; // Number of 32-bit integers needed per row
    assert(wordsPerMatrixRow <=  operandMaxNumberOfWords && "Number of vertices cannot exceed (words_per_device / 2)");
    int oneCount = 0;
    uint64_t words = 0;
    std::vector<unsigned int> src1;
    std::vector<unsigned int> src2;
    int step = V / 10; // Each 10 percent of the total iterations
    uint16_t iterations = 0;

    std::deque <PIMColListElem> pimColList;
    //allocate row
    PimObjId srcObj0 = pimAlloc(PIM_ALLOC_AUTO, wordsPerMatrixRow, PIM_INT32);
    if (srcObj0 == -1)
    {
        std::cout << "src0: pimAlloc" << std::endl;
        return -2;
    }
    
    PimObjId dstObj = pimAllocAssociated(srcObj0, PIM_INT32);
    if (dstObj == -1)
    {
        std::cout << "dst: pimAllocAssociated" << std::endl;
        return -1;
    }
    if (DEBUG) cout << "dst allocated successfully!" << endl;

    PimObjId popCountSrcObj = pimAllocAssociated(srcObj0, PIM_INT32);
    if (popCountSrcObj == -1)
    {
        std::cout << "popCountSrc: pimAllocAssociated" << std::endl;
        return -1;
    }
    if (DEBUG) cout << "popCountSrc allocated successfully!" << endl;


    for (int i = 0; i < V; ++i) {
        PimStatus status = pimCopyHostToDevice((void *)bitAdjMatrix[i].data(), srcObj0);
        if (status != PIM_OK)
        {
            std::cout << "src1: pimCopyHostToDevice Abort" << std::endl;
            return -1;
        }
        if (DEBUG) cout << "src1 copied successfully!" << endl;
        for (int k = 0; k < V; ++k) {
            int j = (i%2)? (V-1-k) : k;
            if (adjMatrix[i][j]) 
            { // If there's an edge between i and j
                ++oneCount;
                // cout << "\nnow in point [" << i << "][" << j << "]" << std::endl;


                auto it = std::find_if(pimColList.begin(), pimColList.end(), [j](const PIMColListElem elem){
                    return elem.ColIdx == j;
                });
                if(it != pimColList.end())
                {
                    // std::cout << "reuse column " << j << std::endl;
                    status = pimAnd(srcObj0, it->ObjId, dstObj);
                    if (status != PIM_OK)
                    {
                        std::cout << "pimAnd Abort" << std::endl;
                        return -1;
                    }
                    if (DEBUG) cout << "pimAnd completed successfully!" << endl;

                    status = pimPopCount(dstObj, popCountSrcObj);
                    if (status != PIM_OK)
                    {
                        std::cout << "pimPopCount Abort" << std::endl;
                        return -1;
                    }
                    if (DEBUG) cout << "pimPopCount completed successfully!" << endl;

                    int64_t sum = 0;
                    status = pimRedSumInt(popCountSrcObj, &sum);
                    if (status != PIM_OK)
                    {
                        std::cout << "pimRedSumInt Abort" << std::endl;
                        return -1;
                    }
                    if (DEBUG) cout << "pimRedSum completed successfully!" << endl; 
                    // cout << "counted: " << sum << " triangles" << std::endl;
                    count += sum;
                }
                else
                {
                    PimObjId srcObjNew = pimAllocAssociated(srcObj0, PIM_INT32);
                    if (srcObjNew == -1)
                    {
                        std::cout << "srcNew: pimAllocAssociated" << pimColList.size() << std::endl;
                        return -1;
                    }
                    if (DEBUG) cout << "srcNew allocated successfully!" << endl;
                    if (pimColList.size()>=60){
                        pimFree(pimColList[0].ObjId);
                        pimColList.pop_front();
                    }
                    PIMColListElem newCol=PIMColListElem(j, srcObjNew);
                    pimColList.push_back(newCol);//add a new column to PIM
                    // std::cout << "add colum "<< j << " to list" << std::endl;

                    //copy new column to PIM
                    // cout << "copying data: " << *(bitAdjMatrix[j].data()) << " to device" << std::endl;
                    PimStatus status = pimCopyHostToDevice((void *)bitAdjMatrix[j].data(), srcObjNew);
                    if (status != PIM_OK)
                    {
                        std::cout << "srcNew: pimCopyHostToDevice Abort" << std::endl;
                        return -1;
                    }
                    if (DEBUG) cout << "srcNew copied successfully!" << endl;

                    status = pimAnd(srcObj0, srcObjNew, dstObj);
                    if (status != PIM_OK)
                    {
                        std::cout << "pimAnd Abort" << std::endl;
                        return -1;
                    }
                    if (DEBUG) cout << "pimAnd completed successfully!" << endl;

                    status = pimPopCount(dstObj, popCountSrcObj);
                    if (status != PIM_OK)
                    {
                        std::cout << "pimPopCount Abort" << std::endl;
                        return -1;
                    }
                    if (DEBUG) cout << "pimPopCount completed successfully!" << endl;

                    int64_t sum = 0;
                    status = pimRedSumInt(popCountSrcObj, &sum);
                    if (status != PIM_OK)
                    {
                        std::cout << "pimRedSumInt Abort" << std::endl;
                        return -1;
                    }
                    if (DEBUG) cout << "pimRedSum completed successfully!" << endl; 
                    // cout << "counted: " << sum << " triangles" << std::endl;
                    count += sum;                    
                }

            }


        }
        if (i % step == 0) {
            std::cout << "Progress: " << (i * 100 / V) << "\% rows completed." << std::endl;
        }
    }
    // Each triangle is counted 6 times (once at each vertex), so divide the count by 6
    pimFree(srcObj0);
    for(size_t i=0; i<pimColList.size(); i++){
        pimFree(pimColList[i].ObjId);
    }
    pimFree(dstObj);
    pimFree(popCountSrcObj);
    return count / 6;


}

int run_adjmatrix_tiled(const vector<vector<bool>>& adjMatrix, const vector<vector<UINT32>>& bitAdjMatrix, uint64_t words_per_device, int tile_size) {
    uint64_t operandsCount = 4; // src1, src2, dst, popCountSrc
    uint64_t operandMaxNumberOfWords = words_per_device / operandsCount;
    int count = 0;
    int V = bitAdjMatrix.size();

    uint64_t wordsPerMatrixRow = (V + BITS_PER_INT - 1) / BITS_PER_INT; // Number of 32-bit integers needed per row
    assert(wordsPerMatrixRow <=  operandMaxNumberOfWords && "Number of vertices cannot exceed (words_per_device / 2)");
    int oneCount = 0;
    uint64_t words = 0;
    std::vector<unsigned int> src1;
    std::vector<unsigned int> src2;

    int tiles = (V + tile_size - 1) / tile_size;
    int batch_size = 1;

    int step = tile_size / 10; // Each 10 percent of the total iterations
    uint16_t iterations = 0;

    std::deque <PIMColListElem> pimColList;
    //allocate row

    if (wordsPerMatrixRow * tiles >= operandMaxNumberOfWords) return 1;
    PimObjId srcObj0 = pimAlloc(PIM_ALLOC_AUTO, wordsPerMatrixRow * tiles, PIM_INT32);
    if (srcObj0 == -1)
    {
        std::cout << "src0: pimAlloc" << std::endl;
        return -2;
    }
    
    PimObjId dstObj = pimAllocAssociated(srcObj0, PIM_INT32);
    if (dstObj == -1)
    {
        std::cout << "dst: pimAllocAssociated" << std::endl;
        return -1;
    }
    if (DEBUG) cout << "dst allocated successfully!" << endl;

    PimObjId popCountSrcObj = pimAllocAssociated(srcObj0, PIM_INT32);
    if (popCountSrcObj == -1)
    {
        std::cout << "popCountSrc: pimAllocAssociated" << std::endl;
        return -1;
    }
    if (DEBUG) cout << "popCountSrc allocated successfully!" << endl;

    //allocate column batch buffer
    vector<PimObjId> srcObj1;
    for (int k = 0; k < batch_size; ++k){
        srcObj1.push_back(pimAllocAssociated(srcObj0, PIM_INT32));
        if (srcObj1.back() == -1)
        {
            std::cout << "srcObj1: pimAllocAssociated" << pimColList.size() << std::endl;
            return -1;
        }
        if (DEBUG) cout << "srcObj1 tile" << k << " allocated successfully!" << endl;
    }


    vector<int> colPtr(tiles, 0);
    vector<int> rowPtr(tiles, 0);
    vector<bool> dsablMask(tiles, false);

    while (std::any_of(colPtr.begin(), colPtr.end(), [V](int val) {return val < V;}) || 
            std::any_of(rowPtr.begin(), rowPtr.end(), [tile_size](int val) {return val < tile_size;}))
    {
        std::vector<std::vector<UINT32>> nextBatch(tiles, std::vector<UINT32>());
        for (int t = 0; t < tiles; ++t){
            int current_row = t*tile_size + rowPtr[t];

            //copy row to device at the beginning
            if (colPtr[t] == 0){
                const vector<UINT32>& copy_data = (rowPtr[t] >= tile_size || current_row >= V)? 
                                                vector<UINT32>(bitAdjMatrix.begin()->size(), 0) : 
                                                bitAdjMatrix[current_row];
                uint64_t idx_begin = t * wordsPerMatrixRow;
                uint64_t idx_end = idx_begin + wordsPerMatrixRow;
                PimStatus status = pimCopyHostToDevice((void *)copy_data.data(), srcObj0, idx_begin, idx_end);
                if (status != PIM_OK)
                {
                    std::cout << "src0: pimCopyHostToDevice Abort" << std::endl;
                    return -1;
                }
                if (DEBUG) cout << "src0 copied successfully!" << endl;
            }


            //if row overflow, disable calculation
            if (rowPtr[t] >= tile_size || current_row >= V){
                rowPtr[t] = tile_size;
                colPtr[t] = V;
                dsablMask[t] = true;
                continue;
                cout << "end of tile" << endl;
            }


            //for each batch, collect active column index
            while(nextBatch[t].size() < batch_size) 
            {
                if (adjMatrix[current_row][colPtr[t]]){
                    nextBatch[t].push_back(colPtr[t]);
                }

                colPtr[t]++;
                
                if (colPtr[t] >= V){
                    colPtr[t] = 0;
                    rowPtr[t]++;
                    //fill the rest of the nextBatch list with "invalid" flag
                    nextBatch[t].insert(nextBatch[t].end(), batch_size-nextBatch[t].size(), V);
                }
            }
        }


        //
        for (int t = 0; t < tiles; ++t){
            if (!dsablMask[t]){
                for (int k = 0; k < batch_size; ++k){
                    const vector<UINT32>& copyCol = (nextBatch[t][k] >= V)? 
                                                        std::vector<UINT32>(bitAdjMatrix.begin()->size(), 0) : 
                                                        bitAdjMatrix[nextBatch[t][k]];

                    uint64_t idx_begin = t * wordsPerMatrixRow;
                    uint64_t idx_end = idx_begin + wordsPerMatrixRow;
                    PimStatus status = pimCopyHostToDevice((void *)copyCol.data(), srcObj1[k], idx_begin, idx_end);
                    if (status != PIM_OK)
                    {
                        std::cout << "srcObj" << k << ": pimCopyHostToDevice Abort" << std::endl;
                        return -1;
                    }
                    if (DEBUG) cout << "srcObj" << k << " copied successfully!" << endl;
                }
            }
        }

        for (auto& list : nextBatch){
            list.clear();
        }

        for (int k = 0; k < batch_size; ++k){
            PimStatus status = pimAnd(srcObj0, srcObj1[k], dstObj);
            if (status != PIM_OK)
            {
                std::cout << "pimAnd Abort" << std::endl;
                return -1;
            }
            if (DEBUG) cout << "pimAnd completed successfully!" << endl;

            status = pimPopCount(dstObj, popCountSrcObj);
            if (status != PIM_OK)
            {
                std::cout << "pimPopCount Abort" << std::endl;
                return -1;
            }
            if (DEBUG) cout << "pimPopCount completed successfully!" << endl;

            int64_t sum = 0;
            status = pimRedSumInt(popCountSrcObj, &sum);
            if (status != PIM_OK)
            {
                std::cout << "pimRedSumInt Abort" << std::endl;
                return -1;
            }
            if (DEBUG) cout << "pimRedSum completed successfully!" << endl; 
            count += sum;
        } 

    }


    // Each triangle is counted 6 times (once at each vertex), so divide the count by 6
    pimFree(srcObj0);
    for(auto& obj : srcObj1){
        pimFree(obj);
    }
    pimFree(dstObj);
    pimFree(popCountSrcObj);
    return count / 6;


}



int run_adjmatrix_sliced(const vector<vector<bool>>& adjMatrix, const vector<vector<UINT32>>& bitAdjMatrix, uint64_t words_per_device, int slices) {
    uint64_t operandsCount = 8; // src1, src2, dst, popCountSrc
    uint64_t operandMaxNumberOfWords = words_per_device / operandsCount;
    int count = 0;
    int V = bitAdjMatrix.size();
    uint64_t wordsPerMatrixRow = (V + BITS_PER_INT - 1) / BITS_PER_INT; // Number of 32-bit integers needed per row
    assert(wordsPerMatrixRow <=  operandMaxNumberOfWords && "Number of vertices cannot exceed (words_per_device / 2)");
    
    uint64_t sliceSize = (wordsPerMatrixRow + slices -1) / slices; // Size of each slice
    assert(sliceSize > 0 && "sliceSize must be greater than 0");

    uint64_t words = 0;
    std::vector<unsigned int> src1;
    std::vector<unsigned int> src2;
    int step = V / 10; // Each 10 percent of the total iterations
    uint16_t iterations = 0;

    for (int i = 0; i < V; ++i) {
        for (int j = 0; j < V; ++j) {
            if (adjMatrix[i][j]) { 
                std::vector<unsigned int> rowSlices; 
                std::vector<unsigned int> colSlices; 

                
                for (uint64_t slice = 0; slice * sliceSize < wordsPerMatrixRow; ++slice) {
                    uint64_t start = slice * sliceSize;
                    uint64_t end = std::min((slice + 1) * sliceSize, wordsPerMatrixRow);

                    
                    bool rowNonZero = false;
                    bool colNonZero = false;
                    for (uint64_t k = start; k < end; ++k) {
                        if (bitAdjMatrix[i][k] != 0) rowNonZero = true;
                        if (bitAdjMatrix[j][k] != 0) colNonZero = true;
                        if (rowNonZero && colNonZero) break; 
                    }

                    
                    if (rowNonZero && colNonZero) {
                        for (uint64_t k = start; k < end; ++k) {
                            rowSlices.push_back(bitAdjMatrix[i][k]);
                            colSlices.push_back(bitAdjMatrix[j][k]);
                        }
                    }
                }

               if (!rowSlices.empty() && !colSlices.empty()) {
                    src1.insert(src1.end(), rowSlices.begin(), rowSlices.end()); 
                    src2.insert(src2.end(), colSlices.begin(), colSlices.end()); 
                    words += rowSlices.size(); 
                }
            }

            if((words + wordsPerMatrixRow > operandMaxNumberOfWords) || ((i == V - 1) && (j == V - 1) && words > 0)){
                cout << "-------------itr[" << iterations << "]-------------" << endl;
                cout << "Number of words that are processed in this iteration: " << words << endl;
                std::vector<unsigned int> dst(words);
                std::vector<unsigned int> popCountSrc(words);
                int sum = vectorAndPopCntRedSum((uint64_t) words, src1, src2, dst, popCountSrc);

                if(sum < 0)
                    return -1;
                
                words = 0;
                iterations++;
                src1.clear();
                src2.clear();
                dst.clear();
                count += sum;
            }

        }
        if (i % step == 0) {
            std::cout << "Progress: " << (i * 100 / V) << "\% rows completed." << std::endl;
        }
    }
    // Each triangle is counted 6 times (once at each vertex), so divide the count by 6
    return count / 6;
}



// Function to convert the adjacency list of a node to bitmap format
void convertToBitMap(const unordered_set<int>& neighbors, int start, vector<uint32_t>& bitMap) {
    for (int neighbor : neighbors) {
        // Set the corresponding bit for each neighbor
        bitMap[start + (neighbor / BITS_PER_INT)] |= (1 << (neighbor % BITS_PER_INT));
    }
}

int run_adjlist(const unordered_map<int, unordered_set<int>>& adjList, uint64_t words_per_device) {
    uint64_t operandsCount = 4; // src1, src2, dst, popCountSrc
    uint64_t operandMaxNumberOfWords = words_per_device / operandsCount;
    cout << "operandMaxNumberOfWords: " << operandMaxNumberOfWords << endl;
    int count = 0;
    int V = adjList.size();
    uint64_t wordsPerMatrixRow = (V + BITS_PER_INT - 1) / BITS_PER_INT; // Number of 32-bit integers needed per row
    assert(wordsPerMatrixRow <=  operandMaxNumberOfWords && "Number of vertices cannot exceed (words_per_device / 2)");
    uint64_t words = 0;
    int step = V / 10; // Each 10 percent of the total iterations
    uint16_t iterations = 0;
    uint32_t i = 0;
    vector<uint32_t> src1(wordsPerMatrixRow, 0);
    vector<uint32_t> src2(wordsPerMatrixRow, 0);
    vector<uint32_t> opt_src1, opt_src2;
    double host_time_if = 0.0;
    for (auto it_u = adjList.begin(); it_u != adjList.end(); ++it_u) {
        const auto& neighborsU = it_u->second;
        for (auto it_v = neighborsU.begin(); it_v != neighborsU.end(); ++it_v) {
            const auto& neighborsV = adjList.find(*it_v)->second;
            convertToBitMap(neighborsU, 0, src1);
            convertToBitMap(neighborsV, 0, src2);
            for (uint64_t j = 0; j < wordsPerMatrixRow; ++j) {
                if(USE_OPT) {
                    auto ifstart = std::chrono::high_resolution_clock::now();
                    if(src1[j] == 0 || src2[j] == 0) {
                        auto ifend = std::chrono::high_resolution_clock::now();
                        auto ifelapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(ifend - ifstart);
                        host_time_if += ifelapsedTime.count();
                        continue;
                    }
                    auto ifend = std::chrono::high_resolution_clock::now();
                    auto ifelapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(ifend - ifstart);
                    host_time_if += ifelapsedTime.count();
                }
                ++words;
                opt_src1.push_back(src1[j]);
                opt_src2.push_back(src2[j]);
            }
            src1.assign(wordsPerMatrixRow, 0);
            src2.assign(wordsPerMatrixRow, 0);
            bool isLastIteration = (std::next(it_u) == adjList.end()) && (std::next(it_v) == neighborsU.end());
            if((words + wordsPerMatrixRow > operandMaxNumberOfWords) || (isLastIteration && words > 0)){
                cout << "-------------itr[" << iterations << "]-------------" << endl;
                cout << "Number of words that are processed in this iteration: " << words << endl;
                std::vector<unsigned int> dst(operandMaxNumberOfWords);
                std::vector<unsigned int> popCountSrc(operandMaxNumberOfWords);
                int sum = vectorAndPopCntRedSum((uint64_t) words, opt_src1, opt_src2, dst, popCountSrc);

                if(sum < 0)
                    return -1;
                
                words = 0;
                iterations++;
                // Reset the vectors
                opt_src1.clear();
                opt_src2.clear();
                dst.clear();
                popCountSrc.clear();
                count += sum;
            }
        }
        if (i % step == 0) {
            std::cout << "Progress: " << (i * 100 / V) << "\% rows completed." << std::endl;
        }
        i++;
    }
    // Each triangle is counted 6 times (once at each vertex), so divide the count by 6
    cout << "Host time (if): " << std::fixed << std::setprecision(3) << host_time_if << " ms." << endl;
    return count / 6;
}

int cpuTrianglesAdjMatrix(const vector<vector<bool>>& adjMatrix) {
    int count = 0;
    int V = adjMatrix.size();
    for (int i = 0; i < V; ++i) {
        for (int j = 0; j < V; ++j) {
            for (int k = 0; k < V; ++k) {
                if (adjMatrix[i][j] && adjMatrix[j][k] && adjMatrix[k][i]) {
                    ++count;
                }
            }
        }
    }
    return count / 6; // Each triangle is counted 6 times (once at each vertex), so divide the count by 6
}

// Function to convert edge list to adjacency list
unordered_map<int, unordered_set<int>> convertToAdjList(const vector<pair<int, int>>& edgeList) {
    unordered_map<int, unordered_set<int>> adjList;

    for (const auto& edge : edgeList) {
        int u = edge.first;
        int v = edge.second;
        adjList[u].insert(v);
        adjList[v].insert(u); // Assuming the graph is undirected
    }

    return adjList;
}

// Function to count triangles using adjacency list
int countTrianglesAdjList(const unordered_map<int, unordered_set<int>>& adjList) {
    int triangleCount = 0;

    for (const auto& [u, neighborsU] : adjList) {
        for (const int& v : neighborsU) {
            if (u < v) { // Ensure each triangle is counted once
                for (const int& w : adjList.at(v)) {
                    if (v < w && neighborsU.find(w) != neighborsU.end()) {
                        triangleCount++;
                    }
                }
            }
        }
    }

    return triangleCount;
}

int main(int argc, char** argv) {
    try {
        struct Params params = getInputParams(argc, argv);
        cout << "Running triangle count on input graph file: " << params.inputFile << endl;
        // Read edge list from JSON file
        vector<pair<int, int>> edgeList = readEdgeList(params.inputFile);
        
        // Determine the number of nodes
        unordered_set<int> nodes;
        for (const auto& edge : edgeList) {
            nodes.insert(edge.first);
            nodes.insert(edge.second);
        }
        int numNodes = nodes.size();
        cout << "Number of nodes: " << numNodes << endl;

        cout << "-----------createDevice-----------" << endl;
        // create device
        if (!createDevice(params.configFile))
            return 1;

        // get device parameters 
        PimDeviceProperties deviceProp;
        PimStatus status = pimGetDeviceProperties(&deviceProp);
        assert(status == PIM_OK);
        uint64_t words_per_device = (uint64_t) deviceProp.numRanks * deviceProp.numBankPerRank * deviceProp.numSubarrayPerBank * deviceProp.numRowPerSubarray * deviceProp.numColPerSubarray / BITS_PER_INT;
        
        int pimTriCount = 0;
        if (INPUT_MODE == ADJ_MATRIX) {
            // Convert edge list to adjacency matrix
            cout << "-----------convertToAdjMatrix-----------" << endl;
            vector<vector<bool>> adjMatrix = edgeListToAdjMatrix(edgeList, numNodes);
            // Convert adjacency matrix to bitwise adjacency matrix
            cout << "-----------convertToBitwiseAdjMatrix-----------" << endl;
            vector<vector<UINT32>> bitAdjMatrix = convertToBitwiseAdjMatrix(adjMatrix);
            cout << "-----------Start running on PIM-----------" << endl;
            if (params.algorithm == "BASELINE"){
                pimTriCount = run_adjmatrix(adjMatrix, bitAdjMatrix, words_per_device);
            }
            else if (params.algorithm == "SERIAL"){
                pimTriCount = run_adjmatrix_serial(adjMatrix, bitAdjMatrix, words_per_device);
            }
            else if (params.algorithm == "TILED"){
                pimTriCount = run_adjmatrix_tiled(adjMatrix, bitAdjMatrix, words_per_device, params.tileSize);
            }
            else if (params.algorithm == "SLICED"){
                pimTriCount = run_adjmatrix_sliced(adjMatrix, bitAdjMatrix, words_per_device, params.sliceNum);
            }
            else {
                fprintf(stderr, "\nINVALID ALGORITHM!!");
            }
            if (params.shouldVerify){
                //run on cpu
                cout << "-----------Triangle Count Verification-----------" << endl;
                int cpuTriCount = cpuTrianglesAdjMatrix(adjMatrix);
                if (cpuTriCount != pimTriCount)
                    printf("PIM count (%d) does not match CPU count(%d)\n", pimTriCount, cpuTriCount);
                else
                    printf("PIM count (%d) matches CPU count(%d)\n", pimTriCount, cpuTriCount);
            }
        } else if (INPUT_MODE == ADJ_LIST) {
            // Convert edge list to adjacency list
            cout << "-----------convertToAdjList-----------" << endl;
            unordered_map<int, unordered_set<int>> adjList = convertToAdjList(edgeList);
            cout << "-----------Start running on PIM-----------" << endl;
            pimTriCount = run_adjlist(adjList, words_per_device);
            if (params.shouldVerify){
                //run on cpu
                cout << "-----------Triangle Count Verification-----------" << endl;
                int cpuTriCount = countTrianglesAdjList(adjList);
                if (cpuTriCount != pimTriCount)
                    printf("PIM count (%d) does not match CPU count(%d)\n", pimTriCount, cpuTriCount);
                else
                    printf("PIM count (%d) matches CPU count(%d)\n", pimTriCount, cpuTriCount);
            }
        }

        //stats
        pimShowStats();

    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
    }
    return 0;
}
