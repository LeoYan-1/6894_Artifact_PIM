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
#if defined(_OPENMP)
#include <omp.h>
#endif

#include "../util.h"
#include "libpimsim.h"

#define BITS_PER_INT 32

#define NUM_SUBARRAY 512
#define ROW_SIZE 8192
#define COL_SIZE 8192
#define WORDS_PER_ROW (ROW_SIZE * NUM_SUBARRAY)

typedef uint32_t UINT32;


using namespace std;


// Params ---------------------------------------------------------------------
typedef struct Params
{
  uint64_t vectorLength;
  char *configFile;
  char *inputFile;
  bool shouldVerify;
} Params;

void usage()
{
  fprintf(stderr,
          "\nUsage:  ./add [options]"
          "\n"
          "\n    -l    input size (default=8M elements)"
          "\n    -c    dramsim config file"
          "\n    -i    input file containing two vectors (default=generates vector with random numbers)"
          "\n    -v    t = verifies PIM output with host output. (default=false)"
          "\n");
}

template <typename T>
void printNestedVector(const std::vector<std::vector<T>>& nestedVec) {
    for (const std::vector<T>& innerVec : nestedVec) {
        for (const T& element : innerVec) {
            std::cout << element << " ";
        }
        std::cout << std::endl; // Print a newline after each inner vector
    }
}

struct Params getInputParams(int argc, char **argv)
{
  struct Params p;
  p.vectorLength = 65536;
  p.configFile = nullptr;
  p.inputFile = nullptr;
  p.shouldVerify = false;

  int opt;
  while ((opt = getopt(argc, argv, "h:l:c:i:v:")) >= 0)
  {
    switch (opt)
    {
    case 'h':
      usage();
      exit(0);
      break;
    case 'l':
      p.vectorLength = strtoull(optarg, NULL, 0);
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
    default:
      fprintf(stderr, "\nUnrecognized option!\n");
      usage();
      exit(0);
    }
  }
  return p;
}

// Function to convert edge list to adjacency matrix
vector<vector<int>> edgeListToAdjMatrix(const vector<pair<int, int>>& edgeList, int numNodes) {
    vector<vector<int>> adjMatrix(numNodes+1, vector<int>(numNodes+1, 0));

    for (const auto& edge : edgeList) {
        int u = edge.first;
        int v = edge.second;
        adjMatrix[u][v] = adjMatrix[v][u] = 1; // assuming undirected graph
    }

    return adjMatrix;
}

// Function to convert standard adjacency matrix to bitwise adjacency matrix
vector<vector<UINT32>> convertToBitwiseAdjMatrix(const vector<vector<int>>& adjMatrix) {
    int V = adjMatrix.size();
    int numInts = (V + BITS_PER_INT - 1) / BITS_PER_INT; // Number of 32-bit integers needed per row

    vector<vector<UINT32>> bitAdjMatrix(V, vector<UINT32>(numInts, 0));

    for (int i = 0; i < V; ++i) {
        for (int j = 0; j < V; ++j) {
            if (adjMatrix[i][j]) {
                bitAdjMatrix[i][j / BITS_PER_INT] |= (1 << (j % BITS_PER_INT));
            }
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

int vectorAndPopCntRedSum(uint64_t numElements, std::vector<unsigned int> &src1, std::vector<unsigned int> &src2, std::vector<unsigned int> &dst){
    unsigned bitsPerElement = sizeof(int) * 8;

    PimObjId srcObj1 = pimAlloc(PIM_ALLOC_V1, numElements, bitsPerElement, PIM_INT32);
    if (srcObj1 == -1)
    {
        std::cout << "Abort" << std::endl;
        return -1;
    }

    PimObjId srcObj2 = pimAllocAssociated(PIM_ALLOC_V1, numElements, bitsPerElement, srcObj1, PIM_INT32);
    if (srcObj2 == -1)
    {
        std::cout << "Abort" << std::endl;
        return -1;
    }

    PimObjId dstObj = pimAllocAssociated(PIM_ALLOC_V1, numElements, bitsPerElement, srcObj1, PIM_INT32);
    if (dstObj == -1)
    {
        std::cout << "Abort" << std::endl;
        return -1;
    }

    PimObjId popCntResObj = pimAllocAssociated(PIM_ALLOC_V1, numElements, bitsPerElement, srcObj1, PIM_INT32);
    if (popCntResObj == -1)
    {
        std::cout << "Abort" << std::endl;
        return -1;
    }

    PimStatus status = pimCopyHostToDevice(PIM_COPY_V, (void *)src1.data(), srcObj1);
    if (status != PIM_OK)
    {
        std::cout << "Abort" << std::endl;
        return -1;
    }

    status = pimCopyHostToDevice(PIM_COPY_V, (void *)src2.data(), srcObj2);
    if (status != PIM_OK)
    {
        std::cout << "Abort" << std::endl;
        return -1;
    }
    
    status = pimAnd(srcObj1, srcObj2, dstObj);
    if (status != PIM_OK)
    {
        std::cout << "Abort" << std::endl;
        return -1;
    }

    status = pimPopCount(dstObj, popCntResObj);
    if (status != PIM_OK)
    {
        std::cout << "Abort" << std::endl;
        return -1;
    }

    int sum = 0;
    status = pimRedSum(popCntResObj, &sum);
    if (status != PIM_OK)
    {
        std::cout << "Abort" << std::endl;
        return -1;
    }

    pimFree(srcObj1);
    pimFree(srcObj2);
    pimFree(dstObj);

    return sum;
}

int run(const vector<vector<int>>& adjMatrix, const vector<vector<UINT32>>& bitAdjMatrix) {
    int count = 0;
    int V = bitAdjMatrix.size();
    // unsigned numElements = V;
    int numElements = (V + BITS_PER_INT - 1) / BITS_PER_INT; // Number of 32-bit integers needed per row
    cout << "number of ndoes: " << V << endl;
    cout << "numElem: " << numElements << endl;
    assert(numElements <= WORDS_PER_ROW && "Number of vertices cannot exceed WORDS_PER_ROW");

    for (int i = 0; i < V; ++i) {
        for (int j = 0; j < V; ++j) {
            if (adjMatrix[i][j]) { // If there's an edge between i and j
                // int l = j / BITS_PER_INT;
                std::vector<unsigned int> src1(numElements);
                std::vector<unsigned int> src2(numElements);
                std::vector<unsigned int> dest(numElements);
                // cout << "i: " << i << ", j: " << j << endl;
                for (int k = 0; k < numElements; ++k) {
                    // dotProduct += __builtin_popcount(bitAdjMatrix[i][k] & bitAdjMatrix[j][k]);
                    src1[k] = bitAdjMatrix[i][k];
                    src2[k] = bitAdjMatrix[j][k];
                } 
                // cout << "src1: " << bitset<32>(src1[0]) << endl;
                // cout << "src2: " << bitset<32>(src2[0]) << endl;
                int sum = vectorAndPopCntRedSum((uint64_t) numElements, src1, src2, dest);
                // cout << "sum: " << sum << endl;
                //redsum
                count += sum;
            }
        }
    }

    cout << "bit32TriangleCount: " << count / 6 << endl;
    // Each triangle is counted three times (once at each vertex), so divide the count by 3
    return count / 6;
}

int main(int argc, char** argv) {
    try {
        struct Params params = getInputParams(argc, argv);
        // Read edge list from JSON file
        string filename = argv[1];
        vector<pair<int, int>> edgeList = readEdgeList(filename);
        
        // Determine the number of nodes
        unordered_set<int> nodes;
        for (const auto& edge : edgeList) {
            nodes.insert(edge.first);
            nodes.insert(edge.second);
        }
        int numNodes = nodes.size();
        cout << "Number of nodes: " << numNodes << endl;

        // Convert edge list to adjacency matrix
        vector<vector<int>> adjMatrix = edgeListToAdjMatrix(edgeList, numNodes);
        cout << "Adjacency Matrix size:" << adjMatrix.size() << endl;

        // cout << "-----edgelist-----\n";
        // printNestedVector(adjMatrix);
        // cout << "-----edgelist-----\n";

        vector<vector<UINT32>> bitAdjMatrix = convertToBitwiseAdjMatrix(adjMatrix);

        //  cout << "-----bitAdjMatrix-----\n";
        // for (const auto& row : bitAdjMatrix) {
        //     for (UINT32 val : row) {
        //         cout << bitset<32>(val) << " ";
        //     }
        //     cout << endl;
        // }
        // cout << "-----bitAdjMatrix-----\n";

        if (!createDevice(params.configFile))
            return 1;
        //run
        run(adjMatrix, bitAdjMatrix);

        //stats
        pimShowStats();

    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
    }
    return 0;
}