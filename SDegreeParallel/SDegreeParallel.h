#ifndef SDEGREEPARALLEL_H
#define SDEGREEPARALLEL_H
#include <vector>
#include <algorithm>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <fstream>
#include "set_intersection.h"
#include <sys/time.h>
using namespace std;
class SDegreeParallel {
public:
    int n, m;
    vector<int> csr_deg;
    vector<int> csr_edge;
    int max_deg;
public:
    SDegreeParallel(string graph_path);
    long long kclique(int k,int threads);
    long long kclique_main(int k, int* candidate,int node_cnt, int* vertex_list,int* edge_list);
    void SDegreeList(const int* vertex_list, const int* edge_list, const int* candidate, int node_cnt, int level, long long& res, int** sub_candidate);
};
#endif //SDEGREEPARALLEL_H
