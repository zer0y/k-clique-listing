#ifndef SDEGREE_H
#define SDEGREE_H
#include <vector>
#include <algorithm>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <fstream>
#include "set_intersection.h"
#include <sys/time.h>
using namespace std;
class SDegree {
public:
    int n, m;
    vector<int> csr_deg;
    vector<int> csr_edge;
    int** sub_candidate;
public:
    SDegree(string graph_path);
    long long kclique(int k);
    long long kclique_main(int k, int* candidate,int node_cnt, int* vertex_list,int* edge_list);
    void SDegreeList(const int* vertex_list, const int* edge_list, const int* candidate, int node_cnt, int level, long long& res);
};
#endif //SDEGREE_H
