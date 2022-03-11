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
    int max_deg;
    vector<int> csr_deg;
    vector<int> csr_edge;
    
public:
    SDegree(string graph_path);
    long long kclique(int k);
    long long kclique_main(int k, int* vertex_list,int* edge_list,int** sub_candidate);
    void SDegreeList(int** sub_candidate,const int* vertex_list, const int* edge_list, const int* candidate, int node_cnt, int level, long long& res);
};
    
#endif //SDEGREE_H
