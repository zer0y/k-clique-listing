#include "SDegree.h"


SDegree::SDegree(string graph_path) {

    int  u, v;
    char line[512];
    char spliter = ' ';
    vector<vector<int>> data;
    FILE* data_file;
    data_file = fopen(graph_path.c_str(), "r");
    while (!feof(data_file)) {
        if (fgets(line, 512, data_file) == NULL) continue;
        if (!isdigit(line[0])) continue;
        if (strlen(line) < 3) continue;
        sscanf(line, "%d%c%d", &u, &spliter, &v);
        data.push_back({ u,v });
        n=max(n,max(u,v));
    }
    fclose(data_file);
    std::sort(data.begin(), data.end());
    data.erase(std::unique(data.begin(), data.end()), data.end());
    m = data.size();
    n = n + 1;
    csr_deg.resize(n + 1, 0);
    csr_edge.resize(m);
    int ne = 0;
    for (auto& e : data) {
        csr_deg[e[0] + 1]++;
        csr_edge[ne++] = e[1];
    }
    data.clear();
}

void SDegree::SDegreeList(const int* vertex_list, const int* edge_list, const int* candidate, int node_cnt, int level, long long& res) {
    for (int i = 0; i < node_cnt; i++) {
        int u = candidate[i];
        int start_idx = vertex_list[u], nbr_cnt = vertex_list[u + 1] - vertex_list[u];
        if(nbr_cnt<=level-2) continue;
        if (level == 2) {
            int itscCnt = intersect_simd4x_count(edge_list + start_idx, nbr_cnt, candidate, node_cnt);
            res += (long long)itscCnt;
        }
        else {
            
            //if(edge_list[start_idx+nbr_cnt-(level - 1)]<candidate[0] ) continue;
            int itscCnt = intersect_simd4x(edge_list + start_idx, nbr_cnt, candidate, node_cnt, sub_candidate[level-2]);
            if (itscCnt > level - 2) {
                SDegreeList(vertex_list, edge_list, sub_candidate[level-2], itscCnt, level - 1, res);
            }
        }
    }
}

long long SDegree::kclique_main(int k, int* candidate,int node_cnt, int* vertex_list,int* edge_list){
    long long res=0;
    for (int i = 0; i < node_cnt; i++) {
      int u = candidate[i];
      int start_idx = vertex_list[u], end_idx = vertex_list[u + 1];
      SDegreeList(vertex_list, edge_list, edge_list + start_idx, end_idx - start_idx, k - 1, res);
    }

    return res;
}
long long SDegree::kclique(int k) {
    struct timeval time_start;
    struct timeval time_end;
    vector<int> candidate;
    int max_deg = 0;
    for (int i = 1; i <= n; i++) {
        max_deg = max(max_deg, csr_deg[i]);
        if (csr_deg[i] >= k - 1) {
            candidate.push_back(i - 1);
        }
        csr_deg[i] += csr_deg[i - 1];
    }
    sub_candidate=(int**) malloc(k*sizeof(int*));
    for(int i=0;i<k;i++){
        sub_candidate[i]=(int*)malloc((max_deg+1)*sizeof(int));
    }
    printf("max degree:%d\n",max_deg);
    gettimeofday(&time_start, NULL);
    long long res = kclique_main(k,candidate.data(),candidate.size(),csr_deg.data(),csr_edge.data());
    gettimeofday(&time_end, NULL);
    double list_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;
    printf("Number of %d-cliques: %lld\n", k, res);
	  printf("- List Time = %fs\n", list_time / 1000);
    for(int i=0;i<k;i++){
        free(sub_candidate[i]);
    }
    free(sub_candidate);
    return res;
}


