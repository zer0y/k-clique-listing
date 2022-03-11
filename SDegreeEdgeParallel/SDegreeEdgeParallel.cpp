#include "SDegreeEdgeParallel.h"

SDegreeParallel::SDegreeParallel(string graph_path) {

    int  u, v;
    char line[512];
    char spliter = ' ';
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
    for (auto& e: data) {
        csr_deg[e[0] + 1]++;
        csr_edge[ne++] = e[1];
    }
    //data.clear();
}

void SDegreeParallel::SDegreeList(const int* vertex_list, const int* edge_list, const int* candidate, int node_cnt, int level, long long& res, int** sub_candidate) {
    if(level == 1){
        res+=node_cnt;
        return;
    }
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
                SDegreeList(vertex_list, edge_list, sub_candidate[level-2], itscCnt, level - 1, res, sub_candidate);
            }
        }
    }
}

long long SDegreeParallel::kclique_main(int k, int* vertex_list,int* edge_list){
    long long res=0;
    long long res_t=0;
    int** sub_candidate;
    #pragma omp parallel private(res_t,sub_candidate) reduction(+:res)
    {
        sub_candidate=(int**) malloc(k*sizeof(int*));
        for(int i=0;i<k;i++){
            sub_candidate[i]=(int*)calloc((max_deg+1),sizeof(int));
        }
        #pragma omp for schedule(dynamic, 1) nowait
        for (int i=0;i<data.size();i++) {
          int u = data[i][0],v=data[i][1];
          int begin_idx_u = vertex_list[u], end_idx_u = vertex_list[u + 1];
          int begin_idx_v = vertex_list[v], end_idx_v = vertex_list[v + 1];
          int itscCnt = intersect_simd4x(edge_list + begin_idx_v, end_idx_v-begin_idx_v, edge_list + begin_idx_u, end_idx_u-begin_idx_u, sub_candidate[k-3]);
            if (itscCnt > k-3) {
                res_t=0;
                SDegreeList(vertex_list, edge_list, sub_candidate[k-3], itscCnt, k - 2, res_t, sub_candidate);
                res+=res_t;
            }
        }
    }
    return res;
}
long long SDegreeParallel::kclique(int k,int threads) {
    struct timeval time_start;
    struct timeval time_end;
    vector<int> candidate;
    max_deg = 0;
    for (int i = 1; i <= n; i++) {
        max_deg = max(max_deg, csr_deg[i]);
        //if (csr_deg[i] >= k - 1) {
        //candidate.push_back(i - 1);
        //}
        csr_deg[i] += csr_deg[i - 1];
    }
    printf("max degree:%d\n",max_deg);
    omp_set_num_threads(threads);
    gettimeofday(&time_start, NULL);
    long long res = kclique_main(k,csr_deg.data(),csr_edge.data());
    gettimeofday(&time_end, NULL);
    double list_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;
    printf("Number of %d-cliques: %lld\n", k, res);
	  printf("- List Time (%d threads)= %fs\n", threads, list_time / 1000);
    return res;
}


