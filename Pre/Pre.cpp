#include <iostream>
#include "Graph.h"
#include <sys/time.h>
void preSDegree(int k, Graph& g) {
    struct timeval time_start;
    struct timeval time_end;
    double pre_time;
    printf("Before preprocess: #node: %d, #edge:%d\n", g.n,g.m);
    gettimeofday(&time_start, NULL);
    g.preCore(k - 1);
    long long pre_cnt = g.preList(k);
    g.degreeReorder();
    gettimeofday(&time_end, NULL);
    pre_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;
    printf("After preprocess: #node: %d, #edge:%d\n", g.n,g.m);
    printf("Preprocess %d-clique: %lld\n", k, pre_cnt);
    printf("Preprocess time:%.3fms\n", pre_time);
}
int main(int argc, char** argv)
{
    if (argc != 3) {
        printf("Please input the command as <./preprocess>  <k>  <graph_path>\n");
        return -1;
    }
    int k = atoi(argv[1]);
    string graph_path = argv[2];
    Graph g(graph_path);
    //Preprocessing for SDegree
    preSDegree(k, g);
    g.writeDegreeOrder(graph_path);
    return 0;
}

