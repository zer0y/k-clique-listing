#include <iostream>
#include "SDegreeParallel.h"
using namespace std;
int main(int argc,char** argv) {
    if(argc!=4){
        printf("Please input the command as </.SDegreeParallel>  <#thread>  <k>  <graph_path>\n");
        return -1;
    }
    int threads=atoi(argv[1]);
    int k=atoi(argv[2]);
    string graph_path=argv[3];
    SDegreeParallel g(graph_path);
    g.kclique(k,threads);
    return 0;
}
