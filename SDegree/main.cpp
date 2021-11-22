#include <iostream>
#include "SDegree.h"
using namespace std;
int main(int argc,char** argv) {
    if(argc!=3){
        printf("Please input the command as </.SDegree>  <k>  <graph_path>\n");
        return -1;
    }
      
    int k=atoi(argv[1]);
    string graph_path=argv[2];
    SDegree g(graph_path);
    g.kclique(k);
    return 0;
}
