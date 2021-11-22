
#ifndef CPU_GRAPH_H
#define CPU_GRAPH_H
#include <string>
#include <vector>
#include <ctime>
#include <stdio.h>
#include <queue>
#include <fstream>
#include <iostream>
#include <string.h>
#include <algorithm>

using namespace std;
class Graph {
public:
    int n, m, k;
    vector<vector<int>> edges;
    vector<int> degrees;
    vector<vector<int>> nodes;
    vector<int> new_id;
    vector<int> org2newid;
public:
    Graph(string data_path);
    void preCore(int k);
    long long combo(long long n, long long m);
    long long preList(int k);
    void degreeReorder();
    void writeDegreeOrder(string order_path);

};
#endif //CPU_GRAPH_H
#pragma once
