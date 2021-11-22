#include "Graph.h"

Graph::Graph(string data_path) {
    int  u, v;
    n = 0;
    char line[512];
    char spliter = ' ';
    vector<pair<int, int>> data;
    FILE* data_file;
    data_file = fopen(data_path.c_str(), "r");
    while (!feof(data_file)) {
        if (fgets(line, 512, data_file) == NULL) continue;
        if (!isdigit(line[0])) continue;
        if (strlen(line) < 3) continue;
        sscanf(line, "%d%c%d", &u, &spliter, &v);
        if (u > v)
            data.push_back({ u,v });
        else
            data.push_back({ v,u });
        n = max(n, max(u, v));
    }
    n = n + 1;
    fclose(data_file);
    sort(data.begin(), data.end());
    data.erase(std::unique(data.begin(), data.end()), data.end());
    edges.resize(n);
    degrees.resize(n, 0);
    for (auto& e : data) {
        degrees[e.first]++;
        degrees[e.second]++;
    }

    for (auto& e : data) {
        edges[e.first].push_back(e.second);
        edges[e.second].push_back(e.first);
    }
    m = data.size();
    data.clear();
    nodes.resize(n,vector<int>(2,-1));
    new_id.resize(n,-1);
    org2newid.resize(n);
    for(int i=0;i<n;i++){
        org2newid[i]=i;
    }
}


long long Graph::combo(long long a, long long b) {
    long long  ans = 1;
    for (long long i = 0; i < b; i++)
        ans *= (a - i);
    for (long long i = 1; i <= b; i++) ans /= i;
    return ans;
}

void Graph::preCore(int k) {
    queue<int> q;
    vector<bool> inq(n, false);
    for (int i = n - 1; i >= 0; i--) {
        if (degrees[i] < k) {
            q.push(i);
            inq[i] = true;
        }
    }
    while (!q.empty()) {
        int id = q.front();
        q.pop();
        for (int u : edges[id]) {
            degrees[u]--;
            if (degrees[u] < k && !inq[u]) {
                inq[u] = true;
                q.push(u);
            }
        }
    }
    vector<int> node_map(n, 0);
    int new_node_cnt = 0;
    for (int i = 0; i < n; i++) {
        if (!inq[i]) {
            node_map[i] = new_node_cnt;
            new_node_cnt++;
        }
    }
    new_node_cnt=0;
    unsigned int new_edge_cnt=0;
    for (int i = 0; i < n; i++) {
        if (!inq[i]) { 
            edges[new_node_cnt].resize(edges[i].size());
            unsigned int  tmp_edge_cnt = 0;
            for (int v : edges[i]) {
                if (!inq[v]) {
                    edges[new_node_cnt][tmp_edge_cnt++] = node_map[v];
                }
            }
            degrees[new_node_cnt] = (int)tmp_edge_cnt;
            edges[new_node_cnt].resize(tmp_edge_cnt);
            new_node_cnt++;
            new_edge_cnt += tmp_edge_cnt;
        }
    }
    edges.resize(new_node_cnt);
    degrees.resize(new_node_cnt);
    n = new_node_cnt;
    m = (int)(new_edge_cnt/2);
}

long long Graph::preList(int k) {
    vector<int> label(n, 0);
    int totcc = 0;
    long long res = 0;
    vector<int> valid(n, 0);
    for (int i = 0; i < n; i++) {
        if (label[i] == 0) {
            ++totcc;
            label[i] = totcc;
            queue<int> q;
            q.push(i);
            long long ccsize = 1, ccdeg = (long long)degrees[i];
            while (!q.empty()) {
                int t = q.front(); q.pop();
                for (int nbr : edges[t]) {
                    if (label[nbr] == 0) {
                        label[nbr] = totcc;
                        ccsize++;
                        ccdeg += (long long)degrees[nbr];
                        q.push(nbr);
                    }
                }
               
            }
            if (ccdeg == (ccsize - 1) * ccsize) {
                if (ccsize >= k) {//contain k-clique 
                    res += combo(ccsize, (long long)k);
                }
                q.push(i);
                valid[i] = -1;
                while (!q.empty()) {
                    int t = q.front(); q.pop();
                    for (int nbr : edges[t]) {
                        if (valid[nbr] == -1)
                            continue;
                        valid[nbr] = -1;
                        q.push(nbr);
                    }
                }
            }
            /*else if(ccdeg < k*(k-1)+2*(k-1)*(ccsize-k)){
                q.push(i);
                valid[i] = -1;
                while (!q.empty()) {
                    int t = q.front(); q.pop();
                    for (int nbr : edges[t]) {
                        if (valid[nbr] == -1)
                            continue;
                        valid[nbr] = -1;
                        q.push(nbr);
                    }
                }
            }*/
            
        }
    }
    vector<int> node_map(n, 0);
    int new_node_cnt = 0;
    for (int i = 0; i < n; i++) {
        //if (valid[label[i]] != -1) {
        if (valid[i] != -1) {
            node_map[i] = new_node_cnt;
            new_node_cnt++;
        }
    }

    new_node_cnt = 0;
    unsigned int new_edge_cnt = 0;
    for (int i = 0; i < n; i++) {
        if (valid[i] != -1) {
            unsigned int  tmp_edge_cnt = 0;
            edges[new_node_cnt].resize(degrees[i]);
            for (int v : edges[i]) {
                if (valid[v] != -1) {
                    edges[new_node_cnt][tmp_edge_cnt++] = node_map[v];
                }
            }
            degrees[new_node_cnt] = tmp_edge_cnt;
            edges[new_node_cnt].resize(tmp_edge_cnt);
            new_node_cnt++;
            new_edge_cnt += tmp_edge_cnt;
        }
    }
    edges.resize(new_node_cnt);
    degrees.resize(new_node_cnt);

    n = new_node_cnt;
    m = (int)(new_edge_cnt/2);
    return res;    
}

void Graph::degreeReorder(){
    
    for(int i=0;i<n;i++){
        nodes[i][0]=i;
        nodes[i][1]=degrees[i];
    }
    stable_sort(nodes.begin(),nodes.end(),[](const std::vector<int> &a, const std::vector<int> &b) {return a[1] > b[1]; });
    for (int i = 0; i < n; ++i){
        new_id[nodes[i][0]]=i;
    }
    for (auto& idx : org2newid) idx = new_id[idx];
    /*

    for(int i=0;i<n;i++){
        deg_map[i]=i;
    }
    sort(deg_map.begin(),deg_map.end(),[&](const int& a, const int& b) -> bool {
        int deg_a = degrees[a];
        int deg_b = degrees[b];
        if (deg_a == deg_b) return a < b;
        return deg_a > deg_b;
    });
    
    for(int i=0;i<n;i++){
        new_id[deg_map[i]]=i;
    }*/
}


void Graph::writeDegreeOrder(string graph) {
    ofstream fout;
    string direced_graph_path = graph.substr(0, graph.find_last_of('.')) + "_sdegree.txt";
    fout.open(direced_graph_path);
    for (int i = 0; i < n; i++) {
        for (int v : edges[i]) {
            if (org2newid[i] > org2newid[v])
                fout << org2newid[i] << " " << org2newid[v] << endl;
        }
    }
    fout.close();
}
