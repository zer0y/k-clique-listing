#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <sys/time.h>

#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed
//unsigned** mask;
unsigned block_num;
const unsigned block_size = 24;
const unsigned** construct_mask() {
	unsigned** mask_t;
	unsigned mask_size = (1 << block_size);
	mask_t = (unsigned**)malloc(mask_size * sizeof(unsigned*));
	unsigned i;
	for (i = 0; i < mask_size; i++) {
		mask_t[i] = (unsigned*)malloc((block_size + 1) * sizeof(unsigned));
	}
	for (i = 0; i < mask_size; i++) {
		unsigned t = i, bit = 0, cnt = 1;
		while (t != 0) {
			if (t & 1) {
				mask_t[i][cnt] = bit;
				cnt++;
			}
			bit++;
			t >>= 1;
		}
		mask_t[i][0] = cnt - 1;
	}
	return (const unsigned**)mask_t;
}
static const unsigned** mask = construct_mask();
typedef struct {
	unsigned s;
	unsigned t;
} edge;

typedef struct {
	unsigned node;
	unsigned deg;
} nodedeg;


typedef struct {

	unsigned n;//number of nodes
	unsigned e;//number of edges
	edge* edges;//list of edges

	unsigned* ns;//ns[l]: number of nodes in G_l
	unsigned** d;//d[l]: degrees of G_l
	unsigned* cd, * cdsub;//cumulative degree: (starts with 0) length=n+1
	unsigned* adj, * adjsub;//truncated list of neighbors
	unsigned char* lab;//lab[i] label of node i
	unsigned** sub;//sub[l]: nodes in G_l
  unsigned* valid;
  unsigned* rank;
} specialsparse;

typedef struct {
	unsigned id;
	unsigned degree;
} iddegree;

iddegree* ig;
specialsparse* subg;
int * color, * ind, * loc, * C;
unsigned* cd0, * adj0, * dsub, * Index, max, K, flag;


int cmp(const void* a, const void* b)
{
	iddegree* x = (iddegree*)a, * y = (iddegree*)b;
	return y->degree - x->degree;
}

int cmpedge(const void* a, const void* b) {
	edge* x = (edge*)a, * y = (edge*)b;
	if (y->s == x->s) {
		return x->t - y->t;
	}
	return x->s - y->s;
}

int cmpadj(const void* a, const void* b)
{
	int* x = (int*)a, * y = (int*)b;
	return color[Index[*y]] - color[Index[*x]];
}

void freespecialsparse(specialsparse* g, unsigned char k) {
	unsigned char i;
	free(g->ns);
	for (i = 2; i < k + 1; i++) {
		free(g->d[i]);
		free(g->sub[i]);
	}
	free(g->d);
	free(g->sub);
	free(g->edges);
	free(g->lab);
	free(g->cd);
	free(g->adj);
	free(g);
}

//Compute the maximum of three unsigned integers.
unsigned int max3(unsigned int a, unsigned int b, unsigned int c);
inline unsigned int max3(unsigned int a, unsigned int b, unsigned int c) {
	a = (a > b) ? a : b;
	return (a > c) ? a : c;
}

specialsparse* readedgelist(char* edgelist) {
	unsigned e1 = NLINKS;
	specialsparse* g = (specialsparse*)malloc(sizeof(specialsparse));
	FILE* file;

	g->n = 0;
	g->e = 0;
	file = fopen(edgelist, "r");
	g->edges = (edge*)malloc(e1 * sizeof(edge));
	while (fscanf(file, "%u %u", &(g->edges[g->e].s), &(g->edges[g->e].t)) == 2) {//Add one edge
		g->n = max3(g->n, g->edges[g->e].s, g->edges[g->e].t);
		g->e++;
		if (g->e == e1) {
			e1 += NLINKS;
			g->edges = (edge*)realloc(g->edges, e1 * sizeof(edge));
		}
	}
	fclose(file);
	g->n++;

	g->edges = (edge*)realloc(g->edges, g->e * sizeof(edge));

	return g;
}

void relabel(specialsparse *g) {
	unsigned i, source, target, tmp;
  unsigned e0 = 0;
	for (i = 0; i < g->e; i++) {
    unsigned s = g->edges[i].s, t = g->edges[i].t;
    if(g->valid[s] == 0 || g->valid[t] == 0) continue;
		source = g->rank[s];
		target = g->rank[t];
    
		if (source < target) {
			tmp = source;
			source = target;
			target = tmp;
		}
		g->edges[e0].s = source;
		g->edges[e0].t = target;
    e0++;
	}
  g->e = e0;

}

///// CORE ordering /////////////////////

typedef struct {
	unsigned key;
	unsigned value;
} keyvalue;

typedef struct {
	unsigned n_max;	// max number of nodes.
	unsigned n;	// number of nodes.
	unsigned *pt;	// pointers to nodes.
	keyvalue *kv; // nodes.
} bheap;


bheap *construct(unsigned n_max) {
	unsigned i;
	bheap *heap = (bheap*)malloc(sizeof(bheap));

	heap->n_max = n_max;
	heap->n = 0;
	heap->pt = (unsigned*)malloc(n_max * sizeof(unsigned));
	for (i = 0; i < n_max; i++) heap->pt[i] = -1;
	heap->kv = (keyvalue*)malloc(n_max * sizeof(keyvalue));
	return heap;
}

void swap(bheap *heap, unsigned i, unsigned j) {
	keyvalue kv_tmp = heap->kv[i];
	unsigned pt_tmp = heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key] = heap->pt[heap->kv[j].key];
	heap->kv[i] = heap->kv[j];
	heap->pt[heap->kv[j].key] = pt_tmp;
	heap->kv[j] = kv_tmp;
}

void bubble_up(bheap *heap, unsigned i) {
	unsigned j = (i - 1) / 2;
	while (i > 0) {
		if (heap->kv[j].value > heap->kv[i].value) {
			swap(heap, i, j);
			i = j;
			j = (i - 1) / 2;
		}
		else break;
	}
}

void bubble_down(bheap *heap) {
	unsigned i = 0, j1 = 1, j2 = 2, j;
	while (j1 < heap->n) {
		j = ((j2 < heap->n) && (heap->kv[j2].value < heap->kv[j1].value)) ? j2 : j1;
		if (heap->kv[j].value < heap->kv[i].value) {
			swap(heap, i, j);
			i = j;
			j1 = 2 * i + 1;
			j2 = j1 + 1;
			continue;
		}
		break;
	}
}

void insert(bheap *heap, keyvalue kv) {
	heap->pt[kv.key] = (heap->n)++;
	heap->kv[heap->n - 1] = kv;
	bubble_up(heap, heap->n - 1);
}

void update(bheap *heap, unsigned key) {
	unsigned i = heap->pt[key];
	if (i != (unsigned)-1) {
		((heap->kv[i]).value)--;
		bubble_up(heap, i);
	}
}

keyvalue popmin(bheap *heap) {
	keyvalue min = heap->kv[0];
	heap->pt[min.key] = -1;
	heap->kv[0] = heap->kv[--(heap->n)];
	heap->pt[heap->kv[0].key] = 0;
	bubble_down(heap);
	return min;
}

//Building the heap structure with (key,value)=(node,degree) for each node
bheap* mkheap(unsigned n, unsigned *v) {
	unsigned i;
	keyvalue kv;
	bheap* heap = construct(n);
	for (i = 0; i < n; i++) {
		kv.key = i;
		kv.value = v[i];
		insert(heap, kv);
	}
	return heap;
}

void freeheap(bheap *heap) {
	free(heap->pt);
	free(heap->kv);
	free(heap);
}

//computing degeneracy ordering and core value
void ord_core(specialsparse* g) {
	unsigned i, j, r = 0, n = g->n;
	keyvalue kv;
	bheap *heap;

	unsigned *d0 = (unsigned*)calloc(g->n, sizeof(unsigned));
	unsigned *cd0 = (unsigned*)malloc((g->n + 1) * sizeof(unsigned));
	unsigned *adj0 = (unsigned*)malloc(2 * g->e * sizeof(unsigned));
	for (i = 0; i < g->e; i++) {
		d0[g->edges[i].s]++;
		d0[g->edges[i].t]++;
	}
	cd0[0] = 0;
	for (i = 1; i < g->n + 1; i++) {
		cd0[i] = cd0[i - 1] + d0[i - 1];
		d0[i - 1] = 0;
	}
	for (i = 0; i < g->e; i++) {
		adj0[cd0[g->edges[i].s] + d0[g->edges[i].s]++] = g->edges[i].t;
		adj0[cd0[g->edges[i].t] + d0[g->edges[i].t]++] = g->edges[i].s;
	}

	heap = mkheap(n, d0);

	g->rank = (unsigned*)malloc(g->n * sizeof(unsigned));
  g->valid = (unsigned*)malloc(g->n * sizeof(unsigned));
	for (i = 0; i < g->n; i++) {
		kv = popmin(heap);
    if(kv.value >= K-1) flag = 1;
    g->valid[kv.key] = flag;
		g->rank[kv.key] = n - (++r);
		for (j = cd0[kv.key]; j < cd0[kv.key + 1]; j++) {
			update(heap, adj0[j]);
		}
	}
	freeheap(heap);
	free(d0);
	free(cd0);
	free(adj0);
}

//////////////////////////
//Building the special graph structure
void mkspecial(specialsparse* g, unsigned char k) {
	unsigned i, ns;
	unsigned* d, * sub;
	unsigned char* lab;

	d = (unsigned*)calloc(g->n, sizeof(unsigned));
	qsort(g->edges, g->e, sizeof(g->edges[0]), cmpedge);
	for (i = 0; i < g->e; i++) {
		d[g->edges[i].s]++;
	}

	g->cd = (unsigned*)malloc((g->n + 1) * sizeof(unsigned));
	ns = 0;
	g->cd[0] = 0;
	max = 0;
	sub = (unsigned*)malloc(g->n * sizeof(unsigned));
	lab = (unsigned char*)malloc(g->n * sizeof(unsigned char));
	for (i = 1; i < g->n + 1; i++) {
		g->cd[i] = g->cd[i - 1] + d[i - 1];
		max = (max > d[i - 1]) ? max : d[i - 1];
		sub[ns++] = i - 1;
		d[i - 1] = 0;
		lab[i - 1] = k;
	}
	printf("max degree = %u\n", max);

	subg = (specialsparse*)malloc(sizeof(specialsparse));
	subg->cd = (unsigned*)malloc((max + 1) * sizeof(unsigned));
	subg->adj = (unsigned*)malloc((max * (max - 1) / 2 + 1) * sizeof(unsigned));
	subg->ns = (unsigned*)malloc((k + 1) * sizeof(unsigned));

	subg->d = (unsigned**)malloc((k + 1) * sizeof(unsigned*));
	subg->sub = (unsigned**)malloc((k + 1) * sizeof(unsigned*));
	
	for (int i = 2; i <= k; i++) {
		subg->sub[i] = (unsigned*)malloc(max * sizeof(unsigned));
	}

	C = (int*)malloc(max * sizeof(int));
	ig = (iddegree*)malloc(max * sizeof(iddegree));
	adj0 = (unsigned*)malloc(2 * (max * (max - 1) / 2) * sizeof(unsigned));
	ind = (int*)malloc(g->n * sizeof(int));
	memset(ind, -1, g->n * sizeof(int));
	loc = (int*)malloc(max * sizeof(int));
	dsub = (unsigned*)calloc(max, sizeof(unsigned));
	Index = (unsigned*)malloc(max * sizeof(unsigned));
	color = (int*)malloc(max * sizeof(int));
	cd0 = (unsigned*)malloc((max + 1) * sizeof(unsigned));

	g->adj = (unsigned*)malloc(g->e * sizeof(unsigned));

	for (i = 0; i < g->e; i++) {
		g->adj[g->cd[g->edges[i].s] + d[g->edges[i].s]++] = g->edges[i].t;
	}
	g->ns = (unsigned*)malloc((k + 1) * sizeof(unsigned));
	g->ns[k] = ns;
  
	g->d = (unsigned**)malloc((k + 1) * sizeof(unsigned*));
	g->sub = (unsigned**)malloc((k + 1) * sizeof(unsigned*));
	for (i = 2; i < k; i++) {
		g->d[i] = (unsigned*)malloc(g->n * sizeof(unsigned));
		g->sub[i] = (unsigned*)malloc(max * sizeof(unsigned));
	}
	g->d[k] = d;
	g->sub[k] = sub;

	g->lab = lab;
}

unsigned block_intersect(unsigned* start_a, unsigned* start_b, unsigned* c) {
	unsigned i;
	unsigned size_c = 0;
	for (i = 0; i < block_num; i++) {
		c[i] = start_a[i] & start_b[i];
		size_c += mask[c[i]][0];
	}
	return size_c;
}

unsigned block_intersect_count(unsigned* start_a, unsigned* start_b) {
	unsigned i, c;
	unsigned size_c = 0;
	for (i = 0; i < block_num; i++) {
		c = start_a[i] & start_b[i];
		size_c += mask[c][0];
	}
	return size_c;
}

void mkspecial_sub(specialsparse* g, unsigned char k) {
	unsigned i, j;
	
	block_num = max / block_size + (max % block_size != 0);

	g->cd[0] = 0;
	
	for (i = 1; i < g->n; i++) {
		g->cd[i] = block_num * i;
	}
	memset(g->adj, 0, sizeof(unsigned) * (block_num * g->n));
	//memset(g->sub[k],0,sizeof(unsigned)*(block_num));
	/*for(i=0;i<g->n;i++){
	  id_num=i/block_size, id_bit=i%block_size;
	  g->sub[k][id_num]|=(1<<id_bit);
	}*/
	j = max / block_size;
	for (i = 0; i < j; i++) {
		g->sub[k][i] = (1 << block_size) - 1;
	}
	j = max % block_size;
	if (j > 0)
		g->sub[k][block_num - 1] = (1 << j) - 1;
  g->ns[k] = g->n;
}


void kclique(unsigned l, specialsparse* g, unsigned long long* n) {
	unsigned i, j, k, end, u, v, w, cnt, id_num, id_bit, s = -1, t = -1;

	if (l == 2) {
		for (i = 0; i < block_num; i++) {
			u = g->sub[2][i];
			cnt = mask[u][0];
			for (j = 1; j < cnt + 1; j++) {
				v = mask[u][j] + block_size * i;
				(*n) += block_intersect_count(g->adj + g->cd[v], g->sub[2]);
			}
		}
		return;

	}

	if (l == K)
	{

		for (i = 0; i < g->ns[l]; i++) {
			u = g->sub[l][i];
			if (g->d[l][u] < K - 1) continue;
			max = g->d[l][u];
			g->ns[l - 1] = 0;
			end = g->cd[u] + g->d[l][u];
			for (j = g->cd[u]; j < end; j++) {//relabeling nodes and forming U'.
				v = g->adj[j];
				g->lab[v] = l - 1;
				g->sub[l - 1][g->ns[l - 1]++] = v;
				g->d[l - 1][v] = 0;//new degrees
			}

			int cnt = -1, edge_num = 0;
			for (j = 0; j < g->ns[l - 1]; j++)
			{//reodering adjacency list and computing new degrees

				v = g->sub[l - 1][j];
				if (ind[v] == -1)
				{
					ind[v] = ++cnt;
					loc[cnt] = v;
					dsub[cnt] = 0;
				}
				end = g->cd[v] + g->d[l][v];
				for (k = g->cd[v]; k < end; k++)
				{
					w = g->adj[k];
					if (g->lab[w] == l - 1)
					{
						if (ind[w] == -1)
						{
							ind[w] = ++cnt;
							loc[cnt] = w;
							dsub[cnt] = 0;
						}
						edge_num++;
						dsub[ind[v]]++;
						dsub[ind[w]]++;
					}
				}


			}

			cd0[0] = 0;
			for (j = 1; j < g->ns[l - 1] + 1; j++) {
				cd0[j] = cd0[j - 1] + dsub[j - 1];
				ig[j - 1].id = j - 1;
				ig[j - 1].degree = dsub[j - 1];
				dsub[j - 1] = 0;
			}

			for (j = 0; j < g->ns[l - 1]; j++)
			{
				color[j] = -1;
				v = g->sub[l - 1][j];
				end = g->cd[v] + g->d[l][v];
				for (k = g->cd[v]; k < end; k++)
				{
					w = g->adj[k];
					if (g->lab[w] == l - 1)
					{
						adj0[cd0[ind[v]] + dsub[ind[v]]++] = ind[w];
						adj0[cd0[ind[w]] + dsub[ind[w]]++] = ind[v];
					}
				}
			}

			qsort(ig, g->ns[l - 1], sizeof(ig[0]), cmp);

			for (j = 0; j < g->ns[l - 1]; j++)
			{
				Index[ig[j].id] = j;
				C[j] = 0;
			}

			color[0] = 0;
			unsigned colorNum = 0;

			for (unsigned i = 1; i < g->ns[l - 1]; i++)
			{
				unsigned tmpdegree = ig[i].degree, tmpid = ig[i].id;

				for (unsigned j = 0; j < tmpdegree; j++)
				{
					int now = Index[adj0[cd0[tmpid] + j]];
					if (color[now] != -1)
						C[color[now]] = 1;
				}
				for (unsigned j = 0; j < ig[0].degree + 1; j++)
					if (C[j] == 0)
					{
						color[i] = j;
						colorNum = j > colorNum ? j : colorNum;
						break;
					}

				for (unsigned j = 0; j < tmpdegree; j++)
				{
					int now = Index[adj0[cd0[tmpid] + j]];
					if (color[now] != -1)
						C[color[now]] = 0;
				}

			}
			subg->n = g->ns[l - 1];

			mkspecial_sub(subg, l - 1);
			for (j = 0; j < g->ns[l - 1]; j++)
			{

				v = g->sub[l - 1][j];
				end = g->cd[v] + g->d[l][v];
				for (k = g->cd[v]; k < end; k++)
				{
					w = g->adj[k];
					if (g->lab[w] == l - 1)
					{
						if (color[Index[ind[v]]] < color[Index[ind[w]]])
						{
							s = ind[w], t = ind[v];
						}
						else if (color[Index[ind[v]]] == color[Index[ind[w]]])
						{
							if (ig[Index[ind[v]]].id < ig[Index[ind[w]]].id)
							{
								s = ind[v], t = ind[w];
							}
							else
							{
								s = ind[w], t = ind[v];
							}
						}
						else if (color[Index[ind[v]]] > color[Index[ind[w]]])
						{
							s = ind[v], t = ind[w];
						}
						id_num = t / block_size, id_bit = t % block_size;
						subg->adj[subg->cd[s] + id_num] |= (1 << id_bit);
					}
				}
			}
			subg->e = edge_num;
			kclique(l - 1, subg, n);

			for (j = 0; j < g->ns[l - 1]; j++) {
				ind[loc[j]] = -1;
				v = g->sub[l - 1][j];
				g->lab[v] = l;
			}

		}
	}

	else
	{
		if (l > g->ns[l])
			return;
		for (i = 0; i < block_num; i++) {
			u = g->sub[l][i];
			cnt = mask[u][0];
			for (j = 1; j < cnt + 1; j++) {
				v = mask[u][j] + block_size * i;;
				if (color[Index[v]] < (int)(l - 1))
					continue;
				g->ns[l - 1] = block_intersect(g->adj + g->cd[v], g->sub[l], g->sub[l - 1]);
				kclique(l - 1, g, n);
			}
		}
	}
}

int main(int argc, char** argv) {

	specialsparse* g;
	unsigned char k = atoi(argv[1]);
	K = k;
	unsigned long long n;
  flag = 0;
	struct timeval time_start;
	struct timeval time_end;
	printf("Reading edgelist from file %s\n", argv[2]);
	gettimeofday(&time_start, NULL);
	g = readedgelist(argv[2]);
	printf("Number of nodes = %u\n", g->n);
	printf("Number of edges = %u\n", g->e);

	gettimeofday(&time_end, NULL);
	double read_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;
	printf("- Read Time = %fs\n", read_time / 1000);
  ord_core(g);
  relabel(g);
  mkspecial(g, k);
	printf("Iterate over all cliques\n");
	gettimeofday(&time_start, NULL);
	n = 0;
	kclique(k, g, &n);
	gettimeofday(&time_end, NULL);
	double list_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;
	printf("Number of %u-cliques: %llu\n", k, n);
	printf("- List Time = %fs\n", list_time / 1000);
	freespecialsparse(g, k);
	free(dsub);
	return 0;
}