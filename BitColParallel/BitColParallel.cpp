#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <sys/time.h>

const unsigned NLINKS = 100000000;//maximum number of edges for memory allocation, will increase if needed
//unsigned** mask;
//unsigned* mask_cnt;
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
	unsigned n;//number of nodes
	unsigned e;//number of edges
	edge* edges;//list of edges
	unsigned* rank;//ranking of the nodes according to degeneracy ordering
  unsigned* valid;
} edgelist;

typedef struct {
	unsigned n;
	unsigned* cd;//cumulative degree: (starts with 0) length=n+1
	unsigned* adj;//truncated list of neighbors
	unsigned core;//core value of the graph
} graph;

typedef struct {
	unsigned* n;//n[l]: number of nodes in G_l
	unsigned** d;//d[l]: degrees of G_l
	unsigned* adj;//truncated list of neighbors
	unsigned** nodes;//sub[l]: nodes in G_l
	unsigned core;
	unsigned* color;
	unsigned block_num;
} subgraph;

typedef struct {
	unsigned id;
	unsigned degree;
} iddegree;

typedef struct {
	unsigned id;
	unsigned color;
} idcolor;

unsigned block_intersect(unsigned* start_a, unsigned* start_b, unsigned* c, unsigned block_num) {
	unsigned i;
	unsigned size_c = 0;
	for (i = 0; i < block_num; i++) {
		c[i] = start_a[i] & start_b[i];
	}
	for (i = 0; i < block_num; i++)
		size_c += mask[c[i]][0];
	//size_c += mask_cnt[c[i]];
	return size_c;
}

unsigned block_intersect_count(unsigned* start_a, unsigned* start_b, unsigned block_num) {
	unsigned i, c;
	unsigned size_c = 0;
	for (i = 0; i < block_num; i++) {
		c = start_a[i] & start_b[i];
		size_c += mask[c][0];
		//size_c += mask_cnt[c];
	}
	return size_c;
}

/*
void construct_mask() {
	unsigned mask_size = (1 << block_size);
	mask = (unsigned**)malloc(mask_size * sizeof(unsigned*));
	//mask_cnt = (unsigned*)malloc(mask_size * sizeof(unsigned));
	unsigned i;
	for (i = 0; i < mask_size; i++) {
		mask[i] = (unsigned*)malloc((block_size + 1) * sizeof(unsigned));
	}
	for (i = 0; i < mask_size; i++) {
		unsigned t = i;
		unsigned short bit = 0, cnt = 1;
		while (t != 0) {
			if (t & 1) {
				mask[i][cnt] = bit;
				cnt++;
			}
			bit++;
			t >>= 1;
		}
		mask[i][0]=cnt-1;
		//mask_cnt[i] = cnt - 1;
	}
}
*/
int cmp(const void* a, const void* b)
{
	iddegree* x = (iddegree*)a, * y = (iddegree*)b;
	return y->degree - x->degree;
}

int cmpadj(const void* a, const void* b)
{
	idcolor* x = (idcolor*)a, * y = (idcolor*)b;
	return y->color - x->color;
}

void free_edgelist(edgelist* el) {
	free(el->edges);
	free(el->rank);
	free(el);
}

void free_graph(graph* g) {
	free(g->cd);
	free(g->adj);
	free(g);
}

void free_subgraph(subgraph* sg, unsigned int k) {
	unsigned int i;
	free(sg->n);
	for (i = 2; i < k; i++) {
		free(sg->d[i]);
		free(sg->nodes[i]);
	}
	free(sg->d);
	free(sg->nodes);
	//free(sg->lab);
	free(sg->adj);
	free(sg);
}


//Compute the maximum of three unsigned integers.
unsigned int max3(unsigned int a, unsigned int b, unsigned int c);
inline unsigned int max3(unsigned int a, unsigned int b, unsigned int c) {
	a = (a > b) ? a : b;
	return (a > c) ? a : c;
}

edgelist* readedgelist(char* input) {
	unsigned e1 = NLINKS;
	edgelist* el = (edgelist*)malloc(sizeof(edgelist));
	FILE* file;

	el->n = 0;
	el->e = 0;
	file = fopen(input, "r");
	el->edges = (edge*)malloc(e1 * sizeof(edge));
	while (fscanf(file, "%u %u", &(el->edges[el->e].s), &(el->edges[el->e].t)) == 2) {//Add one edge
		el->n = max3(el->n, el->edges[el->e].s, el->edges[el->e].t);
		el->e++;
		if (el->e == e1) {
			e1 += NLINKS;
			el->edges = (edge*)realloc(el->edges, e1 * sizeof(edge));
		}
	}
	fclose(file);
	el->n++;

	el->edges = (edge*)realloc(el->edges, el->e * sizeof(edge));

	return el;
}

void relabel(edgelist *el) {
	unsigned i, source, target, tmp, e0 = 0;

	for (i = 0; i < el->e; i++) {
    unsigned s = el->edges[i].s, t = el->edges[i].t;
    if(el->valid[s] == 0 || el->valid[t] == 0) continue;
		source = el->rank[s];
		target = el->rank[t];
		if (source < target) {
			tmp = source;
			source = target;
			target = tmp;
		}
		el->edges[e0].s = source;
		el->edges[e0].t = target;
    e0++;
	}
  el->e = e0;

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
void ord_core(edgelist* el, unsigned k) {
	unsigned i, j, r = 0, n = el->n, e = el->e, flag = 0;
	keyvalue kv;
	bheap *heap;

	unsigned *d0 = (unsigned*)calloc(el->n, sizeof(unsigned));
	unsigned *cd0 = (unsigned*)malloc((el->n + 1) * sizeof(unsigned));
	unsigned *adj0 = (unsigned*)malloc(2 * el->e * sizeof(unsigned));
	for (i = 0; i < e; i++) {
		d0[el->edges[i].s]++;
		d0[el->edges[i].t]++;
	}
	cd0[0] = 0;
	for (i = 1; i < n + 1; i++) {
		cd0[i] = cd0[i - 1] + d0[i - 1];
		d0[i - 1] = 0;
	}
	for (i = 0; i < e; i++) {
		adj0[cd0[el->edges[i].s] + d0[el->edges[i].s]++] = el->edges[i].t;
		adj0[cd0[el->edges[i].t] + d0[el->edges[i].t]++] = el->edges[i].s;
	}

	heap = mkheap(n, d0);

	el->rank = (unsigned*)malloc(n * sizeof(unsigned));
  el->valid = (unsigned*)malloc(n * sizeof(unsigned));
	for (i = 0; i < n; i++) {
		kv = popmin(heap);
    if(kv.value >= k-1) flag = 1;
    el->valid[kv.key] = flag;
		el->rank[kv.key] = n - (++r);
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
//Building the special graph
graph* mkgraph(edgelist* el) {
	unsigned i, max;
	unsigned* d;
	graph* g = (graph*)malloc(sizeof(graph));

	d = (unsigned*)calloc(el->n, sizeof(unsigned));

	for (i = 0; i < el->e; i++) {
		d[el->edges[i].s]++;
	}

	g->cd = (unsigned*)malloc((el->n + 1) * sizeof(unsigned));
	g->cd[0] = 0;
	max = 0;
	for (i = 1; i < el->n + 1; i++) {
		g->cd[i] = g->cd[i - 1] + d[i - 1];
		max = (max > d[i - 1]) ? max : d[i - 1];
		d[i - 1] = 0;
	}

	printf("core value (max truncated degree) = %u\n", max);

	g->adj = (unsigned*)malloc(el->e * sizeof(unsigned));

	for (i = 0; i < el->e; i++) {
		g->adj[g->cd[el->edges[i].s] + d[el->edges[i].s]++] = el->edges[i].t;
	}

	free(d);
	g->core = max;
	g->n = el->n;
	return g;
}


subgraph* allocsub(graph* g, unsigned k) {
	unsigned i;
	subgraph* sg = (subgraph*)malloc(sizeof(subgraph));
	sg->n = (unsigned*)calloc(k, sizeof(unsigned));
	sg->d = (unsigned**)malloc(k * sizeof(unsigned*));
	sg->nodes = (unsigned**)malloc(k * sizeof(unsigned*));
	sg->adj = (unsigned*)malloc(g->core * (g->core / block_size + 1) * sizeof(unsigned));
	for (i = 1; i < k; i++) {
		sg->d[i] = (unsigned*)malloc(g->core * sizeof(unsigned));
		sg->nodes[i] = (unsigned*)malloc(g->core * sizeof(unsigned));
	}
	sg->core = g->core;
	return sg;
}

void mksub(graph* g, unsigned u, subgraph* sg, unsigned k) {
	unsigned i, j, l, v, w, s, t, id_num, id_bit;

	static unsigned* old0 = NULL, * new0 = NULL;//to improve
#pragma omp threadprivate(new0,old0)

	if (old0 == NULL) {
		new0 = (unsigned*)malloc(g->n * sizeof(unsigned));
		old0 = (unsigned*)malloc(g->core * sizeof(unsigned));
		for (i = 0; i < g->n; i++) {
			new0[i] = -1;
		}
	}

	j = g->cd[u + 1] - g->cd[u];
	sg->n[k - 1] = j;
	sg->block_num = j / block_size + (j % block_size != 0);
	memset(sg->nodes[k - 1], 0, sizeof(unsigned) * (sg->block_num));
	j = 0;
	for (i = g->cd[u]; i < g->cd[u + 1]; i++) {
		v = g->adj[i];
		new0[v] = j;
		old0[j] = v;
		sg->d[k - 1][j] = 0;//new degrees
		j++;
	}

	/*for (i = 0; i < sg->n[k - 1]; i++) {
		id_num = i / block_size, id_bit = i % block_size;
		sg->nodes[k - 1][id_num] |= (1 << id_bit);
	}*/
	j = sg->n[k - 1] / block_size;
	for (i = 0; i < j; i++) { 
		sg->nodes[k - 1][i] = (1 << block_size) - 1;
	}
	j = sg->n[k - 1] % block_size;

	if (j > 0)
		sg->nodes[k - 1][sg->block_num - 1] = (1 << j) - 1;
	j = sg->n[k - 1];
	unsigned* d0 = (unsigned*)calloc(j, sizeof(unsigned));
	memset(sg->adj, 0, sizeof(unsigned) * sg->core * sg->block_num);
	for (i = 0; i < sg->n[k - 1]; i++) {//reodering adjacency list and computing new degrees
		v = old0[i];
		unsigned id_numi = i / block_size, id_biti = i % block_size;
		for (l = g->cd[v]; l < g->cd[v + 1]; l++) {
			w = g->adj[l];

			j = new0[w];

			if (j != (unsigned)-1) {
				unsigned id_numj = j / block_size, id_bitj = j % block_size;
				sg->adj[sg->block_num * i + id_numj] |= (1 << id_bitj); //compress the induced subgraphs
				sg->adj[sg->block_num * j + id_numi] |= (1 << id_biti);
				d0[i]++;
				d0[j]++;

			}
		}
	}
	unsigned* C = (unsigned*)calloc(sg->n[k - 1], sizeof(unsigned));
	int* color = (int*)malloc(sg->n[k - 1] * sizeof(int));
	unsigned* Index = (unsigned*)malloc(sg->n[k - 1] * sizeof(unsigned));
	iddegree* ig;
	ig = (iddegree*)malloc(sg->n[k - 1] * sizeof(iddegree));
	for (i = 0; i < sg->n[k - 1]; i++)
	{
		color[i] = -1;
		ig[i].id = i;
		ig[i].degree = d0[i];
	}
	qsort(ig, sg->n[k - 1], sizeof(ig[0]), cmp);

	for (i = 0; i < sg->n[k - 1]; i++)
		Index[ig[i].id] = i;


	//color ordering
	color[0] = 0;
	unsigned colorNum = 0;

	for (unsigned i = 1; i < sg->n[k - 1]; i++)
	{
		int tmpid = ig[i].id;
		for (unsigned h = 0; h < sg->block_num; h++) {
			v = sg->adj[sg->block_num * tmpid + h];
			//unsigned tmp_mask_cnt = mask_cnt[v];
			unsigned tmp_mask_cnt = mask[v][0];
			//for (int j = 0; j < tmp_mask_cnt; j++) {
			for (unsigned j = 1; j <= tmp_mask_cnt; j++) {
				int now = Index[block_size * h + mask[v][j]];
				if (color[now] != -1)
					C[color[now]] = 1;
			}
		}
		for (unsigned j = 0; j < ig[0].degree + 1; j++)
			if (C[j] == 0)
			{
				color[i] = j;
				colorNum = j > colorNum ? j : colorNum;
				break;
			}
		for (unsigned h = 0; h < sg->block_num; h++) {
			v = sg->adj[sg->block_num * tmpid + h];
			//unsigned tmp_mask_cnt = mask_cnt[v];
			unsigned tmp_mask_cnt = mask[v][0];
			//for (int j = 0; j < tmp_mask_cnt; j++) {
			for (unsigned j = 1; j <= tmp_mask_cnt; j++) {
				int now = Index[block_size * h + mask[v][j]];
				if (color[now] != -1)
					C[color[now]] = 0;
			}
		}
	}

	sg->color = (unsigned*)malloc(sg->n[k - 1] * sizeof(unsigned));
	memset(sg->adj, 0, sizeof(unsigned) * sg->n[k - 1] * sg->block_num);
	for (unsigned i = 0; i < sg->n[k - 1]; i++)
	{
		sg->color[i] = color[Index[i]];
	}

	for (i = 0; i < sg->n[k - 1]; i++) {
		v = old0[i];
		for (l = g->cd[v]; l < g->cd[v + 1]; l++) {
			w = g->adj[l];
			j = new0[w];
			if (j != (unsigned)-1) {

				if (color[Index[i]] > color[Index[j]])
				{
					s = i, t = j;
				}
				else
				{
					s = j, t = i;
				}
				id_num = t / block_size, id_bit = t % block_size;
				sg->adj[sg->block_num * s + id_num] |= (1 << id_bit);
			}
		}
	}


	for (i = g->cd[u]; i < g->cd[u + 1]; i++) {
		new0[g->adj[i]] = -1;
	}
}
void kclique_thread(unsigned l, subgraph* sg, unsigned long long* n) {
	unsigned i, j, u, v, cnt;

	if (l == 2) {
		for (i = 0; i < sg->block_num; i++) {
			u = sg->nodes[2][i];

			cnt = mask[u][0];
			//cnt = mask_cnt[u];
			for (j = 1; j <= cnt; j++) {
				v = mask[u][j] + block_size * i;
				(*n) += block_intersect_count(sg->adj + sg->block_num * v, sg->nodes[2], sg->block_num);
			}
		}
		return;
	}

	if (l > sg->n[l])
		return;
	for (i = 0; i < sg->block_num; i++) {
		u = sg->nodes[l][i];
		cnt = mask[u][0];
		//cnt = mask_cnt[u];
		for (j = 1; j <= cnt; j++) {
			v = mask[u][j] + block_size * i;
			if (sg->color[v] < l - 1)
				continue;
			sg->n[l - 1] = block_intersect(sg->adj + sg->block_num * v, sg->nodes[l], sg->nodes[l - 1], sg->block_num);
			kclique_thread(l - 1, sg, n);
		}
	}
}

unsigned long long kclique_main(unsigned k, graph* g) {
	unsigned u;
	unsigned long long n = 0;
	subgraph* sg;
#pragma omp parallel private(sg,u) shared(mask) reduction(+:n)
	{
		sg = allocsub(g, k);
#pragma omp for schedule(dynamic, 1) nowait
		for (u = 0; u < g->n; u++) {
			if (g->cd[u + 1] - g->cd[u] >= k - 1) {
				mksub(g, u, sg, k);
				kclique_thread(k - 1, sg, &n);
			}

		}

	}
	return n;
}

int main(int argc, char** argv) {
	struct timeval time_start;
	struct timeval time_end;
	edgelist* el;
	graph* g;
	unsigned k = atoi(argv[2]);
	unsigned long long n;
	omp_set_num_threads(atoi(argv[1]));

	printf("Reading edgelist from file %s\n", argv[3]);
	gettimeofday(&time_start, NULL);
	el = readedgelist(argv[3]);
	printf("Number of nodes = %u\n", el->n);
	printf("Number of edges = %u\n", el->e);
	gettimeofday(&time_end, NULL);
	double read_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;
	printf("- Read Time = %fs\n", read_time / 1000);
	ord_core(el,k);
	relabel(el);
	g = mkgraph(el);
	free_edgelist(el);
	gettimeofday(&time_end, NULL);
	printf("Iterate over all cliques\n");
	gettimeofday(&time_start, NULL);
	n = kclique_main(k, g);
	gettimeofday(&time_end, NULL);
	double list_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;
	printf("Number of %u-cliques: %llu\n", k, n);
	printf("- List Time = %fs\n", list_time / 1000);
	free_graph(g);
	return 0;
}
