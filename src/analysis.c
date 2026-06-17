#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdarg.h>
#include "nucesGraph.h"

struct label_index {
	int label;
	int index;
};

struct graph_adj {
	int n;
	int m;
	int *out_start;
	int *out;
	int *in_start;
	int *in;
	int *labels;
	struct label_index *map;
};

static int cmp_label_index(const void *a, const void *b)
{
	const struct label_index *la = (const struct label_index *)a;
	const struct label_index *lb = (const struct label_index *)b;
	if (la->label < lb->label) return -1;
	if (la->label > lb->label) return 1;
	return 0;
}

static int map_label_to_index(const struct label_index *map, int n, int label)
{
	int lo = 0;
	int hi = n - 1;
	while (lo <= hi) {
		int mid = lo + (hi - lo) / 2;
		if (map[mid].label == label) {
			return map[mid].index;
		}
		if (map[mid].label < label) {
			lo = mid + 1;
		} else {
			hi = mid - 1;
		}
	}
	return -1;
}

static int build_label_map(const struct nGraph *G, int **labels_out, struct label_index **map_out)
{
	int n = G->V->count;
	int *labels = NULL;
	struct label_index *map = NULL;
	struct vertex *tmp = NULL;
	int i = 0;

	if (n <= 0) {
		return -1;
	}

	labels = (int *)calloc((size_t)n, sizeof(int));
	map = (struct label_index *)calloc((size_t)n, sizeof(struct label_index));
	if (labels == NULL || map == NULL) {
		free(labels);
		free(map);
		return -1;
	}

	tmp = G->V->head;
	while (tmp != NULL && i < n) {
		labels[i] = tmp->label;
		map[i].label = tmp->label;
		map[i].index = i;
		i++;
		tmp = tmp->next;
	}

	if (i != n) {
		free(labels);
		free(map);
		return -1;
	}

	qsort(map, (size_t)n, sizeof(struct label_index), cmp_label_index);

	*labels_out = labels;
	*map_out = map;
	return n;
}

static void free_graph_adj(struct graph_adj *adj)
{
	if (adj == NULL) {
		return;
	}
	free(adj->out_start);
	free(adj->out);
	free(adj->in_start);
	free(adj->in);
	free(adj->labels);
	free(adj->map);
	memset(adj, 0, sizeof(*adj));
}

static int build_adjacency(const struct nGraph *G, int undirected_mode, struct graph_adj *adj)
{
	int n = 0;
	int *out_deg = NULL;
	int *in_deg = NULL;
	int *out_cursor = NULL;
	int *in_cursor = NULL;
	struct edge *e = NULL;
	int i = 0;

	memset(adj, 0, sizeof(*adj));

	n = build_label_map(G, &adj->labels, &adj->map);
	if (n <= 0) {
		return -1;
	}
	adj->n = n;

	out_deg = (int *)calloc((size_t)n, sizeof(int));
	in_deg = (int *)calloc((size_t)n, sizeof(int));
	if (out_deg == NULL || in_deg == NULL) {
		free(out_deg);
		free(in_deg);
		free_graph_adj(adj);
		return -1;
	}

	e = G->E->head;
	while (e != NULL) {
		int u = map_label_to_index(adj->map, n, e->head);
		int v = map_label_to_index(adj->map, n, e->tail);
		if (u < 0 || v < 0) {
			free(out_deg);
			free(in_deg);
			free_graph_adj(adj);
			return -1;
		}
		if (undirected_mode || e->directed == 0) {
			out_deg[u]++;
			out_deg[v]++;
			in_deg[u]++;
			in_deg[v]++;
		} else {
			out_deg[u]++;
			in_deg[v]++;
		}
		e = e->next;
	}

	adj->out_start = (int *)calloc((size_t)(n + 1), sizeof(int));
	adj->in_start = (int *)calloc((size_t)(n + 1), sizeof(int));
	if (adj->out_start == NULL || adj->in_start == NULL) {
		free(out_deg);
		free(in_deg);
		free_graph_adj(adj);
		return -1;
	}

	for (i = 0; i < n; i++) {
		adj->out_start[i + 1] = adj->out_start[i] + out_deg[i];
		adj->in_start[i + 1] = adj->in_start[i] + in_deg[i];
	}
	adj->m = adj->out_start[n];
	if (adj->m > 0) {
		adj->out = (int *)calloc((size_t)adj->m, sizeof(int));
	}
	if (adj->in_start[n] > 0) {
		adj->in = (int *)calloc((size_t)adj->in_start[n], sizeof(int));
	}
	if ((adj->m > 0 && adj->out == NULL) || (adj->in_start[n] > 0 && adj->in == NULL)) {
		free(out_deg);
		free(in_deg);
		free_graph_adj(adj);
		return -1;
	}

	out_cursor = (int *)calloc((size_t)n, sizeof(int));
	in_cursor = (int *)calloc((size_t)n, sizeof(int));
	if (out_cursor == NULL || in_cursor == NULL) {
		free(out_deg);
		free(in_deg);
		free(out_cursor);
		free(in_cursor);
		free_graph_adj(adj);
		return -1;
	}
	for (i = 0; i < n; i++) {
		out_cursor[i] = adj->out_start[i];
		in_cursor[i] = adj->in_start[i];
	}

	e = G->E->head;
	while (e != NULL) {
		int u = map_label_to_index(adj->map, n, e->head);
		int v = map_label_to_index(adj->map, n, e->tail);
		if (undirected_mode || e->directed == 0) {
			adj->out[out_cursor[u]++] = v;
			adj->out[out_cursor[v]++] = u;
			adj->in[in_cursor[u]++] = v;
			adj->in[in_cursor[v]++] = u;
		} else {
			adj->out[out_cursor[u]++] = v;
			adj->in[in_cursor[v]++] = u;
		}
		e = e->next;
	}

	free(out_deg);
	free(in_deg);
	free(out_cursor);
	free(in_cursor);
	return 0;
}

int getVertexLabels(struct nGraph *G, int *labels, int max_labels)
{
	int i = 0;
	if (G == NULL || G->V == NULL) {
		return -1;
	}
	if (labels == NULL || max_labels <= 0) {
		return -1;
	}
	struct vertex *tmp = G->V->head;
	while (tmp != NULL && i < max_labels) {
		labels[i++] = tmp->label;
		tmp = tmp->next;
	}
	return i;
}

int graphNodeCount(struct nGraph *G)
{
	if (G == NULL || G->V == NULL) return -1;
	return G->V->count;
}

int graphEdgeCount(struct nGraph *G)
{
	if (G == NULL || G->E == NULL) return -1;
	return G->E->count;
}

static int count_weak_components(const struct graph_adj *adj)
{
	int n = adj->n;
	int *visited = (int *)calloc((size_t)n, sizeof(int));
	int *queue = (int *)calloc((size_t)n, sizeof(int));
	int components = 0;
	int i;

	if (visited == NULL || queue == NULL) {
		free(visited);
		free(queue);
		return -1;
	}

	for (i = 0; i < n; i++) {
		if (visited[i]) continue;
		components++;
		int qh = 0;
		int qt = 0;
		visited[i] = 1;
		queue[qt++] = i;
		while (qh < qt) {
			int v = queue[qh++];
			int start = adj->out_start[v];
			int end = adj->out_start[v + 1];
			int j;
			for (j = start; j < end; j++) {
				int w = adj->out[j];
				if (!visited[w]) {
					visited[w] = 1;
					queue[qt++] = w;
				}
			}
		}
	}

	free(visited);
	free(queue);
	return components;
}

int cyclomaticComplexity(struct nGraph *G)
{
	struct graph_adj adj;
	int components = 0;
	int complexity = 0;
	int n = 0;

	if (G == NULL || G->V == NULL || G->E == NULL) {
		return -1;
	}
	if (build_adjacency(G, 1, &adj) != 0) {
		return -1;
	}
	components = count_weak_components(&adj);
	n = adj.n;

	if (components < 0) {
		free_graph_adj(&adj);
		return -1;
	}
	/* M = E - N + 2P: count each edge once regardless of directed flag */
	complexity = G->E->count - n + 2 * components;

	free_graph_adj(&adj);
	return complexity;
}

int cfgDfsOrder(struct nGraph *G, int start_label, int *order, int max_order)
{
	struct graph_adj adj;
	int *visited = NULL;
	int *stack = NULL;
	int count = 0;
	int start = -1;
	if (G == NULL || G->V == NULL || G->E == NULL) {
		return -1;
	}
	if (order == NULL || max_order <= 0) {
		return -1;
	}

	if (build_adjacency(G, 0, &adj) != 0) {
		return -1;
	}
	start = map_label_to_index(adj.map, adj.n, start_label);
	if (start < 0) {
		free_graph_adj(&adj);
		return -1;
	}

	visited = (int *)calloc((size_t)adj.n, sizeof(int));
	stack = (int *)calloc((size_t)adj.n, sizeof(int));
	if (visited == NULL || stack == NULL) {
		free(visited);
		free(stack);
		free_graph_adj(&adj);
		return -1;
	}

	/* Mark on push so each node is pushed at most once — stack bounded by n */
	int top = 0;
	visited[start] = 1;
	stack[top++] = start;
	while (top > 0) {
		int v = stack[--top];
		if (count < max_order) {
			order[count++] = adj.labels[v];
		}

		int s = adj.out_start[v];
		int e = adj.out_start[v + 1];
		int i;
		for (i = e - 1; i >= s; i--) {
			int w = adj.out[i];
			if (!visited[w]) {
				visited[w] = 1;
				stack[top++] = w;
			}
		}
	}

	free(visited);
	free(stack);
	free_graph_adj(&adj);
	return count;
}

int cfgBfsOrder(struct nGraph *G, int start_label, int *order, int max_order)
{
	struct graph_adj adj;
	int *visited = NULL;
	int *queue = NULL;
	int count = 0;
	int start = -1;
	if (G == NULL || G->V == NULL || G->E == NULL) {
		return -1;
	}
	if (order == NULL || max_order <= 0) {
		return -1;
	}

	if (build_adjacency(G, 0, &adj) != 0) {
		return -1;
	}
	start = map_label_to_index(adj.map, adj.n, start_label);
	if (start < 0) {
		free_graph_adj(&adj);
		return -1;
	}

	visited = (int *)calloc((size_t)adj.n, sizeof(int));
	queue = (int *)calloc((size_t)adj.n, sizeof(int));
	if (visited == NULL || queue == NULL) {
		free(visited);
		free(queue);
		free_graph_adj(&adj);
		return -1;
	}

	int qh = 0;
	int qt = 0;
	visited[start] = 1;
	queue[qt++] = start;
	while (qh < qt) {
		int v = queue[qh++];
		if (count < max_order) {
			order[count++] = adj.labels[v];
		}

		int s = adj.out_start[v];
		int e = adj.out_start[v + 1];
		int i;
		for (i = s; i < e; i++) {
			int w = adj.out[i];
			if (!visited[w]) {
				visited[w] = 1;
				queue[qt++] = w;
			}
		}
	}

	free(visited);
	free(queue);
	free_graph_adj(&adj);
	return count;
}

static int dom_word_count(int n)
{
	int bits = (int)(sizeof(unsigned long) * CHAR_BIT);
	return (n + bits - 1) / bits;
}

static void dom_set_all(unsigned long *set, int words)
{
	int i;
	for (i = 0; i < words; i++) {
		set[i] = ~0UL;
	}
}

static void dom_set_zero(unsigned long *set, int words)
{
	int i;
	for (i = 0; i < words; i++) {
		set[i] = 0UL;
	}
}

static void dom_set_bit(unsigned long *set, int idx)
{
	int bits = (int)(sizeof(unsigned long) * CHAR_BIT);
	set[idx / bits] |= (1UL << (unsigned int)(idx % bits));
}

static int dom_test_bit(const unsigned long *set, int idx)
{
	int bits = (int)(sizeof(unsigned long) * CHAR_BIT);
	return (set[idx / bits] >> (unsigned int)(idx % bits)) & 1UL;
}

static void dom_set_and(unsigned long *dst, const unsigned long *src, int words)
{
	int i;
	for (i = 0; i < words; i++) {
		dst[i] &= src[i];
	}
}

static int dom_set_equal(const unsigned long *a, const unsigned long *b, int words)
{
	int i;
	for (i = 0; i < words; i++) {
		if (a[i] != b[i]) return 0;
	}
	return 1;
}

int dominatorWordCount(int n)
{
	return dom_word_count(n);
}

static int compute_dominators_internal(const struct graph_adj *adj, int start, unsigned long *dom, int words, int use_preds)
{
	int n = adj->n;
	int changed = 1;
	int v;
	unsigned long *newset = NULL;

	newset = (unsigned long *)calloc((size_t)words, sizeof(unsigned long));
	if (newset == NULL) {
		return -1;
	}

	for (v = 0; v < n; v++) {
		dom_set_all(&dom[v * words], words);
	}
	dom_set_zero(&dom[start * words], words);
	dom_set_bit(&dom[start * words], start);

	while (changed) {
		changed = 0;
		for (v = 0; v < n; v++) {
			unsigned long *cur = &dom[v * words];
			int has_pred = 0;
			int i;
			if (v == start) {
				continue;
			}

			dom_set_all(newset, words);
			if (use_preds) {
				int s = adj->in_start[v];
				int e = adj->in_start[v + 1];
				for (i = s; i < e; i++) {
					int p = adj->in[i];
					dom_set_and(newset, &dom[p * words], words);
					has_pred = 1;
				}
			} else {
				int s = adj->out_start[v];
				int e = adj->out_start[v + 1];
				for (i = s; i < e; i++) {
					int p = adj->out[i];
					dom_set_and(newset, &dom[p * words], words);
					has_pred = 1;
				}
			}

			if (!has_pred) {
				dom_set_zero(newset, words);
			}
			dom_set_bit(newset, v);

			if (!dom_set_equal(cur, newset, words)) {
				memcpy(cur, newset, (size_t)words * sizeof(unsigned long));
				changed = 1;
			}
		}
	}

	free(newset);

	return 0;
}

int computeDominators(struct nGraph *G, int start_label, unsigned long *dom, int word_count)
{
	struct graph_adj adj;
	int start = -1;
	int rc = 0;

	if (G == NULL || G->V == NULL || G->E == NULL || dom == NULL) {
		return -1;
	}
	if (build_adjacency(G, 0, &adj) != 0) {
		return -1;
	}
	start = map_label_to_index(adj.map, adj.n, start_label);
	if (start < 0) {
		free_graph_adj(&adj);
		return -1;
	}
	if (word_count < dom_word_count(adj.n)) {
		free_graph_adj(&adj);
		return -1;
	}

	rc = compute_dominators_internal(&adj, start, dom, word_count, 1);
	free_graph_adj(&adj);
	return rc;
}

int computePostDominators(struct nGraph *G, int exit_label, unsigned long *pdom, int word_count)
{
	struct graph_adj adj;
	int exit_index = -1;
	int rc = 0;

	if (G == NULL || G->V == NULL || G->E == NULL || pdom == NULL) {
		return -1;
	}
	if (build_adjacency(G, 0, &adj) != 0) {
		return -1;
	}
	exit_index = map_label_to_index(adj.map, adj.n, exit_label);
	if (exit_index < 0) {
		free_graph_adj(&adj);
		return -1;
	}
	if (word_count < dom_word_count(adj.n)) {
		free_graph_adj(&adj);
		return -1;
	}

	rc = compute_dominators_internal(&adj, exit_index, pdom, word_count, 0);
	free_graph_adj(&adj);
	return rc;
}

static int dominates(const unsigned long *dom, int words, int a, int b)
{
	return dom_test_bit(&dom[b * words], a);
}

/* Count set bits in a bitvector — used for O(n^2*words) idom extraction */
static int dom_popcount(const unsigned long *set, int words)
{
	int count = 0;
	int i;
	for (i = 0; i < words; i++) {
		unsigned long v = set[i];
		while (v) {
			count += (int)(v & 1UL);
			v >>= 1;
		}
	}
	return count;
}

int computeImmediateDominators(struct nGraph *G, int start_label, int *idom)
{
	struct graph_adj adj;
	int words = 0;
	unsigned long *dom = NULL;
	int start = -1;
	int n = 0;
	int v, d;

	if (G == NULL || G->V == NULL || G->E == NULL || idom == NULL) {
		return -1;
	}
	if (build_adjacency(G, 0, &adj) != 0) {
		return -1;
	}
	start = map_label_to_index(adj.map, adj.n, start_label);
	if (start < 0) {
		free_graph_adj(&adj);
		return -1;
	}

	n = adj.n;
	words = dom_word_count(n);
	dom = (unsigned long *)calloc((size_t)n * (size_t)words, sizeof(unsigned long));
	if (dom == NULL) {
		free_graph_adj(&adj);
		return -1;
	}

	if (compute_dominators_internal(&adj, start, dom, words, 1) != 0) {
		free(dom);
		free_graph_adj(&adj);
		return -1;
	}

	/*
	 * idom(v) is the unique dominator d of v where |dom(d)| = |dom(v)| - 1.
	 * Using popcount on the bitset reduces the triple-loop O(n^3) to O(n^2*words).
	 */
	for (v = 0; v < n; v++) {
		int v_count;
		idom[v] = -1;
		if (v == start) {
			idom[v] = adj.labels[start];
			continue;
		}
		v_count = dom_popcount(&dom[v * words], words);
		for (d = 0; d < n; d++) {
			if (d == v) continue;
			if (!dominates(dom, words, d, v)) continue;
			if (dom_popcount(&dom[d * words], words) == v_count - 1) {
				idom[v] = adj.labels[d];
				break;
			}
		}
	}

	free(dom);
	free_graph_adj(&adj);
	return 0;
}

int computeImmediatePostDominators(struct nGraph *G, int exit_label, int *ipdom)
{
	struct graph_adj adj;
	int words = 0;
	unsigned long *pdom = NULL;
	int exit_index = -1;
	int n = 0;
	int v, d;

	if (G == NULL || G->V == NULL || G->E == NULL || ipdom == NULL) {
		return -1;
	}
	if (build_adjacency(G, 0, &adj) != 0) {
		return -1;
	}
	exit_index = map_label_to_index(adj.map, adj.n, exit_label);
	if (exit_index < 0) {
		free_graph_adj(&adj);
		return -1;
	}

	n = adj.n;
	words = dom_word_count(n);
	pdom = (unsigned long *)calloc((size_t)n * (size_t)words, sizeof(unsigned long));
	if (pdom == NULL) {
		free_graph_adj(&adj);
		return -1;
	}

	if (compute_dominators_internal(&adj, exit_index, pdom, words, 0) != 0) {
		free(pdom);
		free_graph_adj(&adj);
		return -1;
	}

	/* Same popcount-based O(n^2*words) approach as computeImmediateDominators */
	for (v = 0; v < n; v++) {
		int v_count;
		ipdom[v] = -1;
		if (v == exit_index) {
			ipdom[v] = adj.labels[exit_index];
			continue;
		}
		v_count = dom_popcount(&pdom[v * words], words);
		for (d = 0; d < n; d++) {
			if (d == v) continue;
			if (!dominates(pdom, words, d, v)) continue;
			if (dom_popcount(&pdom[d * words], words) == v_count - 1) {
				ipdom[v] = adj.labels[d];
				break;
			}
		}
	}

	free(pdom);
	free_graph_adj(&adj);
	return 0;
}

/*
 * Iterative Tarjan SCC — replaces the recursive version to avoid stack
 * overflow on graphs with O(n) chain depth.  Each frame stores the node
 * being processed and the next edge index (ei) to resume from, exactly
 * mirroring what the call stack held in the recursive formulation.
 */
struct scc_frame {
	int v;
	int ei;
};

static void scc_iterative(const struct graph_adj *adj, int root,
	int *index_arr, int *low, int *stk, int *onstack,
	int *sp, int *time_val, int *scc_id, int *scc_count,
	struct scc_frame *frames)
{
	int fsp = 0;

	index_arr[root] = low[root] = (*time_val)++;
	stk[(*sp)++] = root;
	onstack[root] = 1;
	frames[fsp].v = root;
	frames[fsp].ei = adj->out_start[root];
	fsp++;

	while (fsp > 0) {
		struct scc_frame *f = &frames[fsp - 1];
		int v = f->v;

		if (f->ei < adj->out_start[v + 1]) {
			int w = adj->out[f->ei++];
			if (index_arr[w] < 0) {
				/* push new frame — equivalent to recursive call */
				index_arr[w] = low[w] = (*time_val)++;
				stk[(*sp)++] = w;
				onstack[w] = 1;
				frames[fsp].v = w;
				frames[fsp].ei = adj->out_start[w];
				fsp++;
			} else if (onstack[w] && index_arr[w] < low[v]) {
				low[v] = index_arr[w];
			}
		} else {
			/* all edges processed — equivalent to returning from recursion */
			fsp--;
			if (fsp > 0) {
				int parent = frames[fsp - 1].v;
				if (low[v] < low[parent]) {
					low[parent] = low[v];
				}
			}
			if (low[v] == index_arr[v]) {
				int w;
				do {
					w = stk[--(*sp)];
					onstack[w] = 0;
					scc_id[w] = *scc_count;
				} while (w != v);
				(*scc_count)++;
			}
		}
	}
}

int computeSCCs(struct nGraph *G, int *scc_id)
{
	struct graph_adj adj;
	int *index_arr = NULL;
	int *low = NULL;
	int *stack = NULL;
	int *onstack = NULL;
	struct scc_frame *frames = NULL;
	int sp = 0;
	int time_val = 0;
	int scc_count = 0;
	int v;

	if (G == NULL || G->V == NULL || G->E == NULL || scc_id == NULL) {
		return -1;
	}
	if (build_adjacency(G, 0, &adj) != 0) {
		return -1;
	}

	index_arr = (int *)calloc((size_t)adj.n, sizeof(int));
	low       = (int *)calloc((size_t)adj.n, sizeof(int));
	stack     = (int *)calloc((size_t)adj.n, sizeof(int));
	onstack   = (int *)calloc((size_t)adj.n, sizeof(int));
	frames    = (struct scc_frame *)calloc((size_t)adj.n, sizeof(struct scc_frame));
	if (index_arr == NULL || low == NULL || stack == NULL || onstack == NULL || frames == NULL) {
		free(index_arr);
		free(low);
		free(stack);
		free(onstack);
		free(frames);
		free_graph_adj(&adj);
		return -1;
	}

	for (v = 0; v < adj.n; v++) {
		index_arr[v] = -1;
		low[v]       = -1;
		scc_id[v]    = -1;
	}

	for (v = 0; v < adj.n; v++) {
		if (index_arr[v] < 0) {
			scc_iterative(&adj, v, index_arr, low, stack, onstack,
				&sp, &time_val, scc_id, &scc_count, frames);
		}
	}

	free(index_arr);
	free(low);
	free(stack);
	free(onstack);
	free(frames);
	free_graph_adj(&adj);
	return scc_count;
}

int computeLoopNestingDepth(struct nGraph *G, int start_label, int *depth_out)
{
	struct graph_adj adj;
	int words = 0;
	unsigned long *dom = NULL;
	int start = -1;
	int n = 0;
	int u;
	int max_depth = 0;

	if (G == NULL || G->V == NULL || G->E == NULL || depth_out == NULL) {
		return -1;
	}
	if (build_adjacency(G, 0, &adj) != 0) {
		return -1;
	}
	start = map_label_to_index(adj.map, adj.n, start_label);
	if (start < 0) {
		free_graph_adj(&adj);
		return -1;
	}

	n = adj.n;
	words = dom_word_count(n);
	dom = (unsigned long *)calloc((size_t)n * (size_t)words, sizeof(unsigned long));
	if (dom == NULL) {
		free_graph_adj(&adj);
		return -1;
	}
	if (compute_dominators_internal(&adj, start, dom, words, 1) != 0) {
		free(dom);
		free_graph_adj(&adj);
		return -1;
	}

	for (u = 0; u < n; u++) {
		depth_out[u] = 0;
	}

	/* Allocate once outside the loop — reused for every back edge via memset */
	{
		int *loop_stack = (int *)calloc((size_t)n, sizeof(int));
		int *in_loop    = (int *)calloc((size_t)n, sizeof(int));
		if (loop_stack == NULL || in_loop == NULL) {
			free(loop_stack);
			free(in_loop);
			free(dom);
			free_graph_adj(&adj);
			return -1;
		}

		for (u = 0; u < n; u++) {
			int i;
			for (i = adj.out_start[u]; i < adj.out_start[u + 1]; i++) {
				int v = adj.out[i];
				if (dominates(dom, words, v, u)) {
					int top = 0;
					int x;
					/* Reset for this back-edge's loop body — O(n) but no alloc */
					memset(in_loop, 0, (size_t)n * sizeof(int));
					in_loop[v] = 1;
					in_loop[u] = 1;
					loop_stack[top++] = u;
					while (top > 0) {
						x = loop_stack[--top];
						int p;
						for (p = adj.in_start[x]; p < adj.in_start[x + 1]; p++) {
							int pred = adj.in[p];
							if (!in_loop[pred] && dominates(dom, words, v, pred)) {
								in_loop[pred] = 1;
								loop_stack[top++] = pred;
							}
						}
					}

					for (x = 0; x < n; x++) {
						if (in_loop[x]) {
							depth_out[x]++;
							if (depth_out[x] > max_depth) {
								max_depth = depth_out[x];
							}
						}
					}
				}
			}
		}

		free(loop_stack);
		free(in_loop);
	}

	free(dom);
	free_graph_adj(&adj);
	return max_depth;
}

int sliceForward(struct nGraph *dep, int start_label, int *mark, int mark_len)
{
	struct graph_adj adj;
	int *stack = NULL;
	int count = 0;
	int start = -1;
	if (dep == NULL || dep->V == NULL || dep->E == NULL) {
		return -1;
	}
	if (mark == NULL || mark_len <= 0) {
		return -1;
	}

	if (build_adjacency(dep, 0, &adj) != 0) {
		return -1;
	}
	if (mark_len < adj.n) {
		free_graph_adj(&adj);
		return -1;
	}
	start = map_label_to_index(adj.map, adj.n, start_label);
	if (start < 0) {
		free_graph_adj(&adj);
		return -1;
	}

	stack = (int *)calloc((size_t)adj.n, sizeof(int));
	if (stack == NULL) {
		free_graph_adj(&adj);
		return -1;
	}

	/* Mark on push — each node pushed at most once, stack bounded by n */
	memset(mark, 0, (size_t)adj.n * sizeof(int));
	int top = 0;
	mark[start] = 1;
	count = 1;
	stack[top++] = start;
	while (top > 0) {
		int v = stack[--top];
		int i;
		for (i = adj.out_start[v]; i < adj.out_start[v + 1]; i++) {
			int w = adj.out[i];
			if (!mark[w]) {
				mark[w] = 1;
				count++;
				stack[top++] = w;
			}
		}
	}

	free(stack);
	free_graph_adj(&adj);
	return count;
}

int sliceBackward(struct nGraph *dep, int start_label, int *mark, int mark_len)
{
	struct graph_adj adj;
	int *stack = NULL;
	int count = 0;
	int start = -1;
	if (dep == NULL || dep->V == NULL || dep->E == NULL) {
		return -1;
	}
	if (mark == NULL || mark_len <= 0) {
		return -1;
	}

	if (build_adjacency(dep, 0, &adj) != 0) {
		return -1;
	}
	if (mark_len < adj.n) {
		free_graph_adj(&adj);
		return -1;
	}
	start = map_label_to_index(adj.map, adj.n, start_label);
	if (start < 0) {
		free_graph_adj(&adj);
		return -1;
	}

	stack = (int *)calloc((size_t)adj.n, sizeof(int));
	if (stack == NULL) {
		free_graph_adj(&adj);
		return -1;
	}

	/* Mark on push — each node pushed at most once, stack bounded by n */
	memset(mark, 0, (size_t)adj.n * sizeof(int));
	int top = 0;
	mark[start] = 1;
	count = 1;
	stack[top++] = start;
	while (top > 0) {
		int v = stack[--top];
		int i;
		for (i = adj.in_start[v]; i < adj.in_start[v + 1]; i++) {
			int w = adj.in[i];
			if (!mark[w]) {
				mark[w] = 1;
				count++;
				stack[top++] = w;
			}
		}
	}

	free(stack);
	free_graph_adj(&adj);
	return count;
}

static int appendf(char *buf, size_t cap, size_t *len, const char *fmt, ...)
{
	va_list ap;
	int n;

	if (*len >= cap) {
		return -1;
	}

	va_start(ap, fmt);
	n = vsnprintf(buf + *len, cap - *len, fmt, ap);
	va_end(ap);

	if (n < 0 || (size_t)n >= cap - *len) {
		return -1;
	}
	*len += (size_t)n;
	return 0;
}

static void format_cell(char *dst, size_t size, int value)
{
	if (value < 0) {
		snprintf(dst, size, "NA");
	} else {
		snprintf(dst, size, "%d", value);
	}
}

struct analysis_data {
	/* graph-level metrics */
	int    n;
	int    complexity;      /* cyclomatic: M = E - N + 2P               */
	int    scc_count;       /* strongly connected components             */
	int    max_depth;       /* maximum natural loop nesting depth        */
	double density;         /* E / (N*(N-1)) — directed graph density    */
	int    back_edges;      /* back-edge count = natural loop count      */
	int    predicate_count; /* nodes with fan_out >= 2 (decision points) */
	int    join_count;      /* nodes with fan_in  >= 2 (merge points)    */
	int    dead_count;      /* nodes unreachable from entry              */
	int    longest_path;    /* max BFS depth from entry (path proxy)     */

	/* per-node arrays (indexed by vertex-list order, same as labels[]) */
	int   *labels;
	int   *idom;      /* label of immediate dominator, -1 if none   */
	int   *ipdom;     /* label of immediate post-dominator           */
	int   *scc_id;    /* SCC identifier                              */
	int   *depth;     /* natural loop nesting depth                  */
	int   *fan_in;    /* in-degree                                   */
	int   *fan_out;   /* out-degree                                  */
	int   *dom_depth; /* depth in dominator tree (0=entry, -1=N/A)  */
	int   *reachable; /* 1 if reachable from entry, 0 if dead        */
};

static void free_analysis_data(struct analysis_data *data)
{
	free(data->labels);
	free(data->idom);
	free(data->ipdom);
	free(data->scc_id);
	free(data->depth);
	free(data->fan_in);
	free(data->fan_out);
	free(data->dom_depth);
	free(data->reachable);
	memset(data, 0, sizeof(*data));
}

/*
 * Computes all metrics that require a single shared adjacency build:
 *   fan-in, fan-out, predicate/join counts, density,
 *   reachability + longest-path (BFS), back-edge count (via dominators),
 *   and dominator-tree depth (from the already-computed idom array).
 *
 * Must be called AFTER data->labels and data->idom are populated.
 * All per-node arrays are in vertex-list order (matching data->labels[]).
 */
static void analysis_extra_metrics(const struct nGraph *G, int start_label,
	struct analysis_data *data)
{
	struct graph_adj adj;
	unsigned long *dom = NULL;
	int *bfs_dist = NULL;
	int *queue = NULL;
	int n = data->n;
	int words;
	int start;
	int i;

	if (n <= 0) return;   /* guard for GCC's calloc range analysis across inlining */

	if (build_adjacency(G, 0, &adj) != 0) {
		return;
	}

	start = map_label_to_index(adj.map, n, start_label);

	/* ── Fan-in / Fan-out ─────────────────────────────────────────────────
	 * adj is in label-sorted order; data->labels[] is in vertex-list order.
	 * Convert each vertex-list index i → adj index via map_label_to_index.  */
	for (i = 0; i < n; i++) {
		int ai = map_label_to_index(adj.map, n, data->labels[i]);
		if (ai < 0) continue;
		data->fan_out[i] = adj.out_start[ai + 1] - adj.out_start[ai];
		data->fan_in[i]  = adj.in_start[ai + 1]  - adj.in_start[ai];
		if (data->fan_out[i] >= 2) data->predicate_count++;
		if (data->fan_in[i]  >= 2) data->join_count++;
	}

	/* ── Graph density E / (N*(N-1)) ──────────────────────────────────── */
	if (n > 1) {
		data->density = (double)G->E->count / ((double)n * (double)(n - 1));
	}

	/* ── BFS from entry: reachability + longest path ──────────────────── */
	bfs_dist = (int *)calloc((size_t)n, sizeof(int));
	queue    = (int *)calloc((size_t)n, sizeof(int));
	if (bfs_dist != NULL && queue != NULL && start >= 0) {
		int qh = 0, qt = 0;
		for (i = 0; i < n; i++) bfs_dist[i] = -1;
		bfs_dist[start] = 0;
		queue[qt++] = start;
		while (qh < qt) {
			int v = queue[qh++];
			int j;
			for (j = adj.out_start[v]; j < adj.out_start[v + 1]; j++) {
				int w = adj.out[j];
				if (bfs_dist[w] < 0) {
					bfs_dist[w] = bfs_dist[v] + 1;
					queue[qt++] = w;
				}
			}
		}
		/* Map bfs_dist (adj order) → reachable (vertex-list order) */
		for (i = 0; i < n; i++) {
			int ai = map_label_to_index(adj.map, n, data->labels[i]);
			int d  = (ai >= 0) ? bfs_dist[ai] : -1;
			if (d >= 0) {
				data->reachable[i] = 1;
				if (d > data->longest_path) data->longest_path = d;
			} else {
				data->reachable[i] = 0;
				data->dead_count++;
			}
		}
	}
	free(bfs_dist);
	free(queue);

	/* ── Back edges: (u→v) is a back edge iff v dominates u ──────────── */
	words = dom_word_count(n);
	dom   = (unsigned long *)calloc((size_t)n * (size_t)words, sizeof(unsigned long));
	if (dom != NULL && start >= 0) {
		if (compute_dominators_internal(&adj, start, dom, words, 1) == 0) {
			for (i = 0; i < n; i++) {
				int j;
				for (j = adj.out_start[i]; j < adj.out_start[i + 1]; j++) {
					int v = adj.out[j];
					/* v dominates i → edge (i→v) is a back edge */
					if (dom_test_bit(&dom[i * words], v)) {
						data->back_edges++;
					}
				}
			}
		}
	}
	free(dom);

	/* ── Dominator-tree depth from data->idom[] ───────────────────────
	 * data->idom[i] = label of idom of vertex-list node i.
	 * Fixed-point: dom_depth[root] = 0, then propagate through idom links.
	 * Terminates in at most n passes (tree depth bounded by n-1).         */
	{
		int start_vl = -1; /* start node in vertex-list order */
		int changed;
		for (i = 0; i < n; i++) {
			if (data->labels[i] == start_label) { start_vl = i; break; }
		}
		if (start_vl >= 0) {
			data->dom_depth[start_vl] = 0;
			changed = 1;
			while (changed) {
				changed = 0;
				for (i = 0; i < n; i++) {
					int j;
					if (i == start_vl) continue;
					if (data->dom_depth[i] >= 0) continue;
					if (data->idom[i] < 0) continue;
					/* find parent in vertex-list order */
					for (j = 0; j < n; j++) {
						if (data->labels[j] == data->idom[i]
								&& data->dom_depth[j] >= 0) {
							data->dom_depth[i] = data->dom_depth[j] + 1;
							changed = 1;
							break;
						}
					}
				}
			}
		}
	}

	free_graph_adj(&adj);
}

static int build_analysis_data(struct nGraph *G, int start_label, int exit_label,
	struct analysis_data *data)
{
	int n = G->V->count;
	int i;

	memset(data, 0, sizeof(*data));
	if (n <= 0) {
		return -1;
	}

	data->labels    = (int *)calloc((size_t)n, sizeof(int));
	data->idom      = (int *)calloc((size_t)n, sizeof(int));
	data->ipdom     = (int *)calloc((size_t)n, sizeof(int));
	data->scc_id    = (int *)calloc((size_t)n, sizeof(int));
	data->depth     = (int *)calloc((size_t)n, sizeof(int));
	data->fan_in    = (int *)calloc((size_t)n, sizeof(int));
	data->fan_out   = (int *)calloc((size_t)n, sizeof(int));
	data->dom_depth = (int *)calloc((size_t)n, sizeof(int));
	data->reachable = (int *)calloc((size_t)n, sizeof(int));
	if (!data->labels  || !data->idom   || !data->ipdom || !data->scc_id ||
	    !data->depth   || !data->fan_in || !data->fan_out ||
	    !data->dom_depth || !data->reachable) {
		free_analysis_data(data);
		return -1;
	}

	for (i = 0; i < n; i++) {
		data->idom[i]      = -1;
		data->ipdom[i]     = -1;
		data->scc_id[i]    = -1;
		data->depth[i]     = -1;
		data->dom_depth[i] = -1;  /* -1 = unreachable in dom tree */
	}

	if (getVertexLabels(G, data->labels, n) != n) {
		free_analysis_data(data);
		return -1;
	}

	data->n         = n;
	data->complexity = cyclomaticComplexity(G);
	data->scc_count  = computeSCCs(G, data->scc_id);
	data->max_depth  = computeLoopNestingDepth(G, start_label, data->depth);

	if (computeImmediateDominators(G, start_label, data->idom) != 0) {
		for (i = 0; i < n; i++) data->idom[i] = -1;
	}
	if (computeImmediatePostDominators(G, exit_label, data->ipdom) != 0) {
		for (i = 0; i < n; i++) data->ipdom[i] = -1;
	}

	/* Compute the remaining RE metrics (fan-in/out, density, reachability,
	 * back edges, dom-tree depth) in a single shared adjacency pass.
	 * Must be called after idom is populated (dom_depth depends on it). */
	analysis_extra_metrics(G, start_label, data);

	return 0;
}

int printAnalysisTable(struct nGraph *G, int start_label, int exit_label)
{
	struct analysis_data data;
	int i;
	/* Summary: ~30 chars/metric * ~12 metrics; per-node: ~80 chars * n */
	size_t cap = 1024 + (size_t)G->V->count * 192;
	size_t len = 0;
	char *buf;

	if (G == NULL || G->V == NULL || G->E == NULL) {
		return -1;
	}
	buf = (char *)calloc(cap, 1);
	if (buf == NULL) {
		return -1;
	}
	if (build_analysis_data(G, start_label, exit_label, &data) != 0) {
		free(buf);
		return -1;
	}

	/* ── Graph-level summary ─────────────────────────────────────────── */
	appendf(buf, cap, &len, "+---------------------------+---------------+\n");
	appendf(buf, cap, &len, "| Metric                    | Value         |\n");
	appendf(buf, cap, &len, "+---------------------------+---------------+\n");
	appendf(buf, cap, &len, "| Nodes (N)                 | %13d |\n", data.n);
	appendf(buf, cap, &len, "| Edges (E)                 | %13d |\n", G->E->count);
	appendf(buf, cap, &len, "| Cyclomatic Complexity (M) | %13d |\n", data.complexity);
	appendf(buf, cap, &len, "| Strongly Conn. Components | %13d |\n", data.scc_count);
	appendf(buf, cap, &len, "| Max Loop Nesting Depth    | %13d |\n", data.max_depth);
	appendf(buf, cap, &len, "| Graph Density E/(N(N-1))  | %13.4f |\n", data.density);
	appendf(buf, cap, &len, "| Back Edges (Natural Loops)| %13d |\n", data.back_edges);
	appendf(buf, cap, &len, "| Predicate Nodes (fout>=2) | %13d |\n", data.predicate_count);
	appendf(buf, cap, &len, "| Join Nodes (fin>=2)       | %13d |\n", data.join_count);
	appendf(buf, cap, &len, "| Dead/Unreachable Nodes    | %13d |\n", data.dead_count);
	appendf(buf, cap, &len, "| Longest Path from Entry   | %13d |\n", data.longest_path);
	appendf(buf, cap, &len, "+---------------------------+---------------+\n\n");

	/* ── Per-node table ──────────────────────────────────────────────── */
	appendf(buf, cap, &len,
		"+------+-----+------+----------+----------+-----+------+-----+-------+\n");
	appendf(buf, cap, &len,
		"| Node | Fin | Fout |  IDom    |  IPDom   | LpD | DomD | SCC | Reach |\n");
	appendf(buf, cap, &len,
		"+------+-----+------+----------+----------+-----+------+-----+-------+\n");
	for (i = 0; i < data.n; i++) {
		char idom_c[12], ipdom_c[12], lpd_c[8], domd_c[8], scc_c[8];
		format_cell(idom_c,  sizeof(idom_c),  data.idom[i]);
		format_cell(ipdom_c, sizeof(ipdom_c), data.ipdom[i]);
		format_cell(lpd_c,   sizeof(lpd_c),   data.depth[i]);
		format_cell(domd_c,  sizeof(domd_c),  data.dom_depth[i]);
		format_cell(scc_c,   sizeof(scc_c),   data.scc_id[i]);
		appendf(buf, cap, &len,
			"| %-4d | %-3d | %-4d | %-8s | %-8s | %-3s | %-4s | %-3s | %-5d |\n",
			data.labels[i],
			data.fan_in[i], data.fan_out[i],
			idom_c, ipdom_c,
			lpd_c, domd_c, scc_c,
			data.reachable[i]);
	}
	appendf(buf, cap, &len,
		"+------+-----+------+----------+----------+-----+------+-----+-------+\n");

	printf("%s", buf);
	free(buf);
	free_analysis_data(&data);
	return 0;
}

char *analysisTableDotHtml(struct nGraph *G, int start_label, int exit_label)
{
	struct analysis_data data;
	int i;
	size_t cap;
	size_t len = 0;
	char *buf;

	if (G == NULL || G->V == NULL || G->E == NULL) {
		return NULL;
	}
	/* Each per-node row ~220 chars + ~900 for summary header */
	cap = 900 + (size_t)G->V->count * 220;
	buf = (char *)calloc(cap, 1);
	if (buf == NULL) {
		return NULL;
	}
	if (build_analysis_data(G, start_label, exit_label, &data) != 0) {
		free(buf);
		return NULL;
	}

	/* ── Outer table ──────────────────────────────────────────────────── */
	appendf(buf, cap, &len,
		"<TABLE BORDER=\"0\" CELLBORDER=\"1\" CELLSPACING=\"0\" CELLPADDING=\"3\" BGCOLOR=\"#FAFAFA\">");

	/* ── Title ─────────────────────────────────────────────────────────── */
	appendf(buf, cap, &len,
		"<TR><TD COLSPAN=\"6\" BGCOLOR=\"#2C3E50\" ALIGN=\"CENTER\">"
		"<FONT COLOR=\"white\"><B>CFG / RE Analysis</B></FONT></TD></TR>");

	/* ── Graph-level summary (two columns per row) ─────────────────────── */
	appendf(buf, cap, &len,
		"<TR>"
		"<TD BGCOLOR=\"#3498DB\"><FONT COLOR=\"white\"><B>Nodes (N)</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\"><B>%d</B></TD>"
		"<TD BGCOLOR=\"#3498DB\"><FONT COLOR=\"white\"><B>Edges (E)</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\"><B>%d</B></TD>"
		"<TD BGCOLOR=\"#3498DB\"><FONT COLOR=\"white\"><B>Cyclomatic (M=E-N+2P)</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\"><B>%d</B></TD>"
		"</TR>",
		data.n, G->E->count, data.complexity);

	appendf(buf, cap, &len,
		"<TR>"
		"<TD BGCOLOR=\"#3498DB\"><FONT COLOR=\"white\"><B>SCCs</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\">%d</TD>"
		"<TD BGCOLOR=\"#3498DB\"><FONT COLOR=\"white\"><B>Max Loop Depth</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\">%d</TD>"
		"<TD BGCOLOR=\"#3498DB\"><FONT COLOR=\"white\"><B>Graph Density</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\">%.4f</TD>"
		"</TR>",
		data.scc_count, data.max_depth, data.density);

	appendf(buf, cap, &len,
		"<TR>"
		"<TD BGCOLOR=\"#27AE60\"><FONT COLOR=\"white\"><B>Back Edges</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\">%d</TD>"
		"<TD BGCOLOR=\"#27AE60\"><FONT COLOR=\"white\"><B>Predicate Nodes</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\">%d</TD>"
		"<TD BGCOLOR=\"#27AE60\"><FONT COLOR=\"white\"><B>Join Nodes</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\">%d</TD>"
		"</TR>",
		data.back_edges, data.predicate_count, data.join_count);

	appendf(buf, cap, &len,
		"<TR>"
		"<TD BGCOLOR=\"#E74C3C\"><FONT COLOR=\"white\"><B>Dead Nodes</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\">%d</TD>"
		"<TD BGCOLOR=\"#E74C3C\"><FONT COLOR=\"white\"><B>Longest Path</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\">%d</TD>"
		"<TD COLSPAN=\"2\"></TD>"
		"</TR>",
		data.dead_count, data.longest_path);

	/* ── Per-node header ───────────────────────────────────────────────── */
	appendf(buf, cap, &len,
		"<TR BGCOLOR=\"#2C3E50\">"
		"<TD><FONT COLOR=\"white\"><B>Node</B></FONT></TD>"
		"<TD><FONT COLOR=\"white\"><B>Fin|Fout</B></FONT></TD>"
		"<TD><FONT COLOR=\"white\"><B>IDom</B></FONT></TD>"
		"<TD><FONT COLOR=\"white\"><B>IPDom</B></FONT></TD>"
		"<TD><FONT COLOR=\"white\"><B>LpD|DomD</B></FONT></TD>"
		"<TD><FONT COLOR=\"white\"><B>SCC|Reach</B></FONT></TD>"
		"</TR>");

	/* ── Per-node rows ─────────────────────────────────────────────────── */
	for (i = 0; i < data.n; i++) {
		char idom_c[12], ipdom_c[12], lpd_c[8], domd_c[8], scc_c[8];
		const char *bg = (data.reachable[i] == 0) ? " BGCOLOR=\"#FADBD8\"" : "";
		format_cell(idom_c,  sizeof(idom_c),  data.idom[i]);
		format_cell(ipdom_c, sizeof(ipdom_c), data.ipdom[i]);
		format_cell(lpd_c,   sizeof(lpd_c),   data.depth[i]);
		format_cell(domd_c,  sizeof(domd_c),  data.dom_depth[i]);
		format_cell(scc_c,   sizeof(scc_c),   data.scc_id[i]);
		appendf(buf, cap, &len,
			"<TR%s>"
			"<TD ALIGN=\"CENTER\"><B>%d</B></TD>"
			"<TD ALIGN=\"CENTER\">%d|%d</TD>"
			"<TD ALIGN=\"CENTER\">%s</TD>"
			"<TD ALIGN=\"CENTER\">%s</TD>"
			"<TD ALIGN=\"CENTER\">%s|%s</TD>"
			"<TD ALIGN=\"CENTER\">%s|%d</TD>"
			"</TR>",
			bg,
			data.labels[i],
			data.fan_in[i], data.fan_out[i],
			idom_c, ipdom_c,
			lpd_c, domd_c,
			scc_c, data.reachable[i]);
	}

	appendf(buf, cap, &len, "</TABLE>");

	free_analysis_data(&data);
	return buf;
}
