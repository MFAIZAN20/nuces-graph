#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
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

struct loop_mark_state {
	const struct graph_adj *adj;
	const unsigned long *dom;
	int capacity;
	int words;
	int header;
	int *loop_stack;
	int *in_loop;
	int top;
};

static int compute_cfg_order_internal(const struct nGraph *G, int start_label,
	int *order, int max_order, int breadth_first);
static int compute_dominator_sets_internal(const struct nGraph *G, int root_label,
	unsigned long *dom, int word_count, int use_preds);
static int compute_immediate_dominators_internal(const struct nGraph *G, int root_label,
	int *idom, int idom_len, int use_preds);
static int compute_slice_internal(const struct nGraph *dep, int start_label, int *mark,
	int mark_len, int use_predecessors);

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

static int edge_vertex_indices(const struct graph_adj *adj,
	const struct edge *edge, int *head, int *tail)
{
	*head = map_label_to_index(adj->map, adj->n, edge->head);
	*tail = map_label_to_index(adj->map, adj->n, edge->tail);
	return (*head >= 0 && *tail >= 0) ? 0 : -1;
}

static int count_adjacency_degrees(const struct nGraph *G, int undirected_mode,
	const struct graph_adj *adj, int *out_degree, int *in_degree)
{
	for (const struct edge *edge = G->E->head; edge != NULL; edge = edge->next) {
		int head;
		int tail;
		if (edge_vertex_indices(adj, edge, &head, &tail) != 0) {
			return -1;
		}
		out_degree[head]++;
		in_degree[tail]++;
		if (undirected_mode || edge->directed == 0) {
			out_degree[tail]++;
			in_degree[head]++;
		}
	}
	return 0;
}

static void build_adjacency_offsets(struct graph_adj *adj,
	const int *out_degree, const int *in_degree)
{
	for (int i = 0; i < adj->n; i++) {
		adj->out_start[i + 1] = adj->out_start[i] + out_degree[i];
		adj->in_start[i + 1] = adj->in_start[i] + in_degree[i];
	}
	adj->m = adj->out_start[adj->n];
}

static int fill_adjacency_edges(const struct nGraph *G, int undirected_mode,
	struct graph_adj *adj, int *out_cursor, int *in_cursor)
{
	for (const struct edge *edge = G->E->head; edge != NULL; edge = edge->next) {
		int head;
		int tail;
		if (edge_vertex_indices(adj, edge, &head, &tail) != 0) {
			return -1;
		}
		adj->out[out_cursor[head]++] = tail;
		adj->in[in_cursor[tail]++] = head;
		if (undirected_mode || edge->directed == 0) {
			adj->out[out_cursor[tail]++] = head;
			adj->in[in_cursor[head]++] = tail;
		}
	}
	return 0;
}

static int build_adjacency(const struct nGraph *G, int undirected_mode, struct graph_adj *adj)
{
	int n = 0;
	int *out_deg = NULL;
	int *in_deg = NULL;
	int *out_cursor = NULL;
	int *in_cursor = NULL;

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

	if (count_adjacency_degrees(G, undirected_mode, adj,
	    out_deg, in_deg) != 0) {
		free(out_deg);
		free(in_deg);
		free_graph_adj(adj);
		return -1;
	}

	adj->out_start = (int *)calloc((size_t)(n + 1), sizeof(int));
	adj->in_start = (int *)calloc((size_t)(n + 1), sizeof(int));
	if (adj->out_start == NULL || adj->in_start == NULL) {
		free(out_deg);
		free(in_deg);
		free_graph_adj(adj);
		return -1;
	}

	build_adjacency_offsets(adj, out_deg, in_deg);
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
	if (G->E->head != NULL && (adj->out == NULL || adj->in == NULL)) {
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
	for (int i = 0; i < n; i++) {
		out_cursor[i] = adj->out_start[i];
		in_cursor[i] = adj->in_start[i];
	}

	if (fill_adjacency_edges(G, undirected_mode, adj,
	    out_cursor, in_cursor) != 0) {
		free(out_deg);
		free(in_deg);
		free(out_cursor);
		free(in_cursor);
		free_graph_adj(adj);
		return -1;
	}

	free(out_deg);
	free(in_deg);
	free(out_cursor);
	free(in_cursor);
	return 0;
}

int getVertexLabels(const struct nGraph *G, int *labels, int max_labels)
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

int graphNodeCount(const struct nGraph *G)
{
	if (G == NULL || G->V == NULL) return -1;
	return G->V->count;
}

int graphEdgeCount(const struct nGraph *G)
{
	if (G == NULL || G->E == NULL) return -1;
	return G->E->count;
}

static int enqueue_component_successor(const struct graph_adj *adj,
	int successor, int *visited, int *queue, int *tail)
{
	if (successor < 0 || successor >= adj->n) {
		return -1;
	}
	if (visited[successor]) {
		return 0;
	}
	if (*tail >= adj->n) {
		return -1;
	}
	visited[successor] = 1;
	queue[*tail] = successor;
	(*tail)++;
	return 0;
}

static int visit_weak_component(const struct graph_adj *adj, int root,
	int *visited, int *queue)
{
	int head = 0;
	int tail = 0;

	if (root < 0 || root >= adj->n) {
		return -1;
	}
	visited[root] = 1;
	queue[tail++] = root;
	while (head < tail) {
		if (head >= adj->n) {
			return -1;
		}
		int vertex = queue[head++];
		for (int i = adj->out_start[vertex];
		     i < adj->out_start[vertex + 1]; i++) {
			int successor = adj->out[i];
			if (enqueue_component_successor(adj, successor, visited,
			    queue, &tail) != 0) {
				return -1;
			}
		}
	}
	return 0;
}

static int count_weak_components(const struct graph_adj *adj)
{
	int n = adj->n;
	int *visited = (int *)calloc((size_t)n, sizeof(int));
	int *queue = (int *)calloc((size_t)n, sizeof(int));
	int components = 0;

	if (visited == NULL || queue == NULL) {
		free(visited);
		free(queue);
		return -1;
	}

	for (int i = 0; i < n; i++) {
		if (visited[i]) continue;
		components++;
		if (visit_weak_component(adj, i, visited, queue) != 0) {
			components = -1;
			break;
		}
	}

	free(visited);
	free(queue);
	return components;
}

int cyclomaticComplexity(const struct nGraph *G)
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

int cfgDfsOrder(const struct nGraph *G, int start_label, int *order, int max_order)
{
	return compute_cfg_order_internal(G, start_label, order, max_order, 0);
}

int cfgBfsOrder(const struct nGraph *G, int start_label, int *order, int max_order)
{
	return compute_cfg_order_internal(G, start_label, order, max_order, 1);
}

static int push_cfg_successors(const struct graph_adj *adj, int vertex,
	int breadth_first, int *visited, int *work, int *tail)
{
	int begin = adj->out_start[vertex];
	int end = adj->out_start[vertex + 1];
	int index = breadth_first ? begin : end - 1;
	int step = breadth_first ? 1 : -1;

	while ((breadth_first && index < end) ||
	       (!breadth_first && index >= begin)) {
		int successor = adj->out[index];
		if (!visited[successor]) {
			if (*tail >= adj->n) {
				return -1;
			}
			visited[successor] = 1;
			work[*tail] = successor;
			(*tail)++;
		}
		index += step;
	}
	return 0;
}

static int compute_cfg_order_internal(const struct nGraph *G, int start_label,
	int *order, int max_order, int breadth_first)
{
	struct graph_adj adj;
	int *visited = NULL;
	int *work = NULL;
	int count = 0;
	int start = -1;
	int head = 0;
	int tail = 0;

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
	work = (int *)calloc((size_t)adj.n, sizeof(int));
	if (visited == NULL || work == NULL) {
		free(visited);
		free(work);
		free_graph_adj(&adj);
		return -1;
	}

	visited[start] = 1;
	work[tail++] = start;
	while (head < tail) {
		int v;

		if (breadth_first) {
			v = work[head++];
		} else {
			v = work[--tail];
			head = 0;
		}
		if (count < max_order) {
			order[count++] = adj.labels[v];
		}
		if (push_cfg_successors(&adj, v, breadth_first,
		    visited, work, &tail) != 0) {
			count = -1;
			break;
		}
	}

	free(visited);
	free(work);
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
	for (int i = 0; i < words; i++) {
		set[i] = ~0UL;
	}
}

static void dom_set_zero(unsigned long *set, int words)
{
	for (int i = 0; i < words; i++) {
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
	for (int i = 0; i < words; i++) {
		dst[i] &= src[i];
	}
}

static int dom_set_equal(const unsigned long *a, const unsigned long *b, int words)
{
	for (int i = 0; i < words; i++) {
		if (a[i] != b[i]) return 0;
	}
	return 1;
}

int dominatorWordCount(int n)
{
	return dom_word_count(n);
}

static int intersect_predecessor_sets(const struct graph_adj *adj, int vertex,
	const unsigned long *dom, int words, int use_preds, unsigned long *result)
{
	const int *offsets = use_preds ? adj->in_start : adj->out_start;
	const int *edges = use_preds ? adj->in : adj->out;
	int begin = offsets[vertex];
	int end = offsets[vertex + 1];

	dom_set_all(result, words);
	if (begin == end) {
		dom_set_zero(result, words);
		return 0;
	}
	for (int i = begin; i < end; i++) {
		dom_set_and(result, &dom[edges[i] * words], words);
	}
	return 1;
}

static int compute_dominators_internal(const struct graph_adj *adj, int start,
	unsigned long *dom, int words, int use_preds)
{
	int n = adj->n;
	int changed = 1;
	unsigned long *newset = NULL;

	if (dom == NULL || n <= 0 || start < 0 || start >= n ||
	    words < dom_word_count(n)) {
		return -1;
	}
	newset = (unsigned long *)calloc((size_t)words, sizeof(unsigned long));
	if (newset == NULL) {
		return -1;
	}

	for (int v = 0; v < n; v++) {
		dom_set_all(&dom[v * words], words);
	}
	dom_set_zero(&dom[start * words], words);
	dom_set_bit(&dom[start * words], start);

	while (changed) {
		changed = 0;
		for (int v = 0; v < n; v++) {
			unsigned long *cur = &dom[v * words];
			if (v == start) {
				continue;
			}

			intersect_predecessor_sets(adj, v, dom, words, use_preds, newset);
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

int computeDominators(const struct nGraph *G, int start_label, unsigned long *dom, int word_count)
{
	return compute_dominator_sets_internal(G, start_label, dom, word_count, 1);
}

int computePostDominators(const struct nGraph *G, int exit_label, unsigned long *pdom, int word_count)
{
	return compute_dominator_sets_internal(G, exit_label, pdom, word_count, 0);
}

static int compute_dominator_sets_internal(const struct nGraph *G, int root_label,
	unsigned long *dom, int word_count, int use_preds)
{
	struct graph_adj adj;
	int root = -1;
	int rc = 0;

	if (G == NULL || G->V == NULL || G->E == NULL || dom == NULL) {
		return -1;
	}
	if (build_adjacency(G, 0, &adj) != 0) {
		return -1;
	}
	root = map_label_to_index(adj.map, adj.n, root_label);
	if (root < 0) {
		free_graph_adj(&adj);
		return -1;
	}
	if (word_count < dom_word_count(adj.n)) {
		free_graph_adj(&adj);
		return -1;
	}

	rc = compute_dominators_internal(&adj, root, dom, word_count, use_preds);
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
	for (int i = 0; i < words; i++) {
		unsigned long v = set[i];
		while (v) {
			count += (int)(v & 1UL);
			v >>= 1;
		}
	}
	return count;
}

int computeImmediateDominators(const struct nGraph *G, int start_label, int *idom,
	int idom_len)
{
	return compute_immediate_dominators_internal(G, start_label, idom, idom_len, 1);
}

int computeImmediatePostDominators(const struct nGraph *G, int exit_label, int *ipdom,
	int ipdom_len)
{
	return compute_immediate_dominators_internal(G, exit_label, ipdom, ipdom_len, 0);
}

static int compute_immediate_dominators_internal(const struct nGraph *G, int root_label,
	int *idom, int idom_len, int use_preds)
{
	struct graph_adj adj;
	unsigned long *dom = NULL;
	int words = 0;
	int root = -1;
	int n = 0;

	if (G == NULL || G->V == NULL || G->E == NULL || idom == NULL) {
		return -1;
	}
	if (build_adjacency(G, 0, &adj) != 0) {
		return -1;
	}
	root = map_label_to_index(adj.map, adj.n, root_label);
	if (root < 0) {
		free_graph_adj(&adj);
		return -1;
	}

	n = adj.n;
	if (idom_len < n) {
		free_graph_adj(&adj);
		return -1;
	}
	words = dom_word_count(n);
	dom = (unsigned long *)calloc((size_t)n * (size_t)words, sizeof(unsigned long));
	if (dom == NULL) {
		free_graph_adj(&adj);
		return -1;
	}

	if (compute_dominators_internal(&adj, root, dom, words, use_preds) != 0) {
		free(dom);
		free_graph_adj(&adj);
		return -1;
	}

	/*
	 * idom(v) is the unique dominator d of v where |dom(d)| = |dom(v)| - 1.
	 * Using popcount on the bitset reduces the triple-loop O(n^3) to O(n^2*words).
	 */
	for (int v = 0; v < n; v++) {
		int v_count;
		idom[v] = -1;
		if (v == root) {
			idom[v] = adj.labels[root];
			continue;
		}
		v_count = dom_popcount(&dom[v * words], words);
		if (v_count <= 1) {
			continue;
		}
		for (int d = 0; d < n; d++) {
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

struct scc_state {
	int *index;
	int *low;
	int *stack;
	int *onstack;
	int *scc_id;
	struct scc_frame *frames;
	int stack_size;
	int next_index;
	int scc_count;
	int capacity;
};

static int scc_push_vertex(const struct graph_adj *adj, struct scc_state *state,
	int vertex, int *frame_count)
{
	if (vertex < 0 || vertex >= state->capacity ||
	    state->stack_size >= state->capacity ||
	    *frame_count >= state->capacity) {
		return -1;
	}
	state->index[vertex] = state->next_index;
	state->low[vertex] = state->next_index;
	state->next_index++;
	state->stack[state->stack_size++] = vertex;
	state->onstack[vertex] = 1;
	state->frames[*frame_count].v = vertex;
	state->frames[*frame_count].ei = adj->out_start[vertex];
	(*frame_count)++;
	return 0;
}

static int scc_finish_component(struct scc_state *state, int root)
{
	int vertex;

	do {
		if (state->stack_size <= 0) {
			return -1;
		}
		state->stack_size--;
		vertex = state->stack[state->stack_size];
		state->onstack[vertex] = 0;
		state->scc_id[vertex] = state->scc_count;
	} while (vertex != root);
	state->scc_count++;
	return 0;
}

static int scc_process_successor(const struct graph_adj *adj,
	struct scc_state *state, int vertex, int successor, int *frame_count)
{
	if (successor < 0 || successor >= state->capacity) {
		return -1;
	}
	if (state->index[successor] < 0) {
		return scc_push_vertex(adj, state, successor, frame_count);
	}
	if (state->onstack[successor] &&
	    state->index[successor] < state->low[vertex]) {
		state->low[vertex] = state->index[successor];
	}
	return 0;
}

static int scc_finish_frame(struct scc_state *state, int vertex,
	int frame_count)
{
	if (frame_count > 0) {
		int parent = state->frames[frame_count - 1].v;
		if (state->low[vertex] < state->low[parent]) {
			state->low[parent] = state->low[vertex];
		}
	}
	if (state->low[vertex] == state->index[vertex]) {
		return scc_finish_component(state, vertex);
	}
	return 0;
}

static int scc_iterative(const struct graph_adj *adj, int root,
	struct scc_state *state)
{
	int frame_count = 0;

	if (scc_push_vertex(adj, state, root, &frame_count) != 0) {
		return -1;
	}

	while (frame_count > 0) {
		struct scc_frame *f = &state->frames[frame_count - 1];
		int v = f->v;

		if (f->ei < adj->out_start[v + 1]) {
			int w = adj->out[f->ei++];
			if (scc_process_successor(adj, state, v, w, &frame_count) != 0) {
				return -1;
			}
		} else {
			frame_count--;
			if (scc_finish_frame(state, v, frame_count) != 0) {
				return -1;
			}
		}
	}
	return 0;
}

int computeSCCs(const struct nGraph *G, int *scc_id, int scc_id_len)
{
	struct graph_adj adj;
	int *index_arr = NULL;
	int *low = NULL;
	int *stack = NULL;
	int *onstack = NULL;
	struct scc_frame *frames = NULL;
	struct scc_state state;

	if (G == NULL || G->V == NULL || G->E == NULL || scc_id == NULL) {
		return -1;
	}
	if (build_adjacency(G, 0, &adj) != 0) {
		return -1;
	}
	if (scc_id_len < adj.n) {
		free_graph_adj(&adj);
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

	for (int v = 0; v < adj.n; v++) {
		index_arr[v] = -1;
		low[v]       = -1;
		scc_id[v]    = -1;
	}

	state.index = index_arr;
	state.low = low;
	state.stack = stack;
	state.onstack = onstack;
	state.scc_id = scc_id;
	state.frames = frames;
	state.stack_size = 0;
	state.next_index = 0;
	state.scc_count = 0;
	state.capacity = adj.n;

	for (int v = 0; v < adj.n; v++) {
		if (index_arr[v] < 0 &&
		    scc_iterative(&adj, v, &state) != 0) {
			state.scc_count = -1;
			break;
		}
	}

	free(index_arr);
	free(low);
	free(stack);
	free(onstack);
	free(frames);
	free_graph_adj(&adj);
	return state.scc_count;
}

static int add_loop_predecessor(struct loop_mark_state *state, int predecessor)
{
	if (predecessor < 0 || predecessor >= state->capacity) {
		return -1;
	}
	if (state->in_loop[predecessor] ||
	    !dominates(state->dom, state->words, state->header, predecessor)) {
		return 0;
	}
	if (state->top < 0 || state->top >= state->capacity) {
		return -1;
	}
	state->in_loop[predecessor] = 1;
	state->loop_stack[state->top] = predecessor;
	state->top++;
	return 0;
}

static int mark_natural_loop(const struct graph_adj *adj,
	const unsigned long *dom, int words, int tail, int header,
	int *loop_stack, int *in_loop)
{
	struct loop_mark_state state;

	if (tail < 0 || tail >= adj->n || header < 0 || header >= adj->n) {
		return -1;
	}

	memset(in_loop, 0, (size_t)adj->n * sizeof(int));
	state.adj = adj;
	state.dom = dom;
	state.capacity = adj->n;
	state.words = words;
	state.header = header;
	state.loop_stack = loop_stack;
	state.in_loop = in_loop;
	state.top = 0;
	in_loop[header] = 1;
	in_loop[tail] = 1;
	if (state.top >= state.capacity) {
		return -1;
	}
	loop_stack[state.top++] = tail;
	while (state.top > 0) {
		int vertex = loop_stack[--state.top];
		for (int i = adj->in_start[vertex]; i < adj->in_start[vertex + 1]; i++) {
			int predecessor = adj->in[i];
			if (add_loop_predecessor(&state, predecessor) != 0) {
				return -1;
			}
		}
	}
	return 0;
}

struct back_edge_state {
	int *depth_out;
	int depth_len;
	int *max_depth;
	int *loop_stack;
	int *in_loop;
};

static void add_loop_depths(int n, const int *in_loop, struct back_edge_state *state)
{
	int limit;
	if (n <= 0 || in_loop == NULL || state == NULL || state->depth_out == NULL || state->max_depth == NULL || state->depth_len <= 0) {
		return;
	}
	limit = n < state->depth_len ? n : state->depth_len;
	for (int i = 0; i < limit; i++) {
		if (!in_loop[i]) {
			continue;
		}
		state->depth_out[i]++;
		if (state->depth_out[i] > *state->max_depth) {
			*state->max_depth = state->depth_out[i];
		}
	}
}

static int process_back_edges(const struct graph_adj *adj,
	const unsigned long *dom, int words, struct back_edge_state *state)
{
	if (adj->out == NULL) {
		return 0;
	}
	for (int tail = 0; tail < adj->n; tail++) {
		for (int i = adj->out_start[tail]; i < adj->out_start[tail + 1]; i++) {
			int header = adj->out[i];
			if (!dominates(dom, words, header, tail)) {
				continue;
			}
			if (mark_natural_loop(adj, dom, words, tail, header,
			    state->loop_stack, state->in_loop) != 0) {
				return -1;
			}
			add_loop_depths(adj->n, state->in_loop, state);
		}
	}
	return 0;
}

int computeLoopNestingDepth(const struct nGraph *G, int start_label, int *depth_out,
	int depth_len)
{
	struct graph_adj adj;
	int words = 0;
	unsigned long *dom = NULL;
	int *loop_stack = NULL;
	int *in_loop = NULL;
	int start = -1;
	int n = 0;
	int max_depth = 0;
	struct back_edge_state state;

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
	if (depth_len < n) {
		free_graph_adj(&adj);
		return -1;
	}
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

	for (int u = 0; u < n; u++) {
		depth_out[u] = 0;
	}

	loop_stack = (int *)calloc((size_t)n, sizeof(int));
	in_loop = (int *)calloc((size_t)n, sizeof(int));
	if (loop_stack == NULL || in_loop == NULL) {
		free(loop_stack);
		free(in_loop);
		free(dom);
		free_graph_adj(&adj);
		return -1;
	}
	state.depth_out = depth_out;
	state.depth_len = depth_len;
	state.max_depth = &max_depth;
	state.loop_stack = loop_stack;
	state.in_loop = in_loop;

	if (process_back_edges(&adj, dom, words, &state) != 0) {
		max_depth = -1;
	}

	free(loop_stack);
	free(in_loop);
	free(dom);
	free_graph_adj(&adj);
	return max_depth;
}

int sliceForward(const struct nGraph *dep, int start_label, int *mark, int mark_len)
{
	return compute_slice_internal(dep, start_label, mark, mark_len, 0);
}

int sliceBackward(const struct nGraph *dep, int start_label, int *mark, int mark_len)
{
	return compute_slice_internal(dep, start_label, mark, mark_len, 1);
}

static int compute_slice_internal(const struct nGraph *dep, int start_label, int *mark,
	int mark_len, int use_predecessors)
{
	struct graph_adj adj;
	int *stack = NULL;
	int count = 0;
	int start = -1;
	const int *start_offsets = NULL;
	const int *edges = NULL;
	int top = 0;

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

	start_offsets = use_predecessors ? adj.in_start : adj.out_start;
	edges = use_predecessors ? adj.in : adj.out;

	/* Mark on push — each node pushed at most once, stack bounded by n */
	memset(mark, 0, (size_t)adj.n * sizeof(int));
	mark[start] = 1;
	count = 1;
	stack[top++] = start;
	while (top > 0) {
		int v = stack[--top];
		for (int i = start_offsets[v]; i < start_offsets[v + 1]; i++) {
			int w = edges[i];
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

struct text_buffer {
	char *data;
	size_t capacity;
	size_t length;
};

static int buffer_append_text(struct text_buffer *buffer, const char *text)
{
	size_t text_length = strlen(text);

	if (text_length >= buffer->capacity - buffer->length) {
		return -1;
	}
	memcpy(buffer->data + buffer->length, text, text_length + 1U);
	buffer->length += text_length;
	return 0;
}

static int buffer_append_int(struct text_buffer *buffer, int value)
{
	char text[32];
	int length = snprintf(text, sizeof(text), "%d", value);

	if (length < 0 || (size_t)length >= sizeof(text)) {
		return -1;
	}
	return buffer_append_text(buffer, text);
}

static int buffer_append_double(struct text_buffer *buffer, double value)
{
	char text[64];
	int length = snprintf(text, sizeof(text), "%.4f", value);

	if (length < 0 || (size_t)length >= sizeof(text)) {
		return -1;
	}
	return buffer_append_text(buffer, text);
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

static void fill_minus_one_ints(int n, int *values)
{
	if (n <= 0 || values == NULL) {
		return;
	}
	memset(values, 0xFF, (size_t)n * sizeof(*values));
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
static void compute_degree_metrics(const struct graph_adj *adj,
	struct analysis_data *data)
{
	int n = data->n;

	for (int i = 0; i < n; i++) {
		int ai = map_label_to_index(adj->map, n, data->labels[i]);
		if (ai < 0) continue;
		data->fan_out[i] = adj->out_start[ai + 1] - adj->out_start[ai];
		data->fan_in[i] = adj->in_start[ai + 1] - adj->in_start[ai];
		if (data->fan_out[i] >= 2) data->predicate_count++;
		if (data->fan_in[i] >= 2) data->join_count++;
	}
}

static void map_reachability(const struct graph_adj *adj, const int *distance,
	struct analysis_data *data)
{
	for (int i = 0; i < data->n; i++) {
		int ai = map_label_to_index(adj->map, data->n, data->labels[i]);
		int node_distance = ai >= 0 ? distance[ai] : -1;

		if (node_distance < 0) {
			data->dead_count++;
			continue;
		}
		data->reachable[i] = 1;
		if (node_distance > data->longest_path) {
			data->longest_path = node_distance;
		}
	}
}

static int enqueue_successors(const struct graph_adj *adj, int vertex,
	int *distance, int *queue, int *tail)
{
	for (int i = adj->out_start[vertex]; i < adj->out_start[vertex + 1]; i++) {
		int successor = adj->out[i];
		if (successor < 0 || successor >= adj->n) {
			return -1;
		}
		if (distance[successor] >= 0) {
			continue;
		}
		if (*tail >= adj->n) {
			return -1;
		}
		distance[successor] = distance[vertex] + 1;
		queue[*tail] = successor;
		(*tail)++;
	}
	return 0;
}

static int dequeue_vertex(const struct graph_adj *adj, const int *queue,
	int *head, int tail)
{
	if (*head < 0 || *head >= tail || *head >= adj->n) {
		return -1;
	}
	return queue[(*head)++];
}

static int traverse_reachability(const struct graph_adj *adj, int start,
	int *distance, int *queue)
{
	int head = 0;
	int tail = 0;

	distance[start] = 0;
	queue[tail++] = start;
	while (head < tail) {
		int vertex = dequeue_vertex(adj, queue, &head, tail);
		if (vertex < 0) {
			return -1;
		}
		if (enqueue_successors(adj, vertex, distance, queue, &tail) != 0) {
			return -1;
		}
	}
	return 0;
}

static void compute_reachability_metrics(const struct graph_adj *adj, int start,
	struct analysis_data *data)
{
	int *distance = NULL;
	int *queue = NULL;

	if (start < 0) {
		data->dead_count = data->n;
		return;
	}
	if (data->n != adj->n) {
		return;
	}
	distance = (int *)malloc((size_t)adj->n * sizeof(int));
	queue = (int *)malloc((size_t)adj->n * sizeof(int));
	if (distance == NULL || queue == NULL) {
		free(distance);
		free(queue);
		return;
	}
	for (int i = 0; i < adj->n; i++) {
		distance[i] = -1;
	}
	if (adj->out != NULL &&
	    traverse_reachability(adj, start, distance, queue) != 0) {
		map_reachability(adj, distance, data);
		free(distance);
		free(queue);
		return;
	}
	map_reachability(adj, distance, data);
	free(distance);
	free(queue);
}

static int count_vertex_back_edges(const struct graph_adj *adj,
	const unsigned long *dom, int words, int tail)
{
	int count = 0;

	for (int i = adj->out_start[tail]; i < adj->out_start[tail + 1]; i++) {
		int header = adj->out[i];
		if (header >= 0 && header < adj->n &&
		    dom_test_bit(&dom[tail * words], header)) {
			count++;
		}
	}
	return count;
}

static void compute_back_edge_metric(const struct graph_adj *adj, int start,
	struct analysis_data *data)
{
	int words = dom_word_count(data->n);
	unsigned long *dom = NULL;

	if (start < 0 || adj->out == NULL) {
		return;
	}
	dom = (unsigned long *)calloc((size_t)data->n * (size_t)words,
		sizeof(unsigned long));
	if (dom == NULL) {
		return;
	}
	if (compute_dominators_internal(adj, start, dom, words, 1) == 0) {
		for (int tail = 0; tail < data->n; tail++) {
			data->back_edges += count_vertex_back_edges(adj, dom, words, tail);
		}
	}
	free(dom);
}

static int find_label_index(const int *labels, int n, int label)
{
	for (int i = 0; i < n; i++) {
		if (labels[i] == label) {
			return i;
		}
	}
	return -1;
}

static int update_dominator_depth(struct analysis_data *data, int root)
{
	int changed = 0;

	for (int i = 0; i < data->n; i++) {
		int parent;
		if (i == root || data->dom_depth[i] >= 0 || data->idom[i] < 0) {
			continue;
		}
		parent = find_label_index(data->labels, data->n, data->idom[i]);
		if (parent >= 0 && data->dom_depth[parent] >= 0) {
			data->dom_depth[i] = data->dom_depth[parent] + 1;
			changed = 1;
		}
	}
	return changed;
}

static void compute_dominator_depths(int start_label, struct analysis_data *data)
{
	int root = find_label_index(data->labels, data->n, start_label);

	if (root < 0) {
		return;
	}
	data->dom_depth[root] = 0;
	while (update_dominator_depth(data, root)) {
		/* Continue until every reachable parent depth has propagated. */
	}
}

static void analysis_extra_metrics(const struct nGraph *G, int start_label,
	struct analysis_data *data)
{
	struct graph_adj adj;
	int start;

	if (data->n <= 0 || build_adjacency(G, 0, &adj) != 0) {
		return;
	}
	start = map_label_to_index(adj.map, adj.n, start_label);
	compute_degree_metrics(&adj, data);
	if (data->n > 1) {
		data->density = (double)G->E->count /
			((double)data->n * (double)(data->n - 1));
	}
	compute_reachability_metrics(&adj, start, data);
	compute_back_edge_metric(&adj, start, data);
	compute_dominator_depths(start_label, data);
	free_graph_adj(&adj);
}

static int build_analysis_data(const struct nGraph *G, int start_label,
	int exit_label,
	struct analysis_data *data)
{
	int n = G->V->count;

	memset(data, 0, sizeof(*data));
	if (n <= 0) {
		return -1;
	}
	data->n = n;

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

	fill_minus_one_ints(n, data->idom);
	fill_minus_one_ints(n, data->ipdom);
	fill_minus_one_ints(n, data->scc_id);
	fill_minus_one_ints(n, data->depth);
	fill_minus_one_ints(n, data->dom_depth);  /* -1 = unreachable in dom tree */

	if (getVertexLabels(G, data->labels, n) != n) {
		free_analysis_data(data);
		return -1;
	}

	data->complexity = cyclomaticComplexity(G);
	data->scc_count  = computeSCCs(G, data->scc_id, n);
	data->max_depth  = computeLoopNestingDepth(G, start_label, data->depth, n);

	if (computeImmediateDominators(G, start_label, data->idom, n) != 0) {
		for (int i = 0; i < n; i++) data->idom[i] = -1;
	}
	if (computeImmediatePostDominators(G, exit_label, data->ipdom, n) != 0) {
		for (int i = 0; i < n; i++) data->ipdom[i] = -1;
	}

	/* Compute the remaining RE metrics (fan-in/out, density, reachability,
	 * back edges, dom-tree depth) in a single shared adjacency pass.
	 * Must be called after idom is populated (dom_depth depends on it). */
	analysis_extra_metrics(G, start_label, data);

	return 0;
}

int printAnalysisTable(const struct nGraph *G, int start_label, int exit_label)
{
	struct analysis_data data;

	if (G == NULL || G->V == NULL || G->E == NULL) {
		return -1;
	}
	if (build_analysis_data(G, start_label, exit_label, &data) != 0) {
		return -1;
	}

	printf("+---------------------------+---------------+\n");
	printf("| Metric                    | Value         |\n");
	printf("+---------------------------+---------------+\n");
	printf("| Nodes (N)                 | %13d |\n", data.n);
	printf("| Edges (E)                 | %13d |\n", G->E->count);
	printf("| Cyclomatic Complexity (M) | %13d |\n", data.complexity);
	printf("| Strongly Conn. Components | %13d |\n", data.scc_count);
	printf("| Max Loop Nesting Depth    | %13d |\n", data.max_depth);
	printf("| Graph Density E/(N(N-1))  | %13.4f |\n", data.density);
	printf("| Back Edges (Natural Loops)| %13d |\n", data.back_edges);
	printf("| Predicate Nodes (fout>=2) | %13d |\n", data.predicate_count);
	printf("| Join Nodes (fin>=2)       | %13d |\n", data.join_count);
	printf("| Dead/Unreachable Nodes    | %13d |\n", data.dead_count);
	printf("| Longest Path from Entry   | %13d |\n", data.longest_path);
	printf("+---------------------------+---------------+\n\n");
	printf("+------+-----+------+----------+----------+-----+------+-----+-------+\n");
	printf("| Node | Fin | Fout |  IDom    |  IPDom   | LpD | DomD | SCC | Reach |\n");
	printf("+------+-----+------+----------+----------+-----+------+-----+-------+\n");
	for (int i = 0; i < data.n; i++) {
		char idom_c[12];
		char ipdom_c[12];
		char lpd_c[8];
		char domd_c[8];
		char scc_c[8];
		format_cell(idom_c,  sizeof(idom_c),  data.idom[i]);
		format_cell(ipdom_c, sizeof(ipdom_c), data.ipdom[i]);
		format_cell(lpd_c,   sizeof(lpd_c),   data.depth[i]);
		format_cell(domd_c,  sizeof(domd_c),  data.dom_depth[i]);
		format_cell(scc_c,   sizeof(scc_c),   data.scc_id[i]);
		printf(
			"| %-4d | %-3d | %-4d | %-8s | %-8s | %-3s | %-4s | %-3s | %-5d |\n",
			data.labels[i],
			data.fan_in[i], data.fan_out[i],
			idom_c, ipdom_c,
			lpd_c, domd_c, scc_c,
			data.reachable[i]);
	}
	printf("+------+-----+------+----------+----------+-----+------+-----+-------+\n");
	free_analysis_data(&data);
	return 0;
}

char *analysisTableDotHtml(const struct nGraph *G, int start_label, int exit_label)
{
	struct analysis_data data;
	size_t cap;
	char *buf;
	struct text_buffer output;

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
	output.data = buf;
	output.capacity = cap;
	output.length = 0;

	buffer_append_text(&output,
		"<TABLE BORDER=\"0\" CELLBORDER=\"1\" CELLSPACING=\"0\" CELLPADDING=\"3\" BGCOLOR=\"#FAFAFA\">");
	buffer_append_text(&output,
		"<TR><TD COLSPAN=\"6\" BGCOLOR=\"#2C3E50\" ALIGN=\"CENTER\">"
		"<FONT COLOR=\"white\"><B>CFG / RE Analysis</B></FONT></TD></TR>");
	buffer_append_text(&output,
		"<TR>"
		"<TD BGCOLOR=\"#3498DB\"><FONT COLOR=\"white\"><B>Nodes (N)</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\"><B>");
	buffer_append_int(&output, data.n);
	buffer_append_text(&output, "</B></TD>"
		"<TD BGCOLOR=\"#3498DB\"><FONT COLOR=\"white\"><B>Edges (E)</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\"><B>");
	buffer_append_int(&output, G->E->count);
	buffer_append_text(&output, "</B></TD>"
		"<TD BGCOLOR=\"#3498DB\"><FONT COLOR=\"white\"><B>Cyclomatic (M=E-N+2P)</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\"><B>");
	buffer_append_int(&output, data.complexity);
	buffer_append_text(&output, "</B></TD></TR>");
	buffer_append_text(&output,
		"<TR>"
		"<TD BGCOLOR=\"#3498DB\"><FONT COLOR=\"white\"><B>SCCs</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\">");
	buffer_append_int(&output, data.scc_count);
	buffer_append_text(&output, "</TD>"
		"<TD BGCOLOR=\"#3498DB\"><FONT COLOR=\"white\"><B>Max Loop Depth</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\">");
	buffer_append_int(&output, data.max_depth);
	buffer_append_text(&output, "</TD>"
		"<TD BGCOLOR=\"#3498DB\"><FONT COLOR=\"white\"><B>Graph Density</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\">");
	buffer_append_double(&output, data.density);
	buffer_append_text(&output, "</TD></TR>");
	buffer_append_text(&output,
		"<TR>"
		"<TD BGCOLOR=\"#27AE60\"><FONT COLOR=\"white\"><B>Back Edges</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\">");
	buffer_append_int(&output, data.back_edges);
	buffer_append_text(&output, "</TD>"
		"<TD BGCOLOR=\"#27AE60\"><FONT COLOR=\"white\"><B>Predicate Nodes</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\">");
	buffer_append_int(&output, data.predicate_count);
	buffer_append_text(&output, "</TD>"
		"<TD BGCOLOR=\"#27AE60\"><FONT COLOR=\"white\"><B>Join Nodes</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\">");
	buffer_append_int(&output, data.join_count);
	buffer_append_text(&output, "</TD></TR>");
	buffer_append_text(&output,
		"<TR>"
		"<TD BGCOLOR=\"#E74C3C\"><FONT COLOR=\"white\"><B>Dead Nodes</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\">");
	buffer_append_int(&output, data.dead_count);
	buffer_append_text(&output, "</TD>"
		"<TD BGCOLOR=\"#E74C3C\"><FONT COLOR=\"white\"><B>Longest Path</B></FONT></TD>"
		"<TD ALIGN=\"RIGHT\">");
	buffer_append_int(&output, data.longest_path);
	buffer_append_text(&output, "</TD><TD COLSPAN=\"2\"></TD></TR>");
	buffer_append_text(&output,
		"<TR BGCOLOR=\"#2C3E50\">"
		"<TD><FONT COLOR=\"white\"><B>Node</B></FONT></TD>"
		"<TD><FONT COLOR=\"white\"><B>Fin|Fout</B></FONT></TD>"
		"<TD><FONT COLOR=\"white\"><B>IDom</B></FONT></TD>"
		"<TD><FONT COLOR=\"white\"><B>IPDom</B></FONT></TD>"
		"<TD><FONT COLOR=\"white\"><B>LpD|DomD</B></FONT></TD>"
		"<TD><FONT COLOR=\"white\"><B>SCC|Reach</B></FONT></TD>"
		"</TR>");

	for (int i = 0; i < data.n; i++) {
		char idom_c[12];
		char ipdom_c[12];
		char lpd_c[8];
		char domd_c[8];
		char scc_c[8];
		const char *bg = (data.reachable[i] == 0) ? " BGCOLOR=\"#FADBD8\"" : "";
		format_cell(idom_c,  sizeof(idom_c),  data.idom[i]);
		format_cell(ipdom_c, sizeof(ipdom_c), data.ipdom[i]);
		format_cell(lpd_c,   sizeof(lpd_c),   data.depth[i]);
		format_cell(domd_c,  sizeof(domd_c),  data.dom_depth[i]);
		format_cell(scc_c,   sizeof(scc_c),   data.scc_id[i]);
		buffer_append_text(&output, "<TR");
		buffer_append_text(&output, bg);
		buffer_append_text(&output, "><TD ALIGN=\"CENTER\"><B>");
		buffer_append_int(&output, data.labels[i]);
		buffer_append_text(&output, "</B></TD><TD ALIGN=\"CENTER\">");
		buffer_append_int(&output, data.fan_in[i]);
		buffer_append_text(&output, "|");
		buffer_append_int(&output, data.fan_out[i]);
		buffer_append_text(&output, "</TD><TD ALIGN=\"CENTER\">");
		buffer_append_text(&output, idom_c);
		buffer_append_text(&output, "</TD><TD ALIGN=\"CENTER\">");
		buffer_append_text(&output, ipdom_c);
		buffer_append_text(&output, "</TD><TD ALIGN=\"CENTER\">");
		buffer_append_text(&output, lpd_c);
		buffer_append_text(&output, "|");
		buffer_append_text(&output, domd_c);
		buffer_append_text(&output, "</TD><TD ALIGN=\"CENTER\">");
		buffer_append_text(&output, scc_c);
		buffer_append_text(&output, "|");
		buffer_append_int(&output, data.reachable[i]);
		buffer_append_text(&output, "</TD></TR>");
	}

	buffer_append_text(&output, "</TABLE>");

	free_analysis_data(&data);
	return buf;
}
