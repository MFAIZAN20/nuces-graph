#include "nucesGraph.h"
#include <stdio.h>
#include <stdlib.h>

static int index_of(const int *labels, int n, int label)
{
	for (int i = 0; i < n; i++) {
		if (labels[i] == label) return i;
	}
	return -1;
}

static int dom_test(const unsigned long *dom, int words, int a_idx, int b_idx)
{
	int bits = (int)(sizeof(unsigned long) * 8U);
	int w = a_idx / bits;
	int b = a_idx % bits;
	return (dom[b_idx * words + w] >> (unsigned int)b) & 1UL;
}

static int check(int cond, const char *msg)
{
	if (!cond) {
		fprintf(stderr, "analysis test failed: %s\n", msg);
		return 0;
	}
	return 1;
}

int main(void)
{
	int ok = 1;
	struct nGraph G = newGraph("CFG");
	addVertex(&G, 10);
	addVertex(&G, 20);
	addVertex(&G, 30);
	addVertex(&G, 40);

	addEdgeDirected(&G, 10, 20, 1);
	addEdgeDirected(&G, 20, 30, 1);
	addEdgeDirected(&G, 10, 40, 1);
	addEdgeDirected(&G, 40, 30, 1);

	ok &= check(graphNodeCount(&G) == 4, "node count");
	ok &= check(graphEdgeCount(&G) == 4, "edge count");
	ok &= check(cyclomaticComplexity(&G) == 2, "cyclomatic complexity");

	int labels[4];
	ok &= check(getVertexLabels(&G, labels, 4) == 4, "getVertexLabels count");
	ok &= check(labels[0] == 10 && labels[1] == 20 && labels[2] == 30 && labels[3] == 40, "label order");

	int dfs[4];
	int bfs[4];
	ok &= check(cfgDfsOrder(&G, 10, dfs, 4) == 4, "dfs count");
	ok &= check(cfgBfsOrder(&G, 10, bfs, 4) == 4, "bfs count");
	ok &= check(dfs[0] == 10 && dfs[1] == 20 && dfs[2] == 30 && dfs[3] == 40, "dfs order");
	ok &= check(bfs[0] == 10 && bfs[1] == 20 && bfs[2] == 40 && bfs[3] == 30, "bfs order");

	int words = dominatorWordCount(4);
	unsigned long *dom = (unsigned long *)calloc((size_t)4 * (size_t)words, sizeof(unsigned long));
	unsigned long *pdom = (unsigned long *)calloc((size_t)4 * (size_t)words, sizeof(unsigned long));
	int *idom = (int *)calloc(4, sizeof(int));
	int *ipdom = (int *)calloc(4, sizeof(int));
	if (!dom || !pdom || !idom || !ipdom) {
		fprintf(stderr, "analysis test failed: alloc\n");
		free(dom);
		free(pdom);
		free(idom);
		free(ipdom);
		nGraphFree(&G);
		return 1;
	}

	ok &= check(computeDominators(&G, 10, dom, words) == 0, "compute dominators");
	ok &= check(computePostDominators(&G, 30, pdom, words) == 0, "compute post dominators");
	ok &= check(computeImmediateDominators(&G, 10, idom, 4) == 0, "compute idom");
	ok &= check(computeImmediatePostDominators(&G, 30, ipdom, 4) == 0, "compute ipdom");

	int a = index_of(labels, 4, 10);
	int b = index_of(labels, 4, 20);
	int c = index_of(labels, 4, 30);
	int d = index_of(labels, 4, 40);
	ok &= check(a >= 0 && b >= 0 && c >= 0 && d >= 0, "label index");

	ok &= check(dom_test(dom, words, a, b), "A dominates B");
	ok &= check(dom_test(dom, words, a, c), "A dominates C");
	ok &= check(dom_test(dom, words, a, d), "A dominates D");
	ok &= check(!dom_test(dom, words, b, c), "B does not dominate C");
	ok &= check(!dom_test(dom, words, d, c), "D does not dominate C");

	ok &= check(dom_test(pdom, words, c, a), "C post-dominates A");
	ok &= check(dom_test(pdom, words, c, b), "C post-dominates B");
	ok &= check(dom_test(pdom, words, c, d), "C post-dominates D");

	ok &= check(idom[a] == 10, "idom A");
	ok &= check(idom[b] == 10, "idom B");
	ok &= check(idom[c] == 10, "idom C");
	ok &= check(idom[d] == 10, "idom D");

	ok &= check(ipdom[a] == 30, "ipdom A");
	ok &= check(ipdom[b] == 30, "ipdom B");
	ok &= check(ipdom[c] == 30, "ipdom C");
	ok &= check(ipdom[d] == 30, "ipdom D");

	int depth[4];
	int max_depth = computeLoopNestingDepth(&G, 10, depth, 4);
	ok &= check(max_depth == 0, "loop depth max");
	ok &= check(depth[0] == 0 && depth[1] == 0 && depth[2] == 0 && depth[3] == 0, "loop depth per node");

	int mark[4];
	ok &= check(sliceForward(&G, 10, mark, 4) == 4, "slice forward from A");
	ok &= check(sliceForward(&G, 20, mark, 4) == 2, "slice forward from B");
	ok &= check(sliceBackward(&G, 30, mark, 4) == 4, "slice backward from C");

	ok &= check(printAnalysisTable(&G, 10, 30) == 0, "print analysis table");

	free(dom);
	free(pdom);
	free(idom);
	free(ipdom);
	nGraphFree(&G);

	struct nGraph H = newGraph("SCC");
	addVertex(&H, 1);
	addVertex(&H, 2);
	addVertex(&H, 3);
	addEdgeDirected(&H, 1, 2, 1);
	addEdgeDirected(&H, 2, 3, 1);
	addEdgeDirected(&H, 3, 1, 1);

	int scc_id[3];
	int scc_count = computeSCCs(&H, scc_id, 3);
	ok &= check(scc_count == 1, "scc count");
	ok &= check(scc_id[0] == scc_id[1] && scc_id[1] == scc_id[2], "scc ids");
	nGraphFree(&H);

	struct nGraph U = newGraph("Unreachable");
	addVertex(&U, 1);
	addVertex(&U, 2);
	addVertex(&U, 3);
	addEdgeDirected(&U, 1, 2, 1);

	int unreachable_idom[3];
	int unreachable_depth[3];
	ok &= check(computeImmediateDominators(&U, 1, unreachable_idom, 3) == 0,
		"unreachable idom computation");
	ok &= check(unreachable_idom[0] == 1 && unreachable_idom[1] == 1 &&
		unreachable_idom[2] == -1, "unreachable idom values");
	ok &= check(computeLoopNestingDepth(&U, 1, unreachable_depth, 3) == 0,
		"unreachable loop depth");
	nGraphFree(&U);

	struct nGraph Z = newGraph("NoEdges");
	addVertex(&Z, 7);
	char *empty_analysis = analysisTableDotHtml(&Z, 7, 7);
	ok &= check(empty_analysis != NULL, "empty-edge analysis");
	free(empty_analysis);
	nGraphFree(&Z);

	if (!ok) return 1;
	printf("analysis tests passed\n");
	return 0;
}
