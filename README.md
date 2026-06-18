# NUCES Graph
## C based Graph Theory Library

**NUCES Graph** is a simple C based library that can be used to analyze small to
medium sized graph networks. The code in the library is continuously improved
and developed by students (as a learning effort) taking the course *CS629
Networks & Graph Theory* at the *National University of Computer and Emerging
Sciences (NUCES), Peshawar, Pakistan*.

[![Documentation](https://img.shields.io/badge/docs-doxygen-blue.svg)](https://omar-khan.github.io/NUCES-Graph/)

### Installation

For compilation of the library, run `./autogen.sh`, `./configure.sh`, followed by `make`. This
will add relevant linking files to the *lib* folder.

For compilation of tests included in the library, run `make tests`, which are
then accessible in the *tests* folder. For compilation of API documentation, run `make docs` 

For a system-wide installation, run `make install` with super user privileges.
The installation defaults to */usr/local/lib* for linking files, and
*/usr/local/include* for the header files. To change the default location, pass 
the new location to configure script as `./configure.sh --prefix=/new/location`.

Lastly, configure the `LD_LIBRARY_PATH` environment variable to point to the
installation location. Otherwise, you have to specify this location every time
you run your code.

### Example Usage

```c
#include <nucesGraph.h>

int main(int argc, char **argv)
{
    int i, def_weight = 1;

    nGraph G = newGraph("G");

    for (i = 0; i < 5; i++) {
        addVertex(&G, i);
    }

    for (i = 0; i < 5; i++) {
        addRandomEdge(&G, def_weight);
    }

    addEdge(&G, 0, 1, def_weight);
}
```

For compiling this code, run `gcc code.c -lNucesGraph`. Specify the include and
linking file locations using `-I` and `-L` if necessary (E.g. `gcc code.c -I/usr/local/include -L/usr/local/lib -lNucesGraph`).

To run the code, run directly as `./a.out`, or if the `LD_LIBRARY_PATH` is not
configured, as `LD_LIBRARY_PATH=/location/of/lib/folder ./a.out`.

### Documentation

The documentation is provided as doxygen generated pdf and html pages. The HTML
pages can be generated using `make doc` while the pdf can be generated using
`make docpdf`.

### Uninstall

The build directory can be cleaned using `make clean`. System wide installation
can be reversed by calling `make uninstall` with super user privileges.

---

### CFG and Reverse Engineering Analysis

The library includes a static analysis module in `src/analysis.c` that treats an
`nGraph` as a control-flow graph and computes common reverse engineering and
software maintenance metrics.

#### Graph-Level Metrics

| Metric | Formula | Description |
|---|---|---|
| **Cyclomatic Complexity** | `M = E - N + 2P` | McCabe complexity. `E` = edges, `N` = nodes, `P` = weakly connected components. |
| **Strongly Connected Components** | Tarjan SCC | Counts SCCs in the graph. |
| **Max Loop Nesting Depth** | Back-edge dominance | Depth of the deepest natural loop. |
| **Graph Density** | `E / (N * (N - 1))` | Directed graph density. |
| **Back Edges** | `(u -> v)` where `v` dominates `u` | Natural loop count proxy. |
| **Predicate Nodes** | `fan_out >= 2` | Decision points in the CFG. |
| **Join Nodes** | `fan_in >= 2` | Control-flow merge points. |
| **Dead / Unreachable Nodes** | BFS from entry | Nodes not reachable from the chosen entry. |
| **Longest Path from Entry** | Max BFS depth | Longest reachable depth from entry. |

#### Per-Node Columns

| Column | Description |
|---|---|
| **Fin** | In-degree |
| **Fout** | Out-degree |
| **IDom** | Immediate dominator label |
| **IPDom** | Immediate post-dominator label |
| **LpD** | Loop nesting depth |
| **DomD** | Dominator-tree depth |
| **SCC** | Strongly connected component identifier |
| **Reach** | `1` if reachable from entry, else `0` |

#### Analysis API

```c
int printAnalysisTable(struct nGraph *G, int start_label, int exit_label);

void showAnalysisPdf(struct nGraph *G, int start_label, int exit_label);
void exportAnalysisDot(struct nGraph *G, int start_label, int exit_label);
void showDotWithAnalysis(struct nGraph *G, int start_label, int exit_label);
void exportDotWithAnalysis(struct nGraph *G, int start_label, int exit_label);

int cfgDfsOrder(struct nGraph *G, int start_label, int *order, int max_order);
int cfgBfsOrder(struct nGraph *G, int start_label, int *order, int max_order);

int dominatorWordCount(int n);
int computeDominators(struct nGraph *G, int start_label,
                      unsigned long *dom, int word_count);
int computePostDominators(struct nGraph *G, int exit_label,
                          unsigned long *pdom, int word_count);
int computeImmediateDominators(const struct nGraph *G, int start_label, int *idom, int idom_len);
int computeImmediatePostDominators(const struct nGraph *G, int exit_label, int *ipdom, int ipdom_len);

int computeSCCs(const struct nGraph *G, int *scc_id, int scc_id_len);
int computeLoopNestingDepth(const struct nGraph *G, int start_label, int *depth_out, int depth_len);

int sliceForward(struct nGraph *dep, int start_label, int *mark, int mark_len);
int sliceBackward(struct nGraph *dep, int start_label, int *mark, int mark_len);
```

#### Example Usage

```c
#include <nucesGraph.h>

int main(void)
{
    nGraph G = newGraph("MyCFG");
    G.directed = 1;

    addVertex(&G, 1);
    addVertex(&G, 2);
    addVertex(&G, 3);
    addVertex(&G, 4);

    addEdgeDirected(&G, 1, 2, 0);
    addEdgeDirected(&G, 1, 3, 0);
    addEdgeDirected(&G, 2, 4, 0);
    addEdgeDirected(&G, 3, 4, 0);

    printAnalysisTable(&G, 1, 4);
    showAnalysisPdf(&G, 1, 4);
    showDotWithAnalysis(&G, 1, 4);

    nGraphFree(&G);
    return 0;
}
```

Compile with: `gcc cfg.c -lNucesGraph`
