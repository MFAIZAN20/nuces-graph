#include "nucesGraph.h"
#include <ctype.h>
#include <fcntl.h>
#include <glpk.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

static void write_graphviz_vertices(FILE *tf, const struct nGraph *B,
                                    int include_colors);
static void write_graphviz_edges(FILE *tf, const struct nGraph *B,
                                 int preserve_undirected);
static void write_analysis_table_dot(FILE *tf, const struct nGraph *G,
                                     const char *table, int compact_graph);

static char *sanitize_filename_component(const char *label) {
  size_t len = 0;
  char *sanitized = NULL;

  if (label != NULL) {
    len = strlen(label);
  }
  if (len == 0) {
    sanitized = (char *)malloc(sizeof("graph"));
    if (sanitized != NULL) {
      strcpy(sanitized, "graph");
    }
    return sanitized;
  }

  sanitized = (char *)malloc(len + 1U);
  if (sanitized == NULL) {
    return NULL;
  }

  for (size_t i = 0; i < len; i++) {
    unsigned char ch = (unsigned char)label[i];
    sanitized[i] =
        (isalnum(ch) || ch == '-' || ch == '_' || ch == '.') ? (char)ch : '_';
  }
  sanitized[len] = '\0';
  return sanitized;
}

static char *build_output_filename(const char *label, const char *suffix) {
  char *safe_label = sanitize_filename_component(label);
  char *filename = NULL;
  size_t size = 0;

  if (safe_label == NULL) {
    return NULL;
  }

  size = strlen(safe_label) + strlen(suffix) + 1U;
  filename = (char *)malloc(size);
  if (filename != NULL) {
    snprintf(filename, size, "%s%s", safe_label, suffix);
  }
  free(safe_label);
  return filename;
}

static int create_temp_file_path(const char *context, char **path_out) {
  const char *directory = P_tmpdir;
  char *path = NULL;
  size_t size;
  int fd;

  size = strlen(directory) + sizeof("/nucesGraphXXXXXX");
  path = (char *)malloc(size);
  if (path == NULL) {
    return -1;
  }
  snprintf(path, size, "%s/nucesGraphXXXXXX", directory);
  fd = mkstemp(path);
  if (fd == -1) {
    fprintf(stderr, "%s: Temp File Creation Error\n", context);
    free(path);
    return -1;
  }
  close(fd);
  *path_out = path;
  return 0;
}

/** Exports a Graph in DIMAC format to a file. Output is written to G.dimac
 * where G is the graph label. At the moment, only network flow graphs are
 * supported.
 * @param G Graph object
 * @param file Filename
 */
void exportDimac(struct nGraph *G) {
  char *filename = NULL;
  char *s = NULL;
  FILE *fd = NULL;
  if (G == NULL || G->E == NULL) {
    return;
  }

  s = (char *)calloc((size_t)(G->E->count * 16 + 64), sizeof(char));
  filename = build_output_filename(G->label, ".dimac");
  if (s == NULL || filename == NULL) {
    free(s);
    free(filename);
    return;
  }

  makeDimac(G, s);
  fd = fopen(filename, "w");
  if (fd == NULL) {
    fprintf(stderr, "Error: Could not open %s for writing\n", filename);
    free(s);
    free(filename);
    return;
  }

  fputs(s, fd);
  fclose(fd);
  fprintf(stdout, "Output written to %s\n", filename);
  free(s);
  free(filename);
}

void makeDimac(struct nGraph *G, char *contents) {
  char returnString[128];
  struct edge *tmp = NULL;

  if (G == NULL || G->V == NULL || G->E == NULL || contents == NULL) {
    return;
  }
  contents[0] = '\0';
  tmp = G->E->head;

  snprintf(returnString, sizeof(returnString), "p max %d %d\n", G->V->count,
           G->E->count);
  strcat(contents, returnString);

  snprintf(returnString, sizeof(returnString), "n %d s\n", 1);
  strcat(contents, returnString);

  snprintf(returnString, sizeof(returnString), "n %d t\n", G->V->count);
  strcat(contents, returnString);

  while (tmp != NULL && G->E->count > 0) {
    if (tmp->directed == 1) {
      snprintf(returnString, sizeof(returnString), "a %d %d %d\n", tmp->head,
               tmp->tail, tmp->weight);
    } else {
      snprintf(returnString, sizeof(returnString), "p %d %d\n", tmp->head,
               tmp->tail);
    }
    strcat(contents, returnString);
    tmp = tmp->next;
  }
}

/** Exports GLPK Compatible Linear Program for solving maxflow problems
 * @param G Graph object
 */
void exportGLPK(struct nGraph *G) {
  // Conversion is not direct. First Dimac output is created. Then that Dimac
  // output is again plugged into the GLPK API to create GLPK compatible LP code

  char *code = NULL;
  char *filename = NULL;
  char *filenamelp = NULL;
  FILE *fd = NULL;
  glp_prob *lp = NULL;
  glp_graph *g = NULL;
  int s;
  int t;

  if (G == NULL || G->E == NULL) {
    return;
  }

  code = (char *)calloc((size_t)(G->E->count * 32 + 256), sizeof(char));
  if (code == NULL ||
      create_temp_file_path("exportGLPK()", &filename) != 0 ||
      create_temp_file_path("exportGLPK()", &filenamelp) != 0) {
    free(code);
    free(filename);
    free(filenamelp);
    return;
  }

  makeDimac(G, code);
  fd = fopen(filename, "w");
  if (fd != NULL) {
    fputs(code, fd);
    fclose(fd);
  } else {
    fprintf(stderr, "Error: Could not open %s for writing\n", filename);
    free(code);
    free(filename);
    free(filenamelp);
    return;
  }
  fprintf(stdout, "Dimac Output written to %s\n", filename);

  g = glp_create_graph(0, sizeof(double));
  glp_read_maxflow(g, &s, &t, 0, filename);
  lp = glp_create_prob();
  glp_maxflow_lp(lp, g, GLP_ON, s, t, 0);
  glp_delete_graph(g);
  glp_write_lp(lp, NULL, filenamelp);
  glp_delete_prob(lp);

  // unlink(filename);  // for deleting dot file
  free(code);
  free(filename);
  free(filenamelp);
}

/** Shows a Graph in DIMAC format. Can be manually copied pasted from the
 * terminal.
 * @param G Graph object
 */
void showDimac(struct nGraph *G) {
  char *s = calloc(1024, sizeof(char));
  if (s == NULL) {
    return;
  }
  makeDimac(G, s);

  printf("%s\n", s);
  free(s);
}

void exportDot(struct nGraph *B) {
  /**
   * Export Dot is this function
   */
  fprintf(stdout, "graph %s {\n", B->label);
  //	printf("\toverlap = false;\n");
  fprintf(stdout, "\tnode [shape=\"circle\"];\n");

  struct vertex *tmpV = B->V->head;
  while (tmpV != NULL) {
    fprintf(stdout, "\t%d", tmpV->label);
    if (tmpV->lblString != NULL) {
      fprintf(stdout, " [label=\"%s\"]", tmpV->lblString);
    }
    fprintf(stdout, ";\n");
    tmpV = tmpV->next;
  }

  struct edge *tmp = B->E->head;
  while (tmp != NULL) {
    if (tmp->weight == 0) {
      fprintf(stdout, "\t%d -- %d;\n", tmp->head, tmp->tail);
    } else {
      fprintf(stdout, "\t%d -- %d[label=\"%d\"];\n", tmp->head, tmp->tail,
              tmp->weight);
    }
    tmp = tmp->next;
  }
  fprintf(stdout, "}\n");
}

void exportDotWithAnalysis(struct nGraph *B, int start_label, int exit_label) {
  char *table = analysisTableDotHtml(B, start_label, exit_label);

  fprintf(stdout, "digraph %s {\n", B->label);
  fprintf(stdout, "\tgraph [labelloc=\"b\", labeljust=\"l\"];\n");
  if (table != NULL) {
    fprintf(stdout, "\tlabel=<%s>;\n", table);
  }
  fprintf(stdout, "\tnode [shape=\"circle\"];\n");
  write_graphviz_vertices(stdout, B, 0);
  write_graphviz_edges(stdout, B, 0);
  fprintf(stdout, "}\n");
  free(table);
}

static int create_render_paths(const char *context, char **filename_out,
                               char **fileimage_out) {
  char *filename = NULL;
  char *fileimage = NULL;

  if (create_temp_file_path(context, &filename) != 0) {
    return -1;
  }
  if (create_temp_file_path(context, &fileimage) != 0) {
    unlink(filename);
    free(filename);
    return -1;
  }
  *filename_out = filename;
  *fileimage_out = fileimage;
  return 0;
}

static void write_graphviz_vertices(FILE *tf, const struct nGraph *B,
                                    int include_colors) {
  struct vertex *tmpV = B->V->head;

  while (tmpV != NULL) {
    fprintf(tf, "\t%d", tmpV->label);
    if (tmpV->lblString != NULL || (include_colors && tmpV->color != -1)) {
      fprintf(tf, " [");
      if (tmpV->lblString != NULL) {
        fprintf(tf, "label=\"%s\"", tmpV->lblString);
      }
      if (tmpV->lblString != NULL && include_colors && tmpV->color != -1) {
        fprintf(tf, ", ");
      }
      if (include_colors && tmpV->color != -1) {
        fprintf(tf, "style=filled, color=\"%s\"", colors[tmpV->color]);
      }
      fprintf(tf, "]");
    }
    fprintf(tf, ";\n");
    tmpV = tmpV->next;
  }
}

static void write_graphviz_edges(FILE *tf, const struct nGraph *B,
                                 int preserve_undirected) {
  struct edge *tmp = B->E->head;

  while (tmp != NULL) {
    if (preserve_undirected && tmp->directed == 0) {
      if (tmp->weight == 0) {
        fprintf(tf, "\t%d -> %d[dir=none];\n", tmp->head, tmp->tail);
      } else {
        fprintf(tf, "\t%d -> %d[dir=none,label=\"%d\"];\n", tmp->head,
                tmp->tail, tmp->weight);
      }
    } else {
      if (tmp->weight == 0) {
        fprintf(tf, "\t%d -> %d;\n", tmp->head, tmp->tail);
      } else {
        fprintf(tf, "\t%d -> %d[label=\"%d\"];\n", tmp->head, tmp->tail,
                tmp->weight);
      }
    }
    tmp = tmp->next;
  }
}

static void choose_display_layout(struct nGraph *B, const char **overlap,
                                  const char **splines) {
  int V = B->V->count;
  int E = B->E->count;
  double density = V > 0 ? (double)E / V : 0.0;
  int maxDegree = 0;
  struct vertex *tmpV = B->V->head;

  *overlap = "false";
  *splines = "curved";

  while (tmpV != NULL) {
    int deg = tmpV->degree_in + tmpV->degree_out;
    if (deg > maxDegree) {
      maxDegree = deg;
    }
    tmpV = tmpV->next;
  }

  if (B->displayType != 0) {
    return;
  }
  if (V > 500) {
    setDisplayType(B, "sfdp");
    *overlap = "true";
    *splines = "false";
  } else if ((maxDegree >= 20 || (V > 0 && maxDegree >= V / 2)) && V > 5) {
    setDisplayType(B, "circular");
    *splines = "line";
  } else if (B->directed && density < 2.0) {
    setDisplayType(B, "dot");
    *splines = "polyline";
  } else if (V > 40 || density >= 2.0) {
    setDisplayType(B, "sfdp");
    *overlap = "scale";
    *splines = "line";
  } else if (B->directed) {
    setDisplayType(B, "dot");
  } else {
    setDisplayType(B, "neato");
  }
}

static int resolve_display_type(int display_type, int fallback_display_type) {
  return display_type == 0 ? fallback_display_type : display_type;
}

static const char *graphviz_binary(int display_type) {
  switch (display_type) {
  case 1:
    return "/usr/bin/sfdp";
  case 2:
    return "/usr/bin/twopi";
  case 3:
    return "/usr/bin/dot";
  default:
    return "/usr/bin/neato";
  }
}

static const char *graphviz_name(int display_type) {
  switch (display_type) {
  case 1:
    return "sfdp";
  case 2:
    return "twopi";
  case 3:
    return "dot";
  default:
    return "neato";
  }
}

static void render_graphviz_pdf(const char *filename, const char *fileimage,
                                int display_type, int fallback_display_type,
                                int print_render_type) {
  int effective_type = resolve_display_type(display_type, fallback_display_type);
  const char *binary = graphviz_binary(effective_type);
  pid_t pid = fork();

  if (pid == 0) {
    if (print_render_type) {
      printf("Source (dot) at: %s\nImage  (pdf) at: %s\nRender Display Type: %s\n\n",
             filename, fileimage, graphviz_name(effective_type));
    } else {
      printf("Source (dot) at: %s\nImage  (pdf) at: %s\n\n", filename,
             fileimage);
    }
    execl(binary, binary, filename, "-T", "pdf", "-o", fileimage, NULL);
    _exit(1);
  } else if (pid > 0) {
    wait(NULL);
  }
}

static void open_rendered_pdf(const char *fileimage) {
  pid_t pid;

  if (!getenv("DISPLAY")) {
    return;
  }
  pid = fork();
  if (pid == 0) {
    execl("/usr/bin/okular", "/usr/bin/okular", fileimage, NULL);
    _exit(1);
  }
}

static void write_analysis_table_dot(FILE *tf, const struct nGraph *G,
                                     const char *table, int compact_graph) {
  fprintf(tf, "digraph \"Analysis_%s\" {\n", G->label);
  fprintf(tf, "\tbgcolor=\"white\";\n");
  if (compact_graph) {
    fprintf(tf, "\tgraph [pad=\"0.5\", nodesep=\"0.5\", ranksep=\"0.5\"];\n");
  } else {
    fprintf(tf, "\tgraph [pad=\"0.5\"];\n");
  }
  fprintf(tf, "\tnode [shape=none, margin=0];\n");
  if (table != NULL) {
    fprintf(tf, "\tanalysis [label=<%s>];\n", table);
  } else {
    fprintf(tf, "\tanalysis [label=\"Analysis unavailable\"];\n");
  }
  fprintf(tf, "}\n");
}

/**
 * Renders a graph $B$ using Graphviz and automatically opens it in a PDF
 * viewer. If displayType of $B$ is set to 0 (which is default), an intelligent
 * auto-configuration engine evaluates the graph's topology and selects the
 * optimal layout algorithm and parameters. This is summarized as:
 * - Massive graphs ($V > 500$): Use 'sfdp', overlaps allowed, straight lines
 * (for performance).
 * - Star/Hub graphs (Max Degree $\geq$ 20 or $\geq V/2$): Use 'twopi'
 * (circular) layout.
 * - Sparse directed graphs (Density $<$ 2.0): Use 'dot' for hierarchical flow.
 * - Dense/Medium graphs (V $>$ 40 or Density $\geq$ 2.0): Uses 'sfdp' with
 * scaled overlap removal.
 * - Small graphs: Uses 'dot' (directed) or 'neato' (undirected) with curved
 * splines.
 *
 * @param B Graph object
 */

void showDot(struct nGraph *B) {
  char *filename = NULL;
  char *fileimage = NULL;
  const char *overlap = "false";
  const char *splines = "curved";
  FILE *tf = NULL;

  if (create_render_paths("showDot()", &filename, &fileimage) != 0) {
    return;
  }

  tf = fopen(filename, "w+");
  if (tf) {
    choose_display_layout(B, &overlap, &splines);
    fprintf(tf, "digraph %s {\n", B->label);
    fprintf(tf, "\tlabel=\"%s\";\n", B->label);
    fprintf(tf, "\tlabelloc=\"t\";\n");
    fprintf(tf, "\tfontsize=20;\n");
    fprintf(tf, "\toverlap = %s;\n", overlap);
    fprintf(tf, "\tsplines = \"%s\";\n", splines);
    fprintf(tf, "\tsep = 3;\n");
    fprintf(tf, "\tnode [shape=\"circle\"];\n");
    write_graphviz_vertices(tf, B, 1);
    write_graphviz_edges(tf, B, 1);
    fprintf(tf, "}\n");
    fclose(tf);
  }

  render_graphviz_pdf(filename, fileimage, B->displayType, 0, 1);
  open_rendered_pdf(fileimage);
  // unlink(filename);  // for deleting dot file
  // unlink(fileimage); // for deleting pdf file
  free(filename);
  free(fileimage);
}

void showDotWithAnalysis(struct nGraph *B, int start_label, int exit_label) {
  char *filename = NULL;
  char *fileimage = NULL;
  char *table = analysisTableDotHtml(B, start_label, exit_label);
  FILE *tf = NULL;

  if (create_render_paths("showDotWithAnalysis()", &filename, &fileimage) !=
      0) {
    free(table);
    return;
  }

  tf = fopen(filename, "w+");
  if (tf) {
    fprintf(tf, "digraph %s {\n", B->label);
    fprintf(tf, "\tgraph [labelloc=\"b\", labeljust=\"l\"];\n");
    if (table != NULL) {
      fprintf(tf, "\tlabel=<%s>;\n", table);
    }
    fprintf(tf, "\toverlap = false;\n");
    fprintf(tf, "\tsplines = \"curved\";\n");
    fprintf(tf, "\tsep = 3;\n");
    fprintf(tf, "\tnode [shape=\"circle\"];\n");
    write_graphviz_vertices(tf, B, 1);
    write_graphviz_edges(tf, B, 1);
    fprintf(tf, "}\n");
    fclose(tf);
  }
  free(table);

  render_graphviz_pdf(filename, fileimage, B->displayType, 0, 0);
  open_rendered_pdf(fileimage);

  free(filename);
  free(fileimage);
}

void show(struct nGraph *B) {
  listVertices(B);
  listEdges(B);
}

void listVerticesAlphabet(struct nGraph *G) {
  struct vertex *tmp = G->V->head;
  printf("%s.V = { ", G->label);
  int count = 0;
  while (tmp != NULL) {
    if (tmp->lblString != NULL) {
      printf("%s", tmp->lblString);
    } else {
      printf("%c", 'a' + tmp->label);
    }
    printf("%s", G->V->count - 1 == count ? " " : ", ");
    tmp = tmp->next;
    count++;
  }
  printf("}\n");
}

void listBK_temp(struct nGraph *G) {
  struct vertex *tmp = G->V->head;
  printf("%s%s = { %s", KRED, G->label, KBLU);
  int count = 0;
  while (tmp != NULL && G->V->count > 0) {
    if (tmp->lblString != NULL && strlen(tmp->lblString) > 0) {
      printf("%s", tmp->lblString);
    } else {
      printf("%d", tmp->label);
    }
    printf("%s", G->V->count - 1 == count ? " " : ", ");
    tmp = tmp->next;
    count++;
  }
  printf("%s}%s", KRED, KWHT);
}

void listVertices(struct nGraph *G) {
  struct vertex *tmp = G->V->head;
  printf("%s%s.V(%d) = { ", KRED, G->label, G->V->count);
  int count = 0;
  while (tmp != NULL && G->V->count > 0) {
    if (tmp->lblString != NULL && strlen(tmp->lblString) > 0) {
      printf("%s", tmp->lblString);
    } else {
      printf("%d", tmp->label);
    }
    printf("%s", G->V->count - 1 == count ? " " : ", ");
    tmp = tmp->next;
    count++;
  }
  printf("}\n%s", KWHT);
}

/**
 * Displays a List of Edges in a Graph. Output shows the Graph Label, followed
 * by a list having format as: (Edge Id) vertex from, vertex to.
 * @param Graph Object
 */

void listEdges(struct nGraph *G) {
  struct edge *tmp = G->E->head;
  printf("%s%s.E(%d) = { ", KRED, G->label, G->E->count);
  int count = 0;
  while (tmp != NULL && G->E->count > 0) {
    if (count == (G->E->count - 1)) {
      printf("(%d) %d->%d ", tmp->label, tmp->head, tmp->tail);
    } else {
      printf("(%d) %d->%d, ", tmp->label, tmp->head, tmp->tail);
    }
    tmp = tmp->next;
    count++;
  }
  printf("}\n%s", KWHT);
}

/**
 * Sets display type of a graph (graphviz)
 * @param Graph Object
 * @param type ["circular", "sfdp", "dot", "neato"]
 */

void setDisplayType(struct nGraph *G, char *type) {
  if (strcmp(type, "circular") == 0) {
    G->displayType = 2;
  } else if (strcmp(type, "sfdp") == 0) {
    G->displayType = 1;
  } else if (strcmp(type, "dot") == 0) {
    G->displayType = 3;
  } else {
    G->displayType = 0;
  }
}

void exportAnalysisDot(const struct nGraph *G, int start_label, int exit_label) {
  char *table = analysisTableDotHtml(G, start_label, exit_label);

  write_analysis_table_dot(stdout, G, table, 0);
  free(table);
}

void showAnalysisPdf(const struct nGraph *G, int start_label, int exit_label) {
  char *filename = NULL;
  char *fileimage = NULL;
  char *table = analysisTableDotHtml(G, start_label, exit_label);
  FILE *tf = NULL;

  if (create_render_paths("showAnalysisPdf()", &filename, &fileimage) != 0) {
    free(table);
    return;
  }

  tf = fopen(filename, "w+");
  if (tf) {
    write_analysis_table_dot(tf, G, table, 1);
    fclose(tf);
  }
  free(table);

  render_graphviz_pdf(filename, fileimage, G->displayType, 3, 0);
  open_rendered_pdf(fileimage);

  free(filename);
  free(fileimage);
}

/**
 * Exports a Graph to TikZ format for use in LaTeX papers.
 * Outputs to a file named after the graph's label (e.g., G.tex).
 * Requires \\usetikzlibrary{graphs,graphdrawing} and LuaLaTeX to render.
 * @param Graph Object
 */
void exportTikZ(struct nGraph *G) {
  char *filename = NULL;
  FILE *fd = NULL;

  if (G == NULL || G->V == NULL || G->E == NULL) {
    return;
  }
  filename = build_output_filename(G->label, ".tex");
  if (filename == NULL) {
    return;
  }
  fd = fopen(filename, "w");
  if (!fd) {
    fprintf(stderr, "Error: Could not open %s for writing\n", filename);
    free(filename);
    return;
  }

  fprintf(fd, "%% Build using LuaLaTeX\n");
  fprintf(fd, "\\documentclass[]{article}\n");
  fprintf(fd, "\\usepackage{tikz}\n");
  fprintf(fd, "\\usetikzlibrary{graphs, graphdrawing}\n");
  fprintf(fd, "\\usegdlibrary{force, layered, trees}\n");
  fprintf(fd, "\\begin{document}\n");

  fprintf(fd, "\\begin{tikzpicture}[>=stealth, every node/.style={circle, "
              "draw, minimum size=0.5cm}]\n");
  fprintf(fd, "\\graph [spring layout] {\n");

  // Process Nodes First
  struct vertex *tmpV = G->V->head;
  while (tmpV != NULL) {
    fprintf(fd, "  %d", tmpV->label);
    if (tmpV->lblString != NULL && strlen(tmpV->lblString) > 0) {
      fprintf(fd, " [as=\"%s\"]", tmpV->lblString);
    }
    fprintf(fd, ";\n");
    tmpV = tmpV->next;
  }
  fprintf(fd, "\n");

  // Edges Now
  struct edge *tmpE = G->E->head;
  while (tmpE != NULL) {
    if (tmpE->directed) {
      fprintf(fd, "  %d -> %d", tmpE->head, tmpE->tail);
    } else {
      fprintf(fd, "  %d -- %d", tmpE->head, tmpE->tail);
    }

    if (tmpE->weight != 0) {
      fprintf(fd, " [edge label=%d]", tmpE->weight);
    }
    fprintf(fd, ";\n");

    tmpE = tmpE->next;
  }

  fprintf(fd, "};\n");
  fprintf(fd, "\\end{tikzpicture}\n");
  fprintf(fd, "\\end{document}\n");
  fclose(fd);
  fprintf(stdout, "TikZ output written to %s\\n", filename);
  free(filename);
}

/**
 * Exports a Graph to GraphML format for loading into Gephi or Cytoscape.
 * Outputs to a file named after the graph's label (e.g., G.graphml).
 * @param Graph Object
 */
void exportGraphML(struct nGraph *G) {
  char *filename = NULL;
  FILE *fd = NULL;

  if (G == NULL || G->V == NULL || G->E == NULL) {
    return;
  }
  filename = build_output_filename(G->label, ".graphml");
  if (filename == NULL) {
    return;
  }
  fd = fopen(filename, "w");
  if (!fd) {
    fprintf(stderr, "Error: Could not open %s for writing\\n", filename);
    free(filename);
    return;
  }

  fprintf(fd, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf(fd, "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\">\n");
  fprintf(fd, "  <key id=\"d0\" for=\"node\" attr.name=\"label\" "
              "attr.type=\"string\" />\n");
  fprintf(fd, "  <key id=\"d1\" for=\"edge\" attr.name=\"weight\" "
              "attr.type=\"int\" />\n");

  // Graph element
  if (G->directed) {
    fprintf(fd, "  <graph id=\"%s\" edgedefault=\"directed\">\n", G->label);
  } else {
    fprintf(fd, "  <graph id=\"%s\" edgedefault=\"undirected\">\n", G->label);
  }

  // Nodes
  struct vertex *tmpV = G->V->head;
  while (tmpV != NULL) {
    fprintf(fd, "    <node id=\"n%d\">\n", tmpV->label);

    char *labelStr = (tmpV->lblString != NULL) ? tmpV->lblString : "";
    if (strlen(labelStr) > 0) {
      fprintf(fd, "      <data key=\"d0\">%s</data>\n", labelStr);
    } else {
      fprintf(fd, "      <data key=\"d0\">%d</data>\n", tmpV->label);
    }

    fprintf(fd, "    </node>\n");
    tmpV = tmpV->next;
  }

  // Edges
  struct edge *tmpE = G->E->head;
  int edge_id = 0;
  while (tmpE != NULL) {
    fprintf(fd, "    <edge id=\"e%d\" source=\"n%d\" target=\"n%d\"", edge_id++,
            tmpE->head, tmpE->tail);

    // If specific edge directedness differs from graph default
    if (tmpE->directed != G->directed) {
      fprintf(fd, " directed=\"%s\"", tmpE->directed ? "true" : "false");
    }

    if (tmpE->weight != 0) {
      fprintf(fd, ">\n      <data key=\"d1\">%d</data>\n    </edge>\n",
              tmpE->weight);
    } else {
      fprintf(fd, " />\n");
    }

    tmpE = tmpE->next;
  }

  fprintf(fd, "  </graph>\n");
  fprintf(fd, "</graphml>\n");

  fclose(fd);
  fprintf(stdout, "GraphML output written to %s\n", filename);
  free(filename);
}
