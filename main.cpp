#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <mpi.h>
#include <metis.h>

using VertexList = std::vector<std::uint32_t>;
using ColorMap = std::unordered_map<std::uint32_t, std::uint32_t>;
using RandomMap = std::unordered_map<std::uint32_t, int>;

struct CSRGraph {
  idx_t num_vertices;
  idx_t num_edges;
  idx_t *adjacent_offset;
  idx_t *adjacent;

  CSRGraph() : num_vertices(0), num_edges(0) {}

  // Return the start of the adjacency list of 'u'
  idx_t const *adj(idx_t u) const { return adjacent + adjacent_offset[u]; }

  // Return the degree of 'u'
  std::size_t degree(idx_t u) const {
    return adjacent_offset[u + 1] - adjacent_offset[u];
  }

  ~CSRGraph() {
    delete[] adjacent_offset;
    delete[] adjacent;
  }
};

struct VertexColorMessage {
  int vertex;
  int color;
};

void ff_color(VertexList::const_iterator begin, VertexList::const_iterator end,
              const CSRGraph& graph, ColorMap& color,
              std::uint32_t& max_color, VertexColorMessage* msg = nullptr) {
  int colored = 0;
  for (auto u = begin; u != end; ++u) {
    std::vector<bool> used_colors(max_color + 1);
    idx_t const * u_adj = graph.adj(*u);
    for (int i = 0; i < graph.degree(*u); ++i) {
      auto found = color.find(u_adj[i]);
      if (found != std::end(color))
        used_colors[found->second] = true;
    }

    std::uint32_t ff = 0;
    for (auto uc: used_colors) {
      if (!uc)
        break;

      ++ff;
    }

    color[*u] = ff;
    max_color = std::max(max_color, ff);
    if (msg)
      msg[colored] = {static_cast<int>(*u), static_cast<int>(ff)};

    ++colored;
  }
}

std::size_t num_superstep(std::uint32_t nc, std::uint32_t superstep_s) {
    std::size_t cs = nc / superstep_s;
    return nc % superstep_s != 0 ? cs + 1 : cs;
}

std::uint32_t spcr_framework(const CSRGraph& graph, idx_t* partitioning,
                             const VertexList& internal, VertexList& boundary,
                             std::uint32_t superstep_s, int rank,
                             int num_procs) {
  ColorMap color;
  // First color internal vertices
  uint32_t max_color = 0;
  ff_color(std::begin(internal), std::end(internal), graph, color, max_color);

  // Initialize random numbers for each boundary vertex not limited to
  // this partition
  RandomMap random;
  for (auto u: boundary) {
    std::srand(u);
    random[u] = std::rand();
    idx_t const * u_adj = graph.adj(u);
    for (int i = 0; i < graph.degree(u); ++i) {
      std::uint32_t v = u_adj[i];
      if (partitioning[u] != partitioning[v]) {
        std::srand(v);
        random[v] = std::rand();
      }
    }
  }

  // Num of uncolored boundary vertices
  std::uint32_t nc = boundary.size();

  // Num of required supersteps at this process
  std::uint32_t ns = num_superstep(nc, superstep_s);

  // Num of max supersteps
  std::uint32_t max_ns = 0;
  MPI_Allreduce(&ns, &max_ns, 1, MPI_UINT32_T, MPI_MAX, MPI_COMM_WORLD);
  auto vc_send_msg = new VertexColorMessage[superstep_s];
  auto vc_recv_msg = new VertexColorMessage*[max_ns];
  for (int s = 0; s < max_ns; ++s)
    vc_recv_msg[s] = new VertexColorMessage[num_procs * superstep_s];

  auto reqs = new MPI_Request[max_ns];

  // Each round
  while (max_ns) {
    // // std::cerr << "Proc " << rank <<
    // //     ": num of boundary vertices to be colored: " << nc << "\n";
    // std::cerr << "Number of supersteps: " << ns << "\n";
    auto beg_it = std::begin(boundary);

    // Each superstep
    for (int s = 0; s < max_ns; ++s) {
      memset(vc_send_msg, -1, sizeof(*vc_send_msg) * superstep_s);
      std::uint32_t beg = s * superstep_s;
      std::uint32_t end = std::min(beg + superstep_s, nc);
      if (beg < nc)
        ff_color(beg_it + beg, beg_it + end, graph, color, max_color,
                 /*msg=*/vc_send_msg);

      MPI_Iallgather(vc_send_msg, 2 * superstep_s, MPI_INT, vc_recv_msg[s],
                     2 * superstep_s, MPI_INT, MPI_COMM_WORLD, reqs + s);
    }

    MPI_Waitall(max_ns, reqs, MPI_STATUS_IGNORE);
    // Process received color information
    for (int s = 0; s < max_ns; ++s) {
      for (int p = 0; p < num_procs; ++p) {
        // No need to process ranks colors
        if (p == rank)
          continue;

        for (int i = 0; i < superstep_s; ++i) {
          VertexColorMessage rm = vc_recv_msg[s][p * superstep_s + i];
          if (rm.vertex == -1)
            break;

          color[rm.vertex] = rm.color;
          max_color = std::max(max_color, color[rm.vertex]);
        }
      }
    }

    // Resolve conflicts and determine vertices to be colored next round
    int new_nc = 0;
    for (int i = 0; i < nc; ++i) {
      std::uint32_t u = boundary[i];
      idx_t const * u_adj = graph.adj(u);
      for (int j = 0; j < graph.degree(u); ++j) {
        std::uint32_t v = u_adj[j];
        if (color[u] == color.at(v) && random.at(v) > random[u]) {
          boundary[new_nc++] = u;
          break;
        }
      }
    }

    nc = new_nc;
    ns = num_superstep(nc, superstep_s);
    MPI_Allreduce(&ns, &max_ns, 1, MPI_UINT32_T, MPI_MAX, MPI_COMM_WORLD);
  }

  delete[] vc_send_msg;
  for (int s = 0; s < max_ns; ++s)
    delete[] vc_recv_msg[s];

  delete[] vc_recv_msg;
  delete[] reqs;
  return max_color;
}

std::pair<VertexList, VertexList> partition_vertices(
    const CSRGraph& graph, idx_t* partitions, int rank) {
  VertexList internal;
  VertexList boundary;

  for (int u = 0; u < graph.num_vertices; ++u) {
    if (partitions[u] != rank)
      continue;

    std::size_t ud = graph.degree(u);
    idx_t const * u_adj = graph.adj(u);
    bool is_internal = true;
    for (int i = 0; i < ud; ++i) {
      if (partitions[u_adj[i]] != rank) {
        is_internal = false;
        break;
      }
    }

    if (is_internal)
      internal.push_back(u);
    else
      boundary.push_back(u);
  }

  auto degree_sorter = [&graph](std::uint32_t u, std::uint32_t v) {
    std::size_t ud = graph.degree(u);
    std::size_t vd = graph.degree(v);
    if (ud == vd)
      return u < v;

    return ud > vd;
  };

  std::sort(std::begin(internal), std::end(internal), degree_sorter);
  std::sort(std::begin(boundary), std::end(boundary), degree_sorter);
  return {internal, boundary};
}

// The input graph is directed, we also convert it to an undirected
// graph here
int read_csr(std::string graph_path, CSRGraph& out_graph) {
  std::ifstream graph_stream(graph_path);
  if (!graph_stream.is_open())
    return 1;

  std::string graph_line;
  std::vector<std::unordered_set<std::uint32_t>> graph;
  while (std::getline(graph_stream, graph_line)) {
    std::stringstream tok_stream(graph_line);
    std::vector<std::string> tokens;
    std::string tok;
    while (tok_stream >> tok)
      tokens.push_back(tok);

    if (tokens.empty() || tokens[0] == "c")
      continue;

    if (tokens[0] == "a") {
      std::uint32_t u = std::stoul(tokens[1]) - 1;
      std::uint32_t v = std::stoul(tokens[2]) - 1;
      bool inserted;
      std::tie(std::ignore, inserted) = graph[u].insert(v);
      if (inserted)
        ++out_graph.num_edges;

      std::tie(std::ignore, inserted) = graph[v].insert(u);
      if (inserted)
        ++out_graph.num_edges;
    } else if (tokens[0] == "p" && tokens[1] == "sp") {
      out_graph.num_vertices = std::stoul(tokens[2]);
      graph.resize(out_graph.num_vertices);
      std::cerr << "Num of vertices: " << out_graph.num_vertices<< "\n";
      out_graph.adjacent_offset = new idx_t[out_graph.num_vertices + 1];
    }
  }

  std::cerr << "Num of edges: " << out_graph.num_edges / 2 << "\n";
  out_graph.adjacent = new idx_t[out_graph.num_edges];
  out_graph.adjacent_offset[0] = 0;
  for (int u = 0; u < out_graph.num_vertices; ++u) {
    out_graph.adjacent_offset[u + 1] = out_graph.adjacent_offset[u] +
                                       graph[u].size();
    int i = 0;
    for (auto& v: graph[u]) {
      out_graph.adjacent[out_graph.adjacent_offset[u] + i] = v;
      ++i;
    }
  }

  return 0;
}

int main(int argc, char* argv[]) {
  if (int err = MPI_Init(&argc, &argv) != 0) {
    std::cerr << "Could not initialize MPI\n";
    return err;
  }

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << "<graph_file> [superstep_size]\n";
    return 1;
  }

  int rank;
  int num_procs;
  int error = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  error = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  CSRGraph graph;
  std::cerr << "reading " << argv[1] << "\n";
  read_csr(argv[1], graph);
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
  options[METIS_OPTION_NUMBERING] = 0;
  idx_t* partitioning = new idx_t[graph.num_vertices];
  int ncon = 1;
  idx_t objval;
  if (num_procs > 1) {
    int metis_err = METIS_PartGraphKway(
        &graph.num_vertices, &ncon, graph.adjacent_offset, graph.adjacent,
        nullptr, nullptr, nullptr, &num_procs, nullptr, nullptr, options,
        &objval, partitioning);

    if (metis_err != METIS_OK)
      std::cerr << "Graph could not be partitioned\n";
  } else {
    std::memset(partitioning, 0, sizeof(idx_t) * graph.num_vertices);
  }

  VertexList internal;
  VertexList boundary;
  std::tie(internal, boundary) = partition_vertices(graph, partitioning, rank);
  std::uint32_t superstep_size = 100;
  if (argc == 3)
    superstep_size = std::stoul(argv[2]);

  std::vector<double> times(10);
  std::uint32_t colors = 0;
  for (int i = 0; i < 10; ++i) {
    VertexList boundary_cp = boundary;
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();
        colors = spcr_framework(graph, partitioning, internal, boundary_cp,
                                superstep_size, rank, num_procs);
    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();
    times[i] = end - start;
  }

  if (rank == 0) {
    double min_time = *std::min_element(std::begin(times), std::end(times));
    std::cerr << "SPCRFramework took " << min_time * 1000.0 << "ms\n";
    std::cerr << "Used " << colors + 1 << " colors\n";
  }

  delete[] partitioning;
  MPI_Finalize();
  return error;
}
