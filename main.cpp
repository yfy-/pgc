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
    std::unordered_set<std::uint32_t> forbidden;
    idx_t const * u_adj = graph.adj(*u);
    for (int i = 0; i < graph.degree(*u); ++i) {
      auto found = color.find(u_adj[i]);
      if (found != std::end(color))
        forbidden.insert(found->second);
    }

    std::uint32_t ff = 0;
    while (true) {
      auto found = forbidden.find(ff);
      if (found == std::end(forbidden))
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

std::size_t num_superstep(std::uint32_t* all_nc_sizes, int num_procs,
                          std::uint32_t superstep_s) {
  std::size_t max_num = 0;
  for (int i = 0; i < num_procs; ++i) {
    std::size_t cs = all_nc_sizes[i] / superstep_s;
    if (all_nc_sizes[i] % superstep_s != 0)
      cs++;

    max_num = std::max(max_num, cs);
  }

  return max_num;
}

std::uint32_t spcr_framework(const CSRGraph& graph, idx_t* partitioning,
                             const VertexList& internal,
                             const VertexList& boundary,
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

  VertexList not_colored = boundary;
  std::uint32_t nc = not_colored.size();
  auto nc_sizes = new std::uint32_t[num_procs];
  MPI_Allgather(&nc, 1, MPI_UINT32_T, nc_sizes, 1, MPI_UINT32_T,
                MPI_COMM_WORLD);
  std::uint32_t ns = num_superstep(nc_sizes, num_procs, superstep_s);
  auto vc_send_msg = new VertexColorMessage[superstep_s];
  auto vc_recv_msg = new VertexColorMessage*[ns];
  for (int s = 0; s < ns; ++s)
    vc_recv_msg[s] = new VertexColorMessage[num_procs * superstep_s];

  auto reqs = new MPI_Request[ns];
  auto stats = new MPI_Status[ns];
  while (ns) {
    std::cerr << "Proc " << rank <<
        ": num of boundary vertices to be colored: " << nc << "\n";
    auto beg_it = std::begin(not_colored);
    for (int s = 0; s < ns; ++s) {
      memset(vc_send_msg, -1, sizeof(*vc_send_msg) * superstep_s);
      std::uint32_t beg = s * superstep_s;
      std::uint32_t end = std::min(beg + superstep_s, nc);
      // std::cerr << "Proc " << rank << ": " << s << "\n";
      // std::cerr << "Proc " << rank << ":" << beg << ", " << end << "\n";
      if (beg < nc)
        ff_color(beg_it + beg, beg_it + end, graph, color, max_color,
                 /*msg=*/vc_send_msg);

      MPI_Iallgather(vc_send_msg, 2 * superstep_s, MPI_INT, vc_recv_msg[s],
                     2 * superstep_s, MPI_INT, MPI_COMM_WORLD, reqs + s);
    }

    MPI_Waitall(ns, reqs, stats);
    for (int s = 0; s < ns; ++s) {
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

    not_colored.clear();
    for (auto u: boundary) {
      idx_t const * u_adj = graph.adj(u);
      for (int i = 0; i < graph.degree(u); ++i) {
        std::uint32_t v = u_adj[i];
        if (color[u] == color.at(v) && random.at(v) > random[u])
          not_colored.push_back(u);
      }
    }

    nc = not_colored.size();
    MPI_Allgather(&nc, 1, MPI_INT, nc_sizes, 1, MPI_INT, MPI_COMM_WORLD);
    ns = num_superstep(nc_sizes, num_procs, superstep_s);
  }

  delete[] nc_sizes;
  delete[] vc_send_msg;
  for (int s = 0; s < ns; ++s)
    delete[] vc_recv_msg[s];

  delete[] vc_recv_msg;
  delete[] reqs;
  delete[] stats;
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
      auto ins_res = graph[u].insert(v);
      if (ins_res.second)
        ++out_graph.num_edges;

      ins_res = graph[v].insert(u);
      if (ins_res.second)
        ++out_graph.num_edges;
    } else if (tokens[0] == "p" && tokens[1] == "sp") {
      out_graph.num_vertices = std::stoul(tokens[2]);
      graph.resize(out_graph.num_vertices);
      std::cerr << "Num of vertices: " << out_graph.num_vertices << "\n";
      out_graph.adjacent_offset = new idx_t[out_graph.num_vertices + 1];
    }
  }

  std::cerr << "Num of edges: " << out_graph.num_edges << "\n";
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
  std::cerr << "Proc " << rank << ": internal vertices size: " <<
      internal.size() << "\n";
  std::cerr << "Proc " << rank <<
      ": boundary vertices size: " << boundary.size() << "\n";
  std::uint32_t superstep_size = 100;
  if (argc == 3)
    superstep_size = std::stoul(argv[2]);

  double start = MPI_Wtime();
  std::uint32_t colors = spcr_framework(graph, partitioning, internal,
                                        boundary, superstep_size, rank, num_procs);
  std::cerr << "Proc " << rank << ": SPCRFramework took " <<
      (MPI_Wtime() - start) * 1000.0 << "ms\n";
  delete[] partitioning;
  std::cerr << "Used " << colors + 1 << " colors\n";
  MPI_Finalize();
  std::cerr << "Proc " << rank << " exiting with " << error << "\n";
  return error;
}
