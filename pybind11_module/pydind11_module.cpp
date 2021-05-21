// #include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <queue>

double EPS = 1e-8;
double MAX_DIST = 1e9;
double MAX_FLOW = 1e9;

using namespace std;

// graph - граф в виде матрицы смежности,
// gr - граф в виде списков исходящих ребер,
// g - в виде списка ребер со всей информацией

struct Edge {
    int to;
    int from;
    double weight;
    double capacity;
};

Edge inverted(Edge edge) {
    std::swap(edge.from, edge.to);
    edge.weight *= -1;
    edge.capacity = 0;
    return edge;
}

// Стандартный алгоритм флойда-уоршелла; принимает граф, возвращает матрицу расстояний между парами вершин
std::vector<std::vector<double>> floyd_warshall(vector<vector<double>> &graph) {
    std::vector<std::vector<double>> A(graph.size(),std::vector<double>(graph.size(), MAX_DIST));
    for (size_t i = 0; i < graph.size(); ++i) {
        A[i][i] = 0;
    }
    for (size_t i = 0; i < graph.size(); ++i) {
        for (size_t j = 0; j < graph.size(); ++j) {
            if (graph[i][j] < EPS) {
                A[i][j] = MAX_DIST;
            } else {
                A[i][j] = graph[i][j];
            }
        }
    }
    for (size_t k = 0; k < graph.size(); ++k)
        for (size_t i = 0; i < graph.size(); ++i)
            for (size_t j = 0; j < graph.size(); ++j)
                if (A[i][k] < MAX_DIST && A[k][j] < MAX_DIST) {
                    if (A[i][k] + A[k][j] < A[i][j] - EPS)
                        A[i][j] = A[i][k] + A[k][j];
                }
    return A;
}

// Записывает в вектор dist кратчайшие расстояния от v0 до остальных вершин,
// в paths - вершину, из которой идет в текущую ребро кратчайшего пути
void shortest_paths(std::vector<std::vector<int>> &adj, std::vector<Edge> &g,
                    int v0, std::vector<double> &dist, std::vector<int> &paths) {
    dist.assign(adj.size(), MAX_DIST);
    dist[v0] = 0;
    std::vector<bool> inq(adj.size(), false);
    std::queue<int> q;
    q.push(v0);
    paths.assign(adj.size(), -1);

    while (!q.empty()) {
        int u = q.front();
        q.pop();
        inq[u] = false;
        for (int e_id : adj[u]) {
            Edge e = g[e_id];
            if (e.capacity > EPS && dist[e.to] - EPS > dist[u] + e.weight) {
                dist[e.to] = dist[u] + e.weight;
                paths[e.to] = e_id;
                if (!inq[e.to]) {
                    inq[e.to] = true;
                    q.push(e.to);
                }
            }
        }
    }
}

// Возвращает минимальную стоимость максимального потока из вершины s в вершину t графа gr
double min_cost_flow(std::vector<std::vector<Edge>> &gr, int s, int t, double max_flow=MAX_FLOW) {
    std::vector<Edge> g;
    std::vector<std::vector<int>> adj(gr.size());
    for (std::vector<Edge> &v : gr) {
        for (Edge &e : v) {
            adj[e.from].push_back(g.size());
            g.push_back(e);
            adj[e.from].push_back(g.size());
            g.push_back(inverted(g.back()));
        }
    }

    double flow = 0;
    double cost = 0;
    std::vector<double> dist;
    std::vector<int> path;
    while (flow < max_flow - EPS) {
        shortest_paths(adj, g, s, dist, path);
        if (dist[t] == MAX_DIST)
            break;

        // find max flow on that path
        double f = max_flow - flow;
        int cur = t;
        while (cur != s) {
            Edge e = g[path[cur]];
            f = fmin(f, e.capacity);
            cur = e.from;
        }

        // apply flow
        flow += f;
        cost += f * dist[t];
        cur = t;
        while (cur != s) {
            g[path[cur]].capacity -= f;
            g[path[cur] ^ 1].capacity += f;
            cur = g[path[cur]].from;
        }
    }

    return cost;
}

// Транстпортное расстояние между mu1 и mu2
double fOT(std::vector<double> &mu1, std::vector<double> &mu2, std::vector<std::vector<double>> &dists) {
    // Создаем граф из двух вершин - истока и стока
    std::vector<std::vector<Edge>> gr(2);
    int s = 0, t = 1;
    // Векторы пар вершин двух долей и соответствующих им вершин исходного графа
    std::vector<std::pair<int, int>> part1, part2;

    // Добавили во вспомогательный граф вершины с ненулевым mu1 или mu2 и соответствующие ребра
    for (int i = 0; i < dists.size(); ++i) {
        if (std::abs(mu1[i]) > EPS) {
            part1.push_back({(int)gr.size(), i});
            Edge e = {s, (int)gr.size(), 0, mu1[i]};
            gr[s].push_back(e);
            gr.push_back({});
        }
        if (std::abs(mu2[i]) > EPS) {
            part2.push_back({(int)gr.size(), i});
            Edge e = {(int)gr.size(), t, 0, mu2[i]};
            gr.push_back({e});
        }
    }

    for (auto ids1 : part1) {
        for (auto ids2 : part2) {
            Edge e = {ids1.first, ids2.first, dists[ids1.second][ids2.second], MAX_FLOW};
            gr[ids1.first].push_back(e);
        }
    }

    double ans = min_cost_flow(gr, s, t, 1);

    return ans;
}

// Строят функцию вероятности перехода из v в смежные вершины (для случая исходящих и входящих вершин)
std::vector<double> _create_mu_from(std::vector<std::vector<double>> &graph, int v,
                                    double idleness, std::vector<std::vector<double>> &dists) {
    double C = 0;
    for (int i = 0; i < graph.size(); ++i) {
        if (graph[i][v] > EPS) {
            C += exp(-dists[i][v]);
        }
    }
    double spread = (1 - idleness) / C;
    std::vector<double> mu(graph.size(), 0);
    for (int i = 0; i < graph.size(); ++i) {
        if (graph[i][v] > EPS) {
            mu[i] = spread * exp(-dists[i][v]);
        }
    }
    mu[v] = idleness;
    return mu;
}
std::vector<double> _create_mu_to(std::vector<std::vector<double>> &graph, int v,
                                  double idleness, std::vector<std::vector<double>> &dists) {
    double C = 0;
    for (int i = 0; i < graph.size(); ++i) {
        if (graph[v][i] > EPS) {
            C += exp(-dists[v][i]);
        }
    }
    double spread = (1 - idleness) / C;
    std::vector<double> mu(graph.size(), 0);
    for (int i = 0; i < graph.size(); ++i) {
        if (graph[v][i] > EPS) {
            mu[i] = spread * exp(-dists[v][i]);
        }
    }
    mu[v] = idleness;
    return mu;
}

// Возвращает вектор кривизн Олливье-Риччи для переданного графа с заданным параметром
std::vector<std::vector<double>> calculate_ollivier(std::vector<std::vector<double>> &graph, double k) {
    std::vector<std::vector<double>> distances = floyd_warshall(graph);
    std::vector<std::vector<double>> curvatures(graph.size(), std::vector<double>(graph.size(), 0));

    for (int i = 0; i < graph.size(); ++i) {
        for (int j = 0; j < graph.size(); ++j) {
            if (graph[i][j] > EPS) {
                std::vector<double> mu1 = _create_mu_from(graph, i, k, distances),
                        mu2 = _create_mu_to(graph, j, k, distances);

                double wd = fOT(mu1, mu2, distances);
                curvatures[i][j] = 1 - wd / distances[i][j];
            }
        }
    }

    return curvatures;
}


// возвращает кривизну Формана-Риччи ребра между вершинами v и u графа graph
double forman_edge(std::vector<std::vector<double>> &graph, int v, int u) {
    double we = graph[v][u];
    double wv1 = 1, wv2 = 1;    // полагаем веса вершин равными 1

    double f = (wv1 + wv2) / we;
    for (int i = 0; i < graph.size(); ++i) {                          // to <-> from ???
        if (graph[v][i] > EPS && i != u) {
            f -= wv1 / sqrt(we * graph[v][i]);
        }
    }
    for (int i = 0; i < graph.size(); ++i) {                          // to <-> from ???
        if (graph[i][u] > EPS && i != v) {
            f -= wv1 / sqrt(we * graph[i][u]);
        }
    }
    f *= we;
    return  f;
}

std::vector<std::vector<double>> calculate_forman(std::vector<std::vector<double>> &graph) {
    std::vector<std::vector<double>> curvatures(graph.size(), std::vector<double>(graph.size(), 0));
    for (int i = 0; i < graph.size(); ++i) {
        for (int j = 0; j < graph.size(); ++j) {
            if (graph[i][j] > EPS) {
                curvatures[i][j] = forman_edge(graph, i, j);
            }
        }
    }

    return curvatures;
}


void dfs(std::vector<std::vector<double>> &graph, int v, int cnt, std::vector<bool> &used) {    // cnt - номер текущей компоненты
    used[v] = true;
    // comp[v] = c_num;

    for (int u: graph[v]) {
        if (!used[u] && (graph[v][u] > EPS || graph[u][v] > EPS)) {
            dfs(graph, u, cnt, used);
        }
    }
}

int connected_components(std::vector<std::vector<double>> &graph) {
    int counter = 0;
    std::vector<bool> used(graph.size());

    for (int i = 0; i < graph.size(); i++) {
        if (!used[i]) {     // если вершина не была достижима ни из одной обработанной
            dfs(graph, i, counter, used);
            counter++;
        }
    }

    return counter;
}


// Одна итерация потока Риччи на графе graph с ideleness и коэффициентом изменения веса mult,
// хирургия происходит при увеличении веса ребра в k раз
void ricci_flow(std::vector<std::vector<double>> &graph, double ideleness=1, double mult=1, double k=1) {
    std::vector<std::vector<double>> curvatures = calculate_ollivier(graph, ideleness);
    for (int i = 0; i < graph.size(); ++i) {
        for (int j = 0; j < graph.size(); ++j) {
            graph[i][j] *= 1 - mult * curvatures[i][j];
            if (1 - mult * curvatures[i][j] > k + EPS) {
                graph[i][j] = 0;    // -1 may be?
            }
        }
    }
}



namespace py = pybind11;

PYBIND11_MODULE(ricci_calculator, m) {
    m.doc() = "pybind11 plugin for compute some curvatures";

    m.def("calculate_ollivier", &calculate_ollivier,
          "Computation ollivier-ricci curvatures for all edges in given graph");
    m.def("calculate_forman", &calculate_forman,
          "Computation forman-ricci curvatures for all edges in given graph");
    m.def("connected_components", &connected_components,
          "Calculating the number of connected components in a given graph");
    m.def("ricci_flow", &ricci_flow,
          "Passing the Ricci flow with surgery on a given graph");
}