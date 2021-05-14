#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <queue>
#include <windef.h>

double EPS = 1e-8;
double MAX_DIST = 1e9;
double MAX_FLOW = 1e9;


struct Node {
    int id;
    double weight;
    std::vector<int> from;
    std::vector<int> to;
};

struct Edge {
    int id;
    int from;
    int to;
    double weight;
    double capacity;
};

Edge inverted(Edge edge) {
    std::swap(edge.from, edge.to);
    edge.id += 1;
    edge.weight *= -1;
    edge.capacity = 0;
    return edge;
}

class Graph;
std::vector<std::vector<double>> floyd_warshall(Graph &gr);

class Graph {
public:
    std::vector<Node> nodes;
    std::vector<Edge> edges;

    std::vector<double> forman_curvature;
    std::vector<double> ollivier_curvature;

    std::vector<double> node_forman_curvature;
    std::vector<double> node_ollivier_curvature;

    std::vector<std::vector<double>> distances;

    bool node_weighted;
    bool edge_weighted;

    void shortest_paths(int v0_id, std::vector<double> &dist, std::vector<Edge> &paths) {
        dist.assign(nodes.size(), MAX_DIST);
        dist[v0_id] = 0;
        std::vector<bool> inq(nodes.size(), false);
        std::queue<int> q;
        q.push(v0_id);
        paths.assign(nodes.size(), {});

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            inq[u] = false;
            for (int &edge_id : nodes[u].from) {
                Edge edge = edges[edge_id];
                if (edge.capacity > EPS && dist[edge.to] - EPS > dist[u] + edge.weight) {
                    dist[edge.to] = dist[u] + edge.weight;
                    paths[edge.to] = edge;
                    if (!inq[edge.to]) {
                        inq[edge.to] = true;
                        q.push(edge.to);
                    }
                }
            }
        }
    }

    double min_cost_flow(int s, int t, double max_flow=MAX_FLOW) {
        Graph gr;
        gr.nodes = std::vector<Node>(nodes.size());
        for (const Node &n : nodes) {
            gr.nodes[n.id] = {n.id, n.weight, {}, {}};
        }
        for (Edge e : edges) {
            gr.edges.push_back(e);
            gr.edges.back().id *= 2;
            gr.nodes[e.from].from.push_back(gr.edges.back().id);
            gr.edges.push_back(inverted(gr.edges.back()));
            gr.nodes[e.to].from.push_back(gr.edges.back().id);
        }

        double flow = 0;
        double cost = 0;
        std::vector<double> dist;
        std::vector<Edge> path;
        while (flow < max_flow - EPS) {
            gr.shortest_paths(s, dist, path);
            if (dist[t] == MAX_DIST)
                break;

            // find max flow on that path
            double f = max_flow - flow;
            int cur = t;
            while (cur != s) {
                f = fmin(f, path[cur].capacity);
                cur = path[cur].from;
            }

            // apply flow
            flow += f;
            cost += f * dist[t];
            cur = t;
            while (cur != s) {
                gr.edges[path[cur].id].capacity -= f;
                gr.edges[path[cur].id ^ 1].capacity += f;
                cur = path[cur].from;
            }
        }

        return cost;
    }

    void forman() {
        _fix_graph();
        forman_curvature = std::vector<double>(edges.size(), 0);
        for (Edge &edge : edges) {
            forman_curvature[edge.id] = _forman_edge(edge.id);
        }
    }

    // Транстпортное расстояние между mu1 и mu2
    double fOT(std::vector<double> &mu1, std::vector<double> &mu2, std::vector<std::vector<double>> &dists) {
        Graph gr;
        Node    s = {0, 1, {}, {}}, // добавляем исток
                t = {1, 1, {}, {}}; // и сток
        std::vector<int> mu1_ids, mu2_ids;
        std::vector<std::pair<int, int>> part1, part2;
        gr.nodes.resize(2);
        // Добавили во вспомогательный граф вершины с ненулевым mu1 или mu2 и соответствующие ребра
        for (int i = 0; i < nodes.size(); ++i) {
            if (std::abs(mu1[i]) > EPS) {
                Node n = {(int)gr.nodes.size(), nodes[i].weight, {}, {}};
                Edge e = {(int)gr.edges.size(), s.id, n.id, 0, mu1[i]};
                s.from.push_back(e.id);
                n.to.push_back(e.id);
                gr.nodes.push_back(n);
                gr.edges.push_back(e);
                part1.push_back({n.id, i});
            }
            if (std::abs(mu2[i]) > EPS) {
                Node n = {(int)gr.nodes.size(), nodes[i].weight, {}, {}};
                Edge e = {(int)gr.edges.size(), n.id, t.id, 0, mu2[i]};
                n.from.push_back(e.id);
                t.to.push_back(e.id);
                gr.nodes.push_back(n);
                gr.edges.push_back(e);
                part2.push_back({n.id, i});
            }
        }
        gr.nodes[0] = s;
        gr.nodes[1] = t;

        //std::vector<std::vector<double>> dists = floyd_warshall(*this);
        for (auto ids1 : part1) {
            for (auto ids2 : part2) {
                Edge e = {(int)gr.edges.size(), ids1.first, ids2.first, dists[ids1.second][ids2.second], MAX_FLOW};
                gr.nodes[ids1.first].from.push_back(e.id);
                gr.nodes[ids2.first].to.push_back(e.id);
                gr.edges.push_back(e);
            }
        }

        double ans = gr.min_cost_flow(s.id, t.id, 1);

        return ans;
    }

    // Вычисление кривизны Олливье-Риччи для всего графа потоковым методом
    void ollivier(double idleness=0) {
        // Строим граф-копию и дополняем его истоком и стоком (нет)
        Graph gr;
        gr.nodes = nodes;
        gr.edges = edges;
        gr._fix_graph();
        for (Edge &edge : edges) {
            int v1_id = edge.from, v2_id = edge.to;
            std::vector<double> mu1 = _create_mu_from(v1_id, idleness),
                                mu2 = _create_mu_to(v2_id, idleness);

            double wd = fOT(mu1, mu2, distances);
            ollivier_curvature[edge.id] = 1 - wd / distances[v1_id][v2_id]; // edge.weight;
        }
    }

    // Кривизна вершин (оба варианта) ((это несложно, просто среднее арифметическре кривизн инцедентных ребер))
    void node_forman(int type=3) { // 1 - учитываются только исходящие ребра, 2 - только входящие, 3 - все инцидентные
        _fix_graph();
        for (Node &n : nodes) {
            double sum = 0;
            int count = 0;
            if (type & 1) {
                count += n.from.size();
                for (int id : n.from) {
                    sum += forman_curvature[id];
                }
            }
            if (type & 2) {
                count += n.to.size();
                for (int id : n.to) {
                    sum += forman_curvature[id];
                }
            }
            node_forman_curvature[n.id] = sum / count;
        }
    }
    void node_ollivier(int type=3) { // 1 - учитываются только исходящие ребра, 2 - только входящие, 3 - все инцидентные
        _fix_graph();
        for (Node &n : nodes) {
            double sum = 0;
            int count = 0;
            if (type & 1) {
                count += n.from.size();
                for (int id : n.from) {
                    sum += ollivier_curvature[id];
                }
            }
            if (type & 2) {
                count += n.to.size();
                for (int id : n.to) {
                    sum += ollivier_curvature[id];
                }
            }
            node_ollivier_curvature[n.id] = sum / count;
        }
    }

    // Одна итерация потока Риччи
    void ricci_flow(double step=1) {
        std::vector<double> new_weights(edges.size());
        if (ollivier_curvature.empty()) {
            std::cout << "please run the \"ollivier()\" with required idleness\n";
            return;
        }
        for (Edge &e : edges) {
            e.weight *= 1 - step * ollivier_curvature[e.id];
        }
    }
private:
    // Если граф невзвешенный, все веса полагаем равными 1
    void _fix_graph() {
        if (forman_curvature.empty()) {
            forman_curvature.resize(edges.size());
        }
        if (ollivier_curvature.empty()) {
            ollivier_curvature.resize(edges.size());
        }
        if (node_forman_curvature.empty()) {
            forman_curvature.resize(nodes.size());
        }
        if (node_ollivier_curvature.empty()) {
            ollivier_curvature.resize(nodes.size());
        }

        if (!node_weighted) {
            for (Node &node : nodes) {
                node.weight = 1;
            }
            node_weighted = 1;
        }
        if (!edge_weighted) {
            for (Edge &edge : edges) {
                edge.weight = 1;
            }
            edge_weighted = 1;
        }

        if (distances.empty()) {
            distances = floyd_warshall(*this);
        }
    }

    // возвращает кривизну Формана-Риччи ребра с переданным id
    double _forman_edge(int edge_id) const {
        Edge edge = edges[edge_id];
        int id1 = edge.from, id2 = edge.to;
        Node node1 = nodes[id1], node2 = nodes[id2];

        double we = edge.weight;
        double wv1 = node1.weight, wv2 = node2.weight;

        double f = (wv1 + wv2) / we;
        for (int ev1 : node1.from) {                          // to <-> from ???
            if (nodes[edges[ev1].to].id != node2.id) {
                f -= wv1 / sqrt(we * edges[ev1].weight);
            }
        }
        for (int ev2 : node2.from) {
            if (nodes[edges[ev2].to].id != node1.id) {
                f -= wv2 / sqrt(we * edges[ev2].weight);
            }
        }
        f *= we;
        return  f;
    }

    std::vector<double> _create_mu_from(int node_id, double idleness) const {
        double C = 0;
        for (int id : nodes[node_id].to) {
            C += exp(-distances[edges[id].from][node_id]);
        }
        double spread = (1 - idleness) / C; // nodes[node_id].to.size();
        std::vector<double> mu(nodes.size(), 0);
        for (int id : nodes[node_id].to) {
            mu[edges[id].from] = spread * exp(-distances[edges[id].from][node_id]);
        }
        mu[node_id] = idleness;
        return mu;
    }
    std::vector<double> _create_mu_to(int node_id, double idleness) const {
        double C = 0;
        for (int id : nodes[node_id].from) {
            C += exp(-distances[node_id][edges[id].to]);
        }
        double spread = (1 - idleness) / C; // nodes[node_id].from.size();
        std::vector<double> mu(nodes.size(), 0);
        for (int id : nodes[node_id].from) {
            mu[edges[id].to] = spread * exp(-distances[node_id][edges[id].to]);
        }
        mu[node_id] = idleness;
        return mu;
    }
};

// Стандартный алгоритм флойда-уоршелла; принимает граф, возвращает матрицу расстояний между парами вершин
std::vector<std::vector<double>> floyd_warshall(Graph &gr) {
    std::vector<std::vector<double>> A(gr.nodes.size(),
            std::vector<double>(gr.nodes.size(), MAX_DIST));
    for (size_t i = 0; i < gr.nodes.size(); ++i) {
        A[i][i] = 0;
    }
    for (Node &v1 : gr.nodes) {
        for (int ev1_id : v1.from) {
            A[v1.id][gr.edges[ev1_id].to] = gr.edges[ev1_id].weight;
        }
    }
    for (size_t k = 0; k < gr.nodes.size(); ++k)
        for (size_t i = 0; i < gr.nodes.size(); ++i)
            for (size_t j = 0; j < gr.nodes.size(); ++j)
                if (A[i][k] < MAX_DIST && A[k][j] < MAX_DIST) {
                    if (A[i][k] + A[k][j] < A[i][j] - EPS)
                        A[i][j] = A[i][k] + A[k][j];
                }
    return A;
}

int main() {
    Graph g;
    g.edge_weighted = 0;
    g.node_weighted = 0;
    int n, m;
    std::cin >> n >> m;
    for (int i = 0; i < n; ++i) {
        g.nodes.push_back({i, 1, std::vector<int>(), std::vector<int>()});
    }
    for (int i = 0; i < m; ++i) {
        int v1, v2;
        std::cin >> v1 >> v2;
        g.edges.push_back({i * 2, v1, v2, 1, MAX_FLOW});
        g.nodes[v1].from.push_back(i * 2);
        g.nodes[v2].to.push_back(i * 2);
        g.edges.push_back({i * 2 + 1, v2, v1, 1, MAX_FLOW});
        g.nodes[v2].from.push_back(i * 2 + 1);
        g.nodes[v1].to.push_back(i * 2 + 1);
    }
    g.forman();
    g.ollivier();

    for (Edge e : g.edges) std::cout << e.from << ' ' << e.to << ' ' << g.ollivier_curvature[e.id] << '\n';
    for (Edge e : g.edges) std::cout << e.from << ' ' << e.to << ' ' << g.forman_curvature[e.id] << '\n';

    g.ricci_flow(0.65);
    g.ollivier();
    for (Edge e : g.edges) std::cout << e.from << ' ' << e.to << ' ' << g.ollivier_curvature[e.id] << '\n';

    return 0;
}
