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
    std::vector<int> neighbors; // edges_from;
    std::vector<int> edges_to;
};

struct Edge {
    int id;
    int from;
    int to;
    double weight;
    double capacity;
    double forman_curvature;
    double ollivier_curvature;

    //Edge(int f) : id(f), from(f), to(f), weight(f), capacity(f) {}

    void invert() {
        std::swap(from, to);
        weight *= -1;
        capacity = 0;
    }
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

    bool node_weighted;
    bool edge_weighted;

    void _shortest_paths(int v0_id, std::vector<double> &dist, std::vector<Edge> &path) {
        dist.assign(nodes.size(), MAX_DIST);
        dist[v0_id] = 0;
        std::vector<bool> inq(nodes.size(), false);
        std::queue<int> q;
        q.push(v0_id);
        path.assign(nodes.size(), {});

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            inq[u] = false;
            for (int &edge_id : nodes[u].neighbors) {
                Edge edge = edges[edge_id];
                if (edge.capacity > EPS && dist[edge.to] - EPS > dist[u] + edge.weight) {
                    dist[edge.to] = dist[u] + edge.weight;
                    path[edge.to] = edge;
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
        for (Node n : nodes) {
            gr.nodes[n.id] = {n.id, n.weight, {}, {}};
        }
        for (Edge e : edges) {
            gr.edges.push_back(e);
            gr.edges.back().id *= 2;
            gr.nodes[e.from].neighbors.push_back(gr.edges.back().id);
            gr.edges.push_back(inverted(gr.edges.back()));
            gr.nodes[e.to].neighbors.push_back(gr.edges.back().id);
        }

        double flow = 0;
        double cost = 0;
        std::vector<double> dist;
        std::vector<Edge> path;
        while (flow < max_flow) {
            gr._shortest_paths(s, dist, path);
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

        //if (flow < max_flow)
        //    return -1;
        //else
            return cost;
    }

    void forman() {
        _fix_graph();
        for (Edge &edge : edges) {
            edge.forman_curvature = _forman_edge(edge.id);
        }
    }

    void ollivier(double idleness=0) {
        Graph gr;
        gr.nodes = nodes;
        gr.edges = edges;
        gr._fix_graph();
        //std::vector<std::vector<double>> paths = floyd_warshall(*this);
        for (Edge &edge : edges) {
            int v1_id = edge.from, v2_id = edge.to;
            std::vector<double> mu1 = _create_mu(v1_id, idleness),
                                mu2 = _create_mu(v2_id, idleness);
            int e_size = edges.size(), n_size = nodes.size();
            Node s = {n_size, 1};
            for (int edge_id : nodes[v1_id].neighbors) {
                Edge e = {(int)gr.edges.size(), s.id, edges[edge_id].to, 0,
                          mu1[edges[edge_id].to]};
                gr.edges.push_back(e);
                s.neighbors.push_back(e.id);
                //gr.nodes[gr.edges[edge_id].to].edges_to.push_back(e.id);
            }
            Node t = {n_size + 1, 1};
            for (int edge_id : nodes[v2_id].neighbors) {
                Edge e = {(int)gr.edges.size(), edges[edge_id].to, t.id, 0,
                          mu2[edges[edge_id].to]};
                gr.edges.push_back(e);
                gr.nodes[edges[edge_id].to].neighbors.push_back(e.id);
                //t.edges_to.push_back(e.id);
            }
            Edge e = {(int)gr.edges.size(), s.id, v1_id, 0, mu1[v1_id]};
            gr.edges.push_back(e);
            s.neighbors.push_back(e.id);
            //gr.nodes[v1_id].edges_to.push_back(e.id);
            gr.nodes.push_back(s);
            e = {(int)gr.edges.size(), v2_id, t.id, 0, mu2[v2_id]};
            gr.edges.push_back(e);
            gr.nodes[v2_id].neighbors.push_back(e.id);
            t.edges_to.push_back(e.id);
            gr.nodes.push_back(t);

            double wd = gr.min_cost_flow(n_size, n_size + 1, 1);
            edge.ollivier_curvature = 1 - wd / edge.weight;

            while (gr.edges.size() > e_size) gr.edges.pop_back();
            while (gr.nodes.size() > n_size) gr.nodes.pop_back();
            //gr.nodes[v1_id].edges_to.pop_back();
            //for (int edge_id : nodes[v1_id].neighbors) {
            //    gr.nodes[edges[edge_id].to].edges_to.pop_back();
            //}
            gr.nodes[v2_id].neighbors.pop_back();
            for (int edge_id : nodes[v2_id].neighbors) {
                gr.nodes[edges[edge_id].to].neighbors.pop_back();
            }
        }
    }
private:
    void _fix_graph() {
        if (!node_weighted) {
            for (Node &node : nodes) {
                node.weight = 1;
            }
        }
        if (!edge_weighted) {
            for (Edge &edge : edges) {
                edge.weight = 1;
            }
        }
    }

    double _forman_edge(int edge_id) const {
        Edge edge = edges[edge_id];
        int v1 = edge.from, v2 = edge.to;
        Node node1 = nodes[v1], node2 = nodes[v2];

        double we = edge.weight;
        double wv1 = node1.weight, wv2 = node2.weight;

        double f = (wv1 + wv2) / we;
        for (int ev1 : node1.neighbors) {
            if (nodes[edges[ev1].to].id != node2.id) {
                f -= wv1 / sqrt(we * edges[ev1].weight);
            }
        }
        for (int ev2 : node2.neighbors) {
            if (nodes[edges[ev2].to].id != node1.id) {
                f -= wv2 / sqrt(we * edges[ev2].weight);
            }
        }
        f *= we;
        return  f;
    }

    std::vector<double> _create_mu(int node_id, double idleness) const {
        double spread = (1 - idleness) / nodes[node_id].neighbors.size();
        std::vector<double> mu(nodes.size(), 0);
        for (int id : nodes[node_id].neighbors) {
            mu[edges[id].to] = spread;
        }
        mu[node_id] = idleness;
        return mu;
    }

    std::vector<std::vector<Node>> data;
};

std::vector<std::vector<double>> floyd_warshall(Graph &gr) {
    std::vector<std::vector<double>> A(gr.nodes.size(),
            std::vector<double>(gr.nodes.size(), MAX_DIST));
    for (size_t i = 0; i < gr.nodes.size(); ++i) {
        A[i][i] = 0;
    }
    for (Node &v1 : gr.nodes) {
        for (int ev1_id : v1.neighbors) {
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
        g.edges.push_back({i * 2, v1, v2, 1, MAX_FLOW, 0, 0});
        g.nodes[v1].neighbors.push_back(i * 2);
        //g.nodes[v2].edges_to.push_back(i * 2);
        g.edges.push_back({i * 2 + 1, v2, v1, 1, MAX_FLOW, 0, 0});
        g.nodes[v2].neighbors.push_back(i * 2 + 1);
        //g.nodes[v1].edges_to.push_back(i * 2 + 1);
    }
    g.forman();
    g.ollivier();

    for (Edge e : g.edges) std::cout << e.from << ' ' << e.to << ' ' << e.ollivier_curvature << '\n';

    return 0;
}
/*
5 5
0 1
0 3
1 2
1 3
3 4
*/
