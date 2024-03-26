// Based on https://www.geeksforgeeks.org/all-topological-sorts-of-a-directed-acyclic-graph/

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <vector>
#include <iostream>
#include <list>
using namespace std;

class Graph {
  int n_items;
  list<int> *adj;
  vector<int> indegree;
  void alltopologicalSortUtil(vector<int>& res, bool visited[]);

public:
  Graph(int n_items);
  void addEdge(int v, int w);
  void alltopologicalSort();
  std::vector<std::vector<int>> m;
};

Graph::Graph(int n_items) : n_items { n_items } {
  adj = new list<int>[n_items];
  for (int i = 0; i < n_items; i++) indegree.push_back(0);
}

void Graph::addEdge(int v, int w) {
  adj[v].push_back(w);
  indegree[w]++;
}

void Graph::alltopologicalSortUtil(vector<int>& res, bool visited[]) {
  bool flag = false;

  for (int i = 0; i < n_items; i++) {
    if (indegree[i] == 0 && !visited[i]) {
      list<int>:: iterator j;
      for (j = adj[i].begin(); j != adj[i].end(); j++)
        indegree[*j]--;

      res.push_back(i);
      visited[i] = true;
      alltopologicalSortUtil(res, visited);

      visited[i] = false;
      res.erase(res.end() - 1);
      for (j = adj[i].begin(); j != adj[i].end(); j++)
        indegree[*j]++;

      flag = true;
    }
  }

  if (!flag){
    m.push_back(res);
  }
}

void Graph::alltopologicalSort() {
  bool *visited = new bool[n_items];
  for (int i = 0; i < n_items; i++)
    visited[i] = false;

  vector<int> res;
  alltopologicalSortUtil(res, visited);
}


// [[Rcpp::export]]
arma::imat all_topological_sorts_cpp(arma::imat prefs, int n_items) {
  Graph g(n_items);
  for(size_t i{}; i < prefs.n_rows; i++) {
    g.addEdge(prefs.at(i, 1) - 1, prefs.at(i, 0) - 1);
  }
  g.alltopologicalSort();

  arma::imat m(g.m.size(), n_items);
  for(size_t i{}; i < m.n_rows; i++) {
    for(size_t j{}; j < m.n_cols; j++) {
      m(i, j) = g.m[i][j] + 1;
    }
  }

  return m;
}
