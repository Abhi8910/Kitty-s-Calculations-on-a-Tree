#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <memory>
#include <list>
using namespace std;

struct Node {
  int data;
  list<shared_ptr<Node>> children;
  shared_ptr<Node> parent;
};

class LCA {
public:
  LCA(shared_ptr<Node> root, int N) {
    m_depth.resize(N);
    m_first.resize(N);
    dfs(root, 0);
    
    cerr << "m_depth: [";
    for (auto el: m_depth)
      cerr << el << ' ';
    cerr << ']' << endl;

    cerr << "m_traversal_order: [";
    for (auto el: m_traversal_order)
      cerr << el->data << ' ';
    cerr << ']' << endl;

    cerr << "m_first: [";
    for (auto el: m_first)
      cerr << el << ' ';
    cerr << ']' << endl;

    auto K = m_traversal_order.size();
    m_segment_tree.resize(K * 4);
    fill_segment_tree(1, 0, K - 1);

    cerr << "m_segment_tree: [";
    for (auto el: m_segment_tree)
      cerr << (el?el->data:-1) << ' ';
    cerr << ']' << endl;
  }

  shared_ptr<Node> operator()(shared_ptr<Node> f, shared_ptr<Node> s) const {
    auto search_l = m_first[f->data-1];
    auto search_r = m_first[s->data-1];
    if (search_l > search_r)
      swap(search_l, search_r);
    auto tree_r = m_traversal_order.size() - 1;
    int idx;
    return find_min(1, 0, tree_r, search_l, search_r, idx);
  }

  shared_ptr<Node> operator()(
    const vector<shared_ptr<Node>>::iterator begin,
    const vector<shared_ptr<Node>>::iterator end,
    vector<shared_ptr<Node>>::iterator &split) const {

    auto search_l = m_first[(*begin)->data-1];
    auto search_r = m_first[(*(end-1))->data-1];
    auto tree_r = m_traversal_order.size() - 1;
    int idx;
    auto res = find_min(1, 0, tree_r, search_l, search_r, idx);
    // split the input vector such that 
    // m_first[split->data-1] >= idx and
    // it does not hold for any preceding node
    split = lower_bound(begin, end, idx,
      [this](shared_ptr<Node> n, int idx){
        return this->m_first[n->data-1] < idx;
      });
    return res;
  }

  void sort_nodes(vector<shared_ptr<Node>>::iterator begin, vector<shared_ptr<Node>>::iterator end) const {
    sort(begin, end,
      [this](shared_ptr<Node> f, shared_ptr<Node> s)
      {return this->m_first[f->data-1] < this->m_first[s->data-1];});
  }
private:
  // height of the node
  vector<int> m_depth;
  // order of traversal
  vector<shared_ptr<Node>> m_traversal_order;
  // an index of the first occurance of a node in the m_traversal_order
  vector<int> m_first;
  // a segment tree
  vector<shared_ptr<Node>> m_segment_tree;

  void dfs(shared_ptr<Node> root, int h) {
    auto idx = root->data-1;
    m_depth[idx] = h;
    m_first[idx] = m_traversal_order.size();
    m_traversal_order.push_back(root);
    for (auto child: root->children) {
      dfs(child, h+1);
      m_traversal_order.push_back(root);
    }
  }

  void fill_segment_tree(int i, int l, int r) {
    if (l == r) {
      m_segment_tree[i] = m_traversal_order[l];
      return;
    }
    auto k = (l + r) / 2;
    auto f = 2 * i;
    auto s = 2 * i + 1;
    fill_segment_tree(f, l, k);
    fill_segment_tree(s, k + 1, r);
    auto f_node = m_segment_tree[f];
    auto s_node = m_segment_tree[s];
    if (m_depth[f_node->data-1] < m_depth[s_node->data-1])
      m_segment_tree[i] = f_node;
    else
      m_segment_tree[i] = s_node;
  }

  shared_ptr<Node> find_min(int i, int tree_l, int tree_r, int search_l, int search_r, int &idx) const {
    if (tree_l == search_l && tree_r == search_r) {
      idx = tree_l;
      return m_segment_tree[i];
    }
    auto k = (tree_l + tree_r) / 2;
    if (k < search_l)
      return find_min(2*i+1, k+1, tree_r, search_l, search_r, idx);
    if (k >= search_r)
      return find_min(2*i, tree_l, k, search_l, search_r, idx);
    int l_idx, r_idx;
    auto l_res = find_min(2*i, tree_l, k, search_l, k, l_idx);
    auto r_res = find_min(2*i+1, k+1, tree_r, k+1, search_r, r_idx);
    if (m_depth[l_res->data-1] < m_depth[r_res->data-1]) {
      idx = l_idx;
      return l_res;
    }
    idx = r_idx;
    return r_res;
  }
};

struct QueryNode {
  QueryNode(shared_ptr<Node> b): base(b) {}
  shared_ptr<Node> base;
  list<shared_ptr<QueryNode>> children;
};

void construct_query_tree(
  const vector<shared_ptr<Node>>::iterator begin,
  const vector<shared_ptr<Node>>::iterator end,
  const LCA &lca,
  shared_ptr<QueryNode> &root
  ) {
  if (begin == end)
    return;
  if (end - begin == 1) { // the single node
    root->children.push_back(make_shared<QueryNode>(*begin));
    return;
  }

  vector<shared_ptr<Node>>::iterator split;
  shared_ptr<Node> cur_lca = lca(begin, end, split);

  shared_ptr<QueryNode> cur_root;
  if (root->base == cur_lca) // the same root
    cur_root = root;
  else {
    cur_root = make_shared<QueryNode>(cur_lca);
    root->children.push_back(cur_root);
  }

  if (*split == cur_lca) { // the query node is in a set of LCA
    construct_query_tree(begin, split, lca, cur_root);
    cur_root->children.push_back(make_shared<QueryNode>(cur_lca));
    construct_query_tree(split+1, end, lca, cur_root);
  } else {
    construct_query_tree(begin, split, lca, cur_root);
    construct_query_tree(split, end, lca, cur_root);
  }
}

int main() {
  int N, Q;
  cin >> N >> Q;
  // read nodes and edges
  vector<shared_ptr<Node>> nodes(N);
  for (auto i = 0; i < N; ++i) {
    nodes[i] = make_shared<Node>();
    nodes[i]->data = i+1;
  }
  for (auto i = 0; i < N - 1; ++i) {
    int f, s;
    cin >> f >> s;
    nodes[f-1]->children.push_back(nodes[s-1]);
    nodes[s-1]->children.push_back(nodes[f-1]);
  }
  // choose root
  shared_ptr<Node> root = nodes[0];
  vector<shared_ptr<Node>> to_process(1, root);
  while (!to_process.empty()) {
    auto c = to_process.back();
    to_process.pop_back();
    for (auto el: c->children) {
      el->parent = c;
      el->children.remove(c);
      to_process.push_back(el);
    }
  }
  LCA lca(root, N);
  while (Q--) {
    int K, i;
    cin >> K;
    vector<shared_ptr<Node>> query_nodes;
    while (K--) {
      cin >> i;
      query_nodes.push_back(nodes[i-1]);
    }
    //construct query tree from query nodes and their LCAs
    lca.sort_nodes(query_nodes.begin(), query_nodes.end());
    auto root = make_shared<QueryNode>(nullptr);
    construct_query_tree(query_nodes.begin(), query_nodes.end(), lca, root);

    root = root->children.front();
    if (root->children.empty()) {
      cout << 0 << endl;
      continue;
    }
    list<shared_ptr<QueryNode>> to_process(1, root);
    while (!to_process.empty()) {
      auto cur_node = to_process.front();
      to_process.pop_front();
      cerr << cur_node->base->data << ": [";
      for (auto c: cur_node->children)
        cerr << c->base->data << ' ';
      cerr << "]" << endl;
    }
  }
  return 0;
}

