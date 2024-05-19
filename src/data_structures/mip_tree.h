// Copyright 2024 the MiniMIP Project
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef SRC_DATA_STRUCTURES_MIP_TREE_H_
#define SRC_DATA_STRUCTURES_MIP_TREE_H_

#include <utility>
#include <vector>

#include "ortools/base/logging.h"
#include "ortools/base/strong_vector.h"
#include "ortools/util/strong_integers.h"
#include "src/data_structures/mip_types.h"
#include "src/lp_interface/lpi.h"

namespace minimip {

// All nodes in the search tree are indexed with a strong index type.
DEFINE_STRONG_INDEX_TYPE(NodeIndex);
constexpr NodeIndex kRootNode(0);
constexpr NodeIndex kInvalidNode(-1);

// Spacial value for depth to indicate a node is free (see NodeData below).
constexpr int kFreeNodeDepth = -1;

// As of 2022/08/29, we use the array-of-struct approach to store the nodes'
// data. In comparison to struct-of-array approach, it incurs a larger memory
// footprint due to memory padding and no packing (of boolean members).
// However, regarding speed performance, it is unclear which one is better
// (array-of-struct has better locality). Thus, we opted for the easier to
// implement (and less bug prone) approach (revisit the decision if needed).
//
// We say a node is active (non-free) if its NodeData stores data of an
// actual node in the search tree. We say a node is free if its NodeData
// doesn't store data for any node in the search tree. We maintain a list of
// free nodes to dynamically reclaim them and overwrite with active nodes
// (i.e., reuse already allocated memory).
struct NodeData {
  // For an active node, this is the index to its parent node (for the root
  // node it is equal to kInvalidNode). For a free node, this is the index
  // to the next node in the free list (for the last free node on the list
  // is equal to kInvalidNode).
  NodeIndex parent = kInvalidNode;

  // Depth of this node. It is 0 for root node, strictly positive for any
  // non-root active node, and `kFreeNodeDepth` for any free node.
  int depth = kFreeNodeDepth;

  // Number of cuts added to the node.
  int n_cuts_added = 0;

  // Number of separation rounds.
  int n_separation_rounds = 0;

  // If this node is *not* a root node then it was created by branching from
  // some parent. We need to include this branching decision before solving the
  // LP relaxation. In particular, we set new bounds on the branched variable.
  // For `branch_down`, we will impose the following *upper* bound:
  //   branch_variable <= std::floor(branch_primal_value_in_parent_lp)
  // wheres for `!branch_down`, we will impose the following *lower* bound:
  //   branch_variable >= std::ceil(branch_primal_value_in_parent_lp)
  //
  // TODO(lpawel): Consider adding `struct BranchingData` to wrap these.
  //
  // Index of the branching variable (while branching from the parent to this
  // node).
  ColIndex branch_variable = kInvalidCol;

  // Branching direction (while branching from the parent to this node).
  bool branch_down = false;

  // The primal value of `branch_variable` in the parent's LP relaxation.
  double branch_primal_value_in_parent_lp = 0.0;

  // The objective value of this node's LP relaxation. Note, we always assume
  // minimization in MiniMip. Hence, if `objective_value` equals -kInfinity then
  // the node is unbounded (should never happen, because if the root relaxation
  // was bounded, so should all nodes). If `objective_value` equals +kInfinity
  // then the node is infeasible.
  double lp_objective_value = -kInfinity;

  // Implied variable bounds newly discovered in this node (usually, computed
  // before the LP relaxation via bound propagator). These bounds are inferred
  // from the branching decisions and the implied bounds of all ancestor nodes.
  // To get variable lower/upper bounds that take into account all branching
  // decisions and all implied bounds in the intermediate nodes (up to root) use
  // `RetrieveLowerBounds()` and `RetrieveUpperBounds()`.
  std::vector<ColAndValue> implied_lower_bounds;
  std::vector<ColAndValue> implied_upper_bounds;
};

// This is the main registry of all nodes in the search tree.
class MipTree {
 public:
  // The tree is initialized with an "empty" root node and one free node.
  MipTree();

  // This adds a single branch to the parent node.
  NodeIndex AddNodeByBranchingFromParent(
      NodeIndex parent, ColIndex branch_variable, bool branch_down,
      double branch_primal_value_in_parent_lp);

  // This marks the node `n` as closed, which means the node has been pruned,
  // its LP relaxation is infeasible, or the two child branches have been added.
  // This function also attempts to reclaim the closed node (and iteratively all
  // ancestors up to the root too). A node can be reclaimed (i.e., turned into
  // a free node to be used again later to store another node data) if all
  // node's children are closed too (if there are any).
  void CloseNodeAndReclaimNodesUpToRootIfPossible(NodeIndex n);

  // If tree is empty it means all nodes have been closed (note -- the tree is
  // initialized with an empty root node, thus it is *not* empty on creation).
  bool TreeIsEmpty() const {
    VLOG(10) << "calling TreeIsEmpty().";
    // If the next free node is the root, it means we have reclaimed the root,
    // (i.e., all nodes have been closed).
    return next_free_node_ == kRootNode;
  }

  // A getter of node `n`.
  // WARNING: the underlying nodes' container is not stable on growth. Thus, the
  // reference may be invalidated by adding new nodes.
  // TODO(lpawel): Implement scoped nodes and CHECK() that no new nodes are
  // added when there exist a scoped reference to a anode.
  const NodeData& node(NodeIndex n) const {
    CHECK_GE(n, kRootNode);
    CHECK_LT(n.value(), nodes_.size());
    CHECK(NodeIsActive(n));
    return nodes_[n];
  }

  // Sets some (needed for later) data about the node's LP relaxation.
  // Typically, this is called when exploring a node.
  // TODO(lpawel): When expanding to store extra LP data in nodes, wrap all of
  // them into `struct LpRelaxationData` and refactor this function.
  void SetLpRelaxationDataInNode(NodeIndex n, double lp_objective_value);

  // Sets implied bounds computed for the node. Typically, this is called when
  // exploring a node, before solving its LP relaxation.
  void SetImpliedVariableBoundsInNode(NodeIndex n,
                                      std::vector<ColAndValue> lower_bounds,
                                      std::vector<ColAndValue> upper_bounds);

  // For each variable, merges its lower bound imposed by all branchings
  // and all implied bounds from node `n` up to root node.
  // TODO(lpawel): Consider a sparse bound API.
  DenseRow RetrieveLowerBounds(NodeIndex n, DenseRow root_lower_bounds) const;

  // For each variable, merges its upper bound imposed by all branchings
  // and all implied bounds from node `n` up to root node.
  DenseRow RetrieveUpperBounds(NodeIndex n, DenseRow root_upper_bounds) const;

  // Number of closed (aka processed, visited) nodes so far. A node is closed if
  // all work for that node has been finished. A closed node may have descendent
  // nodes that are not yet closed (due to branching), and thus a closed node is
  // still stored in the internal storage. If a closed node has no descendants
  // (it is infeasible or pruned) it is *reclaimed* (completely removed from
  // the internal storage). If all descendent nodes of the focal node have been
  // reclaimed, the focal node is also reclaimed. Thus, the number of closed
  // nodes is always strictly larger than the number of reclaimed nodes as long
  // as the search is in progress. The two become equal when the search tree is
  // fully explored.
  // TODO(lpawel): Track other statistics like the number of added nodes,
  // reclaimed nodes, pruned nodes, infeasible nodes, etc. Perhaps create a
  // dedicated "stat"class for this (and for other components)?
  int num_closed_nodes() const { return num_closed_nodes_; }

  // Number of added nodes, excluding the root node.
  int num_added_nodes() const { return nodes_.size() - 1; }

  // Get the number of open child nodes for a given node.
  int num_open_children(NodeIndex n) const {
    VLOG(10) << "calling num_open_children().";
    return number_of_open_child_nodes_[n];
  }

  // Get the next free node index.
  NodeIndex get_next_free_node() const { return next_free_node_; }

  void GenerateDotFile(const NodeIndex& current_node,
                       const std::string& filepath);

 private:
  // A helper function used in checks.
  bool NodeIsActive(NodeIndex n) const {
    VLOG(10) << "calling NodeIsActive().";
    DCHECK_GE(n, kRootNode);
    DCHECK_LT(n.value(), nodes_.size());
    return nodes_[n].depth >= 0;
  }

  // Container of all nodes (both active and free).
  absl::StrongVector<NodeIndex, NodeData> nodes_;

  // Counter of the number of child nodes for each node. This is needed to
  // dynamically reclaim elements from `nodes_` when marking a node as closed.
  absl::StrongVector<NodeIndex, int> number_of_open_child_nodes_;

  // The index to the head of free nodes' list.
  NodeIndex next_free_node_;

  // The number of "closed" nodes. A node is closed if it has been visited and
  // processed. In practice, it means that all work (e.g., LP relaxation, bound
  // propagation, cut generation, heuristics) for a node has been finished, and
  // the node turned out infeasible, prunded, or a branching decision has been
  // made and two child branches have been added. Note, this is *not* the same
  // as the number of nodes that have been "reclaimed" (to reclaim a node, it
  // must be closed and all its descendent nodes must be reclaimed).
  int num_closed_nodes_ = 0;
};

}  // namespace minimip

#endif  // SRC_DATA_STRUCTURES_MIP_TREE_H_
