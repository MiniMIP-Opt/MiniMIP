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

#include "src/data_structures/mip_tree.h"

#include <fstream>

namespace minimip {

MipTree::MipTree() {
  VLOG(10) << "calling MipTree().";
  // Let's reserve some memory to reduce memory re-allocations (an empty node is
  // rather cheap memory-wise, hence this is OK).
  nodes_.reserve(1024);
  nodes_.resize(2);
  nodes_[kRootNode].depth = 0;
  next_free_node_ = NodeIndex(1);
  number_of_open_child_nodes_.resize(2, 0);
}

NodeIndex MipTree::AddNodeByBranchingFromParent(
    NodeIndex parent, ColIndex branch_variable, bool branch_down,
    double branch_primal_value_in_parent_lp) {
  VLOG(10) << "calling AddNodeByBranchingFromParent().";
  CHECK_GE(parent, kRootNode);
  CHECK_LT(parent.value(), nodes_.size());
  CHECK(NodeIsActive(parent));

  if (next_free_node_ == kInvalidNode) {
    nodes_.resize(nodes_.size() + 1);
    number_of_open_child_nodes_.resize(nodes_.size(), 0);
    next_free_node_ = NodeIndex(nodes_.size() - 1);
  }

  const NodeIndex n = next_free_node_;
  next_free_node_ = nodes_[next_free_node_].parent;
  CHECK(!NodeIsActive(n));
  nodes_[n].parent = parent;
  nodes_[n].branch_variable = branch_variable;
  nodes_[n].branch_down = branch_down;
  nodes_[n].branch_primal_value_in_parent_lp = branch_primal_value_in_parent_lp;
  nodes_[n].depth = nodes_[nodes_[n].parent].depth + 1;
  // As of 2022/08/29, we assume a binary search tree.
  CHECK_GE(number_of_open_child_nodes_[nodes_[n].parent], 0);
  CHECK_LE(number_of_open_child_nodes_[nodes_[n].parent], 1);
  ++number_of_open_child_nodes_[nodes_[n].parent];
  number_of_open_child_nodes_[n] = 0;
  return n;
}

void MipTree::CloseNodeAndReclaimNodesUpToRootIfPossible(NodeIndex n) {
  VLOG(10) << "calling CloseNodeAndReclaimNodesUpToRootIfPossible().";
  CHECK_GE(n, kRootNode);
  CHECK_LT(n.value(), nodes_.size());
  CHECK(NodeIsActive(n));

  // As of 2022/08/29, we assume a binary search tree, i.e., when marking the
  // node as closed it must either have 0 children (the node was pruned or
  // infeasible), or 2 children (the node generated up and down branches).
  CHECK(number_of_open_child_nodes_[n] == 0 ||
        number_of_open_child_nodes_[n] == 2);
  CHECK_GE(number_of_open_child_nodes_[nodes_[n].parent], 0);
  CHECK_LE(number_of_open_child_nodes_[nodes_[n].parent], 2);

  while (true) {
    if (number_of_open_child_nodes_[n] > 0) return;

    // Retrieve the parent node index
    const NodeIndex parent = nodes_[n].parent;

    // Reset the current node's data to prevent further access
    nodes_[n] = NodeData();
    nodes_[n].parent = next_free_node_;

    // Ensure `parent` is not invalid before accessing its children count
    if (parent != kInvalidNode) {
      --number_of_open_child_nodes_[parent];
    }

    ++num_closed_nodes_;
    next_free_node_ = n;

    // Move up to the parent node
    n = parent;
    if (parent == kInvalidNode) return;
  }
}

void MipTree::SetLpRelaxationDataInNode(NodeIndex n,
                                        double lp_objective_value) {
  VLOG(10) << "calling SetLpRelaxationDataInNode().";
  CHECK_GE(n, kRootNode);
  CHECK_LT(n.value(), nodes_.size());
  CHECK(NodeIsActive(n));
  // Typically, we should only need to update the relaxation in leaf
  // nodes. Hence the check (remove it, if ever needed).
  CHECK_EQ(number_of_open_child_nodes_[n], 0);
  nodes_[n].lp_objective_value = lp_objective_value;
}

void MipTree::SetImpliedVariableBoundsInNode(
    NodeIndex n, std::vector<ColAndValue> lower_bounds,
    std::vector<ColAndValue> upper_bounds) {
  VLOG(10) << "calling SetImpliedVariableBoundsInNode().";
  CHECK_GE(n, kRootNode);
  CHECK_LT(n.value(), nodes_.size());
  CHECK(NodeIsActive(n));
  // Typically, we should only need to update the implied bounds in leaf
  // nodes. Hence the check (remove it, if ever needed).
  CHECK_EQ(number_of_open_child_nodes_[n], 0);
  nodes_[n].implied_lower_bounds = std::move(lower_bounds);
  nodes_[n].implied_upper_bounds = std::move(upper_bounds);
}

DenseRow MipTree::RetrieveLowerBounds(NodeIndex n,
                                      DenseRow root_lower_bounds) const {
  VLOG(10) << "calling RetrieveLowerBounds().";
  CHECK_GE(n, kRootNode);
  CHECK_LT(n.value(), nodes_.size());
  CHECK(NodeIsActive(n));
  DenseRow lower_bounds = std::move(root_lower_bounds);
  for (; n != kInvalidNode; n = nodes_[n].parent) {
    for (const ColAndValue& lb : nodes_[n].implied_lower_bounds) {
      lower_bounds[lb.col] = std::max(lower_bounds[lb.col], lb.value);
    }
    if (!nodes_[n].branch_down and nodes_[n].branch_variable >= 0) {
      lower_bounds[nodes_[n].branch_variable] =
          std::max(lower_bounds[nodes_[n].branch_variable],
                   std::ceil(nodes_[n].branch_primal_value_in_parent_lp));
    }
  }
  return lower_bounds;
}

DenseRow MipTree::RetrieveUpperBounds(NodeIndex n,
                                      DenseRow root_upper_bounds) const {
  VLOG(10) << "calling RetrieveUpperBounds().";
  CHECK_GE(n, kRootNode);
  CHECK_LT(n.value(), nodes_.size());
  CHECK(NodeIsActive(n));
  DenseRow upper_bounds = std::move(root_upper_bounds);
  for (; n != kInvalidNode; n = nodes_[n].parent) {
    for (const ColAndValue& ub : nodes_[n].implied_upper_bounds) {
      upper_bounds[ub.col] = std::min(upper_bounds[ub.col], ub.value);
    }
    if (nodes_[n].branch_down) {
      upper_bounds[nodes_[n].branch_variable] =
          std::min(upper_bounds[nodes_[n].branch_variable],
                   std::floor(nodes_[n].branch_primal_value_in_parent_lp));
    }
  }
  return upper_bounds;
}

}  // namespace minimip
