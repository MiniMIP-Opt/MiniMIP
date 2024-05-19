// Copyright (2024) the MiniMIP Project
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

#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "unit_tests/utils.h"

namespace minimip {

TEST(MipTreeTest, Initialization) {
  MipTree tree;
  EXPECT_FALSE(tree.TreeIsEmpty());
  // Tree should not be empty after initialization
  EXPECT_EQ(tree.node(kRootNode).depth, 0);
  // Root node depth should be 0
  EXPECT_EQ(tree.node(kRootNode).parent, kInvalidNode);
  // Root should have no parent
}

TEST(MipTreeTest, BasicBranchingFromRoot) {
  MipTree tree;

  // Add a node by branching from the root
  NodeIndex root_node = kRootNode;
  NodeIndex child_node =
      tree.AddNodeByBranchingFromParent(root_node, ColIndex(1), true, 1.5);

  // Verify that the new node has the correct parent
  EXPECT_EQ(tree.node(child_node).parent, root_node);
  EXPECT_EQ(tree.node(child_node).branch_variable, ColIndex(1));
  EXPECT_TRUE(tree.node(child_node).branch_down);
  EXPECT_EQ(tree.node(child_node).branch_primal_value_in_parent_lp, 1.5);

  // Check the depth of the new node (should be one more than the root)
  EXPECT_EQ(tree.node(child_node).depth, tree.node(root_node).depth + 1);

  // Ensure the root's open child count has increased
  EXPECT_EQ(tree.num_open_children(root_node), 1);
}

TEST(MipTreeTest, EnforceBinarySearchTree) {
  MipTree tree;

  // Root node index
  NodeIndex root_node = kRootNode;

  // Create the initial free node
  NodeIndex initial_free_node =
      tree.AddNodeByBranchingFromParent(root_node, ColIndex(2), false, 3.0);
  NodeIndex initial_free_node2 =
      tree.AddNodeByBranchingFromParent(root_node, ColIndex(2), true, 3.0);
  EXPECT_DEATH(
      tree.AddNodeByBranchingFromParent(root_node, ColIndex(2), true, 2.0),
      "Check failed: *");
}

TEST(MipTreeTest, BranchingWithInvalidParent) {
  // Initialize a new MipTree
  MipTree tree;

  // Attempt to branch from an invalid parent node
  EXPECT_DEATH(
      tree.AddNodeByBranchingFromParent(NodeIndex(999), ColIndex(1), true, 1.0),
      "Check failed: *");
}

TEST(MipTreeTest, CloseAndReclaimNodeAtRoot) {
  MipTree tree;
  NodeIndex child_node1 =
      tree.AddNodeByBranchingFromParent(kRootNode, ColIndex(1), true, 1.5);

  NodeIndex child_node2 =
      tree.AddNodeByBranchingFromParent(kRootNode, ColIndex(1), false, 1.5);

  // Initially, ensure the tree is not empty
  ASSERT_FALSE(tree.TreeIsEmpty());

  // Close the child node and fail to reclaim the root node.
  tree.CloseNodeAndReclaimNodesUpToRootIfPossible(child_node2);
  ASSERT_FALSE(tree.TreeIsEmpty());
  ASSERT_EQ(tree.num_closed_nodes(), 1);

  // Close the child node and reclaim the root node.
  tree.CloseNodeAndReclaimNodesUpToRootIfPossible(child_node1);
  EXPECT_TRUE(tree.TreeIsEmpty());
}

// This Tests makes sure that the tree is able to close and
// reclaim multiple nodes including the root node. The number of grandchildren
// is important to check the behaviour of the tree when closing multiple nodes,
// including non-root node parents.
TEST(MipTreeTest, CloseAndReclaimMultipleNodesUpToRoot) {
  // Create the MipTree instance
  MipTree tree;

  // Add a single child node branching from the root node
  NodeIndex child_node =
      tree.AddNodeByBranchingFromParent(kRootNode, ColIndex(1), true, 1.5);

  // Add a grandchild node branching from the child node
  NodeIndex grandchild_node1 =
      tree.AddNodeByBranchingFromParent(child_node, ColIndex(2), true, 2.5);

  // Add a grandchild node branching from the child node
  NodeIndex grandchild_node2 =
      tree.AddNodeByBranchingFromParent(child_node, ColIndex(2), false, 2.5);

  // Add a grand-grandchild node branching from the child node
  NodeIndex grand_grandchild_node = tree.AddNodeByBranchingFromParent(
      grandchild_node1, ColIndex(3), true, 3.5);

  // Add a grand-grand-grandchild node branching from the child node
  NodeIndex grand_grand_grandchild_node = tree.AddNodeByBranchingFromParent(
      grand_grandchild_node, ColIndex(4), true, 3.5);

  // Ensure all nodes are active after creation
  ASSERT_EQ(tree.node(kRootNode).depth, 0);
  ASSERT_EQ(tree.node(child_node).depth, 1);
  ASSERT_EQ(tree.node(grandchild_node1).depth, 2);
  ASSERT_EQ(tree.node(grandchild_node2).depth, 2);
  ASSERT_EQ(tree.node(grand_grandchild_node).depth, 3);
  ASSERT_EQ(tree.node(grand_grand_grandchild_node).depth, 4);

  // Close the grand-grand-grandchild node, which should also reclaim the
  // grand-grandchild node and the grandchild_1 node
  tree.CloseNodeAndReclaimNodesUpToRootIfPossible(grand_grand_grandchild_node);
  ASSERT_EQ(tree.num_closed_nodes(), 3);

  // Now close the grandchild_2 node, which should reclaim both the child
  // and root node as well.
  tree.CloseNodeAndReclaimNodesUpToRootIfPossible(grandchild_node2);

  // Verify that tree is empty except for the root
  ASSERT_TRUE(tree.TreeIsEmpty());
  ASSERT_EQ(tree.num_closed_nodes(), 6);
}

TEST(MipTreeTest, SetLpRelaxationData) {
  MipTree tree;
  NodeIndex new_node =
      tree.AddNodeByBranchingFromParent(kRootNode, ColIndex(1), true, 1.5);

  double lp_value = -123.45;
  tree.SetLpRelaxationDataInNode(new_node, lp_value);
  EXPECT_EQ(tree.node(new_node).lp_objective_value, lp_value);
}

TEST(MipTreeTest, ImpliedBoundsHandling) {
  MipTree tree;
  NodeIndex new_node =
      tree.AddNodeByBranchingFromParent(kRootNode, ColIndex(0), true, 1.5);

  DenseRow root_lbs;
  root_lbs.push_back((ColIndex(0), -kInf));
  DenseRow root_ubs;
  root_ubs.push_back((ColIndex(0), kInf));

  DenseRow retrieved_lbs = tree.RetrieveLowerBounds(new_node, root_lbs);
  DenseRow retrieved_ubs = tree.RetrieveUpperBounds(new_node, root_ubs);

  EXPECT_EQ(retrieved_lbs.at(ColIndex(0)), -kInf);
  EXPECT_EQ(retrieved_ubs.at(ColIndex(0)), 1.0);

  std::vector<ColAndValue> lower_bounds = {{ColIndex(0), 4.0}};
  std::vector<ColAndValue> upper_bounds = {{ColIndex(0), 5.0}};
  tree.SetImpliedVariableBoundsInNode(new_node, lower_bounds, upper_bounds);

  DenseRow retrieved_lbs2 = tree.RetrieveLowerBounds(new_node, root_lbs);
  DenseRow retrieved_ubs2 = tree.RetrieveUpperBounds(new_node, root_ubs);

  EXPECT_EQ(retrieved_lbs2.at(ColIndex(0)), 4.0);
  EXPECT_EQ(retrieved_ubs2.at(ColIndex(0)), 1.0);
}

}  // namespace minimip