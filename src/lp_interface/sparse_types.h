#ifndef SRC_LP_INTERFACE_SPARSE_TYPES_H_
#define SRC_LP_INTERFACE_SPARSE_TYPES_H_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <numeric>
#include <string>
#include <vector>

namespace minimip {

class SparseVector;

// Represents one entry in a sparse one-dimensional structure.
struct SparseEntry {
  int index;
  double value;
};

// Abstract sparse vector represented by its structural nonzeros indices and
// values. Indices are maintained in increasing order.
class AbstractSparseVector {
 public:
  virtual size_t NumNonZeros() const = 0;
  virtual bool empty() const { return NumNonZeros() == 0; };
  virtual const int* IndicesData() const   = 0;
  virtual const double* ValuesData() const = 0;
  // copies the current abstract sparse vector to a new SparseVector
  virtual SparseVector CopyToOwned() const = 0;
  virtual double ValueAt(int idx) const {
    const int* indices   = this->IndicesData();
    const double* values = this->ValuesData();
    for (size_t i = 0; i < this->NumNonZeros(); ++i) {
      if (indices[i] == idx) return values[i];
    }
    return 0.0;
  }
  // Returns the entry at the idx-th nonzero of the sparse vector.
  virtual SparseEntry EntryAt(size_t idx) {
    assert(idx >= 0);
    assert(idx < NumNonZeros());
    return SparseEntry{IndicesData()[idx], ValuesData()[idx]};
  };
};

// Sparse vector type
// Note that it does not specify its true size
// (which could be greater than max(indices)).
class SparseVector : public AbstractSparseVector {
 public:
  SparseVector(std::vector<int> indices, std::vector<double> values) {
#ifndef NDEBUG
    assert(indices.size() == values.size());
    if (!indices.empty()) {
      for (size_t idx = 0; idx < indices.size() - 1; ++idx) {
        assert(indices[idx] < indices[idx + 1]);
      }
    }
#endif
    this->indices_ = indices;
    this->values_  = values;
  };
  SparseVector() : indices_({}), values_({}){};

  // Number of structural nonzeros of the sparse vector
  // Indices are assumed to be unique and sorted.
  size_t NumNonZeros() const override { return indices_.size(); }

  bool empty() const override { return indices_.empty(); }

  const int* IndicesData() const override { return this->indices_.data(); }

  const double* ValuesData() const override { return this->values_.data(); }

  SparseVector CopyToOwned() const override {
    return SparseVector(this->indices_, this->values_);
  }

  // insert an index-value pair while maintaining sorted indices.
  // if the index is already present, the corresponding value is replaced.
  // Inserting a 0.0 is equivalent to removing the corresponding entry if it
  // exists.
  SparseVector& InsertSorted(int index, double value) {
    // only avoid insertion if value is a true 0 (no tolerance)
    if (indices_.empty() || index > indices_[NumNonZeros() - 1]) {
      if (value != 0.0) {
        indices_.push_back(index);
        values_.push_back(value);
      }
      return *this;
    }
    size_t first_greater = 0;
    while (indices_[first_greater] < index) ++first_greater;

    if (indices_[first_greater] == index) {
      if (value == 0.0) {
        indices_.erase(indices_.begin() + first_greater);
        values_.erase(values_.begin() + first_greater);
      } else {
        values_[first_greater] = value;
      }
    } else {
      indices_.insert(indices_.begin() + first_greater, index);
      values_.insert(values_.begin() + first_greater, value);
    }
    return *this;
  }

  // Replaces the nonzero entry at position by (index, value)
  // permutes entries to maintain order
  SparseVector& ReplaceEntry(size_t position, int index, double value) {
    assert(position < NumNonZeros());
    this->indices_[position] = index;
    this->values_[position]  = value;
    Resort();
    return *this;
  }

  SparseVector& reserve(int capacity) {
    this->values_.reserve(capacity);
    this->indices_.reserve(capacity);
    return *this;
  }

 private:
  std::vector<int> indices_;
  std::vector<double> values_;

  // resorts the array when operations broke indices ordering
  void Resort() {
    // create permutation of indices
    std::vector<std::size_t> p(indices_.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(), std::less<int>());

    std::vector<bool> done(p.size());
    for (size_t i = 0; i < p.size(); ++i) {
      if (done[i]) {
        continue;
      }
      done[i]       = true;
      size_t prev_j = i;
      size_t j      = p[i];
      while (i != j) {
        std::swap(indices_[prev_j], indices_[j]);
        std::swap(values_[prev_j], values_[j]);
        done[j] = true;
        prev_j  = j;
        j       = p[j];
      }
    }
  }
};

// A sparse vector structure only holding views
// to the underlying indices and coefficients
class SparseViewVector : public AbstractSparseVector {
 public:
  SparseViewVector(int* indices, double* values, int num_nonzeros) {
    for (int idx = 0; idx < num_nonzeros - 1; ++idx) {
      assert(indices[idx] < indices[idx + 1]);
    }
    this->indices_      = indices;
    this->values_       = values;
    this->num_nonzeros_ = num_nonzeros;
  };

  size_t NumNonZeros() const override { return num_nonzeros_; }

  const int* IndicesData() const override { return indices_; }
  const double* ValuesData() const override { return values_; }
  SparseVector CopyToOwned() const override {
    SparseVector v(
        std::vector<int>(this->indices_, this->indices_ + num_nonzeros_),
        std::vector<double>(this->values_, this->values_ + num_nonzeros_));
    return v;
  }

 private:
  int* indices_;
  double* values_;
  int num_nonzeros_;
};

}  // namespace minimip
#endif  // SRC_LP_INTERFACE_SPARSE_TYPES_H_
