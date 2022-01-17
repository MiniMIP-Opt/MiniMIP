#ifndef SRC_LP_INTERFACE_SPARSE_TYPES_H_
#define SRC_LP_INTERFACE_SPARSE_TYPES_H_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>
#include <cassert>

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

// Abstract class for a sparse matrix representation
class AbstractSparseMatrix {

public:
  const int num_rows_;
  const int num_cols_;

  virtual ~AbstractSparseMatrix() {};
  AbstractSparseMatrix(int num_rows, int num_cols) : num_rows_(num_rows), num_cols_(num_cols) {}
  
  // Number of structural nonzeros in the matrix
  virtual size_t num_nonzeros() const = 0;

  // Accesses element (row_idx, col_idx) in the matrix.
  virtual double at(int row_idx, int col_idx) const = 0;

  // return a copy of all entries in the form (row_idx, col_idx, value)
  virtual std::vector<std::tuple<int, int, double>> all_entries() const = 0;

  virtual AbstractSparseMatrix& insert(int row_idx, int col_idx, double val) = 0;
};


// Sparse matrix in the Compressed Sparse Column (CSC) format
// Column j is kept at indices column_indices[j] to column_indices[j+1]-1
// We refer to "nonzeros" as the element stored in the matrix, even though some of them can have a value of 0.
// Common synonyms are "entries" and "structural nonzeros".
// row_indices contains the row index for each of the nonzeros.
// values contains the value for each of the nonzeros.
// Iterating over columns is efficient and can produce views without copies.
class SparseCompressedMatrix : public AbstractSparseMatrix {

public:
  // these fields remain public because they
  // will be necessary for some operations.
  std::vector<int> column_indices;
  std::vector<int> row_indices;
  std::vector<double> values;

  SparseCompressedMatrix(int num_rows, int num_cols) :
    AbstractSparseMatrix(num_rows, num_cols),
    column_indices(std::vector<int>(num_cols + 1)),
    row_indices({}), values({})
    {};
  
  SparseCompressedMatrix(int num_rows, int num_cols, std::vector<int> column_indices, std::vector<int> row_indices, std::vector<double> values) :
    AbstractSparseMatrix(num_rows, num_cols),
    column_indices(column_indices),
    row_indices(row_indices), values(values)
    {};
  
  SparseCompressedMatrix(const SparseCompressedMatrix& that) : 
    AbstractSparseMatrix(that.num_rows_, that.num_cols_),
    column_indices(that.column_indices),
    row_indices(that.row_indices), values(that.values)
    {};

  size_t num_nonzeros() const override {
    return values.size();
  }

  double at(int row_idx, int col_idx) const override {
    assert(row_idx >= 0);
    assert(col_idx >= 0);
    assert(row_idx < num_rows_);
    assert(col_idx < num_cols_);

    auto col_start = this->column_indices[col_idx];
    auto col_end = this->column_indices[col_idx + 1] - 1;

    if (col_start > col_end) {
      return 0.0;
    }
    int search_idx;
    for (search_idx = col_start; search_idx <= col_end; ++search_idx) {
      if (row_indices[search_idx] == row_idx) {
        break;
      }
      if (search_idx == col_end) {
        ++search_idx;
      }
    }
    if (search_idx > col_end || row_indices[search_idx] != row_idx) {
      return 0.0;
    }
    return values[search_idx];
  };

  // Inserts a value in the matrix.
  AbstractSparseMatrix& insert(int row_idx, int col_idx, double val) override {
    assert(row_idx >= 0);
    assert(col_idx >= 0);
    assert(row_idx < num_rows_);
    assert(col_idx < num_cols_);

    auto col_start = this->column_indices[col_idx];
    auto col_end = this->column_indices[col_idx + 1] - 1;

    int search_idx;
    if (row_indices.size() > 0 && row_idx < row_indices[col_start]) {
      // place new index before the current col_start
      search_idx = col_start;
    } else {
      for (search_idx = col_start; search_idx <= col_end; ++search_idx) {
        if (row_indices[search_idx] == row_idx) {
          break;
        }
        if (search_idx == col_end) {
          ++search_idx;
        }
      }
    }

    // column j exists and contains entry (i,j), replace it
    if (search_idx <= col_end && row_indices[search_idx] == row_idx) {
      values[search_idx] = val;
      return *this;
    }
    // if zero, no need to create new entry
    if (val == 0.0) {
      return *this;
    }
    if (search_idx >= static_cast<int>(row_indices.size())) {
      row_indices.push_back(row_idx);
      values.push_back(val);
    } else {
      row_indices.insert(row_indices.begin() + search_idx, row_idx);
      values.insert(values.begin() + search_idx, val);
    }
    // offset all column starting indices after new element
    for (auto col_offset_idx = col_idx + 1; col_offset_idx <= num_cols_; ++col_offset_idx) {
      ++column_indices[col_offset_idx];
    }
    return *this;
  };

  // TODO add a method returning an iterator over columns?

  // Vector of columns of the matrices.
  // These columns are views over the actual data of the sparse matrix.
  // Mutating the elements of the SparseVectorView will mutate the matrix itself.
  std::vector<SparseViewVector> column_views() {
    std::vector<SparseViewVector> cols;
    cols.reserve(this->num_cols_);
    for (auto col_idx = 0; col_idx < num_cols_; ++col_idx) {
      auto col_begin = this->column_indices[col_idx];
      auto col_end = this->column_indices[col_idx + 1] - 1;
      if (col_end < col_begin) {
        cols.push_back(SparseViewVector({}, {}, 0));
      } else {
        auto nnonzeros = col_end - col_begin + 1;
        cols.push_back(
          SparseViewVector(this->row_indices.data() + col_begin, this->values.data() + col_begin, nnonzeros)
        );
      }
    }
    return cols;
  };

  // return a copy of all entries in the form (row_idx, col_idx, value)
  std::vector<std::tuple<int, int, double>> all_entries() const override {
    std::vector<std::tuple<int, int, double>> entries;
    entries.reserve(this->num_nonzeros());
    for (auto col_idx = 0; col_idx < num_cols_; ++col_idx) {
      auto col_begin = this->column_indices[col_idx];
      auto col_end = this->column_indices[col_idx + 1] - 1;
      if (col_end >= col_begin) {
        auto nnonzeros = col_end - col_begin + 1;
        for (auto elem_idx = 0; elem_idx < nnonzeros; ++elem_idx) {
          entries.push_back(std::tuple<int, int, double>(this->row_indices[col_begin + elem_idx], col_idx, this->values[col_begin + elem_idx]));
        }
      }
    }
    return entries;
  };
};

// The transpose of a SparseCompressedMatrix
// Maintains a pointer to the data of the transposed matrix.
// And whether the matrix allocation was done during the construction.
// If this is the case, the destructor cleans up the transpose. 
class RowSparseMatrix : public AbstractSparseMatrix {
private:
  SparseCompressedMatrix* transposed_;
  bool is_allocated;
public:
  // Allocating constructor creating the underlying transposed matrix.
  RowSparseMatrix(int num_row, int num_col) :
    AbstractSparseMatrix(num_row, num_col),
    transposed_(new SparseCompressedMatrix(num_row, num_col)),
    is_allocated(true) {};
  
  // Zero-allocation constructor, constructs a RowSparseMatrix from the
  // column matrix to transpose
  RowSparseMatrix(SparseCompressedMatrix* transposed) :
    AbstractSparseMatrix(transposed->num_cols_, transposed->num_rows_),
    transposed_(transposed), is_allocated(false) {};

  ~RowSparseMatrix() {
    if (this->is_allocated) {
      delete this->transposed_;
    }
  }
  SparseCompressedMatrix const & transposed_view() const {
    return *(this->transposed_);
  };

  size_t num_nonzeros() const override {
    return transposed_->num_nonzeros();
  };

  // Accesses element (row_idx, col_idx) in the matrix.
  double at(int row_idx, int col_idx) const override {
    return transposed_->at(col_idx, row_idx);
  };

  AbstractSparseMatrix& insert(int row_idx, int col_idx, double val) override {
    assert(row_idx >= 0);
    assert(col_idx >= 0);
    assert(row_idx < num_rows_);
    assert(col_idx < num_cols_);
    this->transposed_->insert(col_idx, row_idx, val);
    return *this;
  };

  // return a copy of all entries in the form (row_idx, col_idx, value)
  std::vector<std::tuple<int, int, double>> all_entries() const override {
    auto entries = transposed_->all_entries();
    for (size_t idx; idx < entries.size(); ++idx) {
      entries[idx] = std::tuple<int, int, double>(
        std::get<1>(entries[idx]),
        std::get<0>(entries[idx]),
        std::get<2>(entries[idx])
      );
    }
    return entries;
  };

  // TODO add a method returning an iterator over rows?

  // Vector of rows of the matrix.
  // See col_views on CompressedSparseMatrix
  std::vector<SparseViewVector> row_views() {
    return this->transposed_->column_views();
  };

};

}  // namespace minimip
#endif  // SRC_LP_INTERFACE_SPARSE_TYPES_H_
