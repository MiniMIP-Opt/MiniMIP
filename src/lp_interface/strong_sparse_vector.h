// Implementation of `SparseRow` and `SparseCol` by templatizing a lazy, strong,
// sparse vector containing doubles.
//
// "Sparse" means that (in the "clean" state) we only store non-zeros. To this
// end, we store pairs of (sparse index, value), which are sorted by sparse
// index, there are no duplicates, and value is a non-zero.
//
// "Lazy" means that the entries are actually appended in O(1), thus - in the
// interim - we may store store zeros, duplicated and/or unsorted entries.
// A client calls `CleanUp()` to sort and remove zeros and duplicates.
//
// "Strong" means we use strongly typed ints to distinguish between col, row and
// entry at compile time (i.e., provide strong type safety guarantees).
//
// Example:
//   StrongRow row({{1, 10.}, {5, 50.}});
//   row.AddEntry(RowIndex(2), 20.);
//   row.AddEntry(RowIndex(1), 0.);
//   row.AddEntry(RowIndex(5), 25.0);
//   row.CleanUp();
//   for (const auto& e : row.entries()) {
//     LOG(INFO) << e;
//   }
//
// This will print (in this order):
//   index=2, value=20.
//   index=5, value=25.

#ifndef SRC_LP_INTERFACE_STRONG_SPARSE_VECTOR_H_
#define SRC_LP_INTERFACE_STRONG_SPARSE_VECTOR_H_

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "ortools/base/logging.h"
#include "ortools/base/strong_vector.h"
#include "ortools/util/strong_integers.h"

namespace minimip {

// These are strong types for columns (aka variables) and rows (aka
// constraints) used everywhere in MiniMip. The individual entries of a
// sparse vector are indexed with `EntryIndex`.
DEFINE_STRONG_INDEX_TYPE(ColIndex);
DEFINE_STRONG_INDEX_TYPE(RowIndex);
DEFINE_STRONG_INDEX_TYPE(EntryIndex);

// Forward declarations.
template <typename SparseIndex>
struct SparseEntry;
template <typename SparseIndex>
class StrongSparseVectorOfDoubles;

// Handy shortcuts, use those in code.
using ColEntry  = SparseEntry<RowIndex>;
using RowEntry  = SparseEntry<ColIndex>;
using SparseRow = StrongSparseVectorOfDoubles<ColIndex>;
using SparseCol = StrongSparseVectorOfDoubles<RowIndex>;

// This is needed for CHECK_EQ(), EXPECT_EQ(), and macros from abseil.
template <typename SparseIndex>
bool operator==(const SparseEntry<SparseIndex>& lhs,
                const SparseEntry<SparseIndex>& rhs) {
  return lhs.index == rhs.index && lhs.value == rhs.value;
}

// This is needed for nice printing in unit tests and via LOG().
template <typename SparseIndex>
std::ostream& operator<<(std::ostream& out, const SparseEntry<SparseIndex>& e) {
  out << "{.index = " << e.index.value() << ", .value = " << e.value << "}";
  return out;
}

// Represents one entry in a sparse one-dimensional vector.
template <typename SparseIndex>
struct SparseEntry {
  // This is needed for `resize()` in STL containers storing SparseEntries.
  SparseEntry() : index(-1), value(0.0) {}

  // This is needed for initializer lists, e.g., `SparseRow r = {{1, 2.0}}`.
  SparseEntry(SparseIndex index, double value) : index(index), value(value) {}

  SparseIndex index;
  double value;
};

template <typename SparseIndex>
class StrongSparseVectorOfDoubles {
 public:
  StrongSparseVectorOfDoubles() = default;

  explicit StrongSparseVectorOfDoubles(
      absl::StrongVector<EntryIndex, SparseEntry<SparseIndex>> entries)
      : entries_(std::move(entries)), needs_cleaning_(true) {
    CleanUp();
    DCHECK(IsClean());  // This also verifies whether indices are >= 0.
  }

  // Use the default copy constructor.
  StrongSparseVectorOfDoubles(const StrongSparseVectorOfDoubles&) = default;

  // Use a default move constructor, but mark it noexcept. The constructor
  // is defined outside of the body, for details on this technique see
  // https://akrzemi1.wordpress.com/2015/09/11/declaring-the-move-constructor/.
  // Thanks to this, C++ will use move constructor on reallocation.
  StrongSparseVectorOfDoubles(StrongSparseVectorOfDoubles&&) noexcept;

  const absl::StrongVector<EntryIndex, SparseEntry<SparseIndex>>& entries()
      const {
    DCHECK(!needs_cleaning_) << "The vector is not clean! Call `CleanUp()`.";
    return entries_;
  }

  // Use with care! Use this to in-place modify values and indices of entries.
  // To add new entries use `AddEntry()`. Regardless of what modifications are
  // introduced, the vector will be marked for cleaning.
  absl::StrongVector<EntryIndex, SparseEntry<SparseIndex>>& mutable_entries() {
    DCHECK(!needs_cleaning_) << "The vector is not clean! Call `CleanUp()`.";
    // To be on the safe side, we unconditionally mark the entries for clean up.
    // TODO(lpawel): Reconsider this, if it turns out inefficient.
    needs_cleaning_ = true;
    return entries_;
  }

  // Adds new entry and marks the vector for cleaning (if needed).
  void AddEntry(SparseIndex index, double value) {
    if (!needs_cleaning_ && !entries_.empty()) {
      // The entries will be unsorted or duplicated and will require cleaning.
      if (entries_.back().index >= index) {
        needs_cleaning_ = true;
      }
    }
    // STL resizes the underlying vector exponentially on growth beyond current
    // capacity (to minimize the number of memory reallocations and moves).
    entries_.emplace_back(index, value);
  }

  // Gets the value for a specific value of sparse index. This assumes the
  // vector is cleaned up and runs in O(log(entries().size()).
  double value(SparseIndex index) const {
    DCHECK(!needs_cleaning_) << "The vector is not clean! Call `CleanUp()`.";
    const auto& ge_entry = std::lower_bound(
        entries_.begin(), entries_.end(), index,
        [](const SparseEntry<SparseIndex>& a, const SparseIndex i) -> bool {
          return a.index < i;
        });
    return (ge_entry == entries_.end() || ge_entry->index != index)
               ? 0.0
               : ge_entry->value;
  }

  // Removes all entries, but does not release memory of the underlying storage.
  void Clear() {
    entries_.clear();
    needs_cleaning_ = false;
  }

  // Removes duplicates (the last entry takes precedence), removes zero entries,
  // and sorts entries. Runs in O(n log(n)), where n = entries().size()).
  void CleanUp() {
    if (!needs_cleaning_) return;
    std::stable_sort(entries_.begin(), entries_.end(),
                     [](const SparseEntry<SparseIndex>& lhs,
                        const SparseEntry<SparseIndex>& rhs) -> bool {
                       return lhs.index < rhs.index;
                     });

    EntryIndex new_pos = 0;
    for (EntryIndex pos = 0; pos < EntryIndex(entries_.size()); ++pos) {
      if (pos < entries_.size() - 1 &&
          entries_[pos].index == entries_[pos + 1].index) {
        continue;
      }
      if (entries_[pos].value == 0.0) {
        continue;
      }
      entries_[new_pos++] = entries_[pos];
    }
    entries_.resize(new_pos.value());
    needs_cleaning_ = false;
  }

  // Verifies whether `entries()` are clean (i.e., sorted, no duplicates, no
  // zeros).
  bool IsClean() const {
    // The indices must be in strictly increasing order, i.e., no duplicates
    // (hence "<=" in the comparison).
    if (!std::is_sorted(entries_.begin(), entries_.end(),
                        [](const auto& lhs, const auto& rhs) {
                          return lhs.index <= rhs.index;
                        })) {
      return false;
    }

    // All entries must be non-zeros and the indices must be non-negative.
    return std::all_of(entries_.begin(), entries_.end(), [](const auto& e) {
      return e.value != 0.0 && e.index >= SparseIndex(0);
    });
  }

 private:
  // With this, the index and the corresponding value are kept next to each
  // other in memory. Thus, likely they end up in the same page while accessing
  // (which is good!). However, the compiler will most likely pad this to 16
  // bytes (see, https://godbolt.org/z/jnvfv3rdn). So, we're over-using memory
  // by 25%.
  // TODO(lpawel): Investigate if this becomes an issue.
  absl::StrongVector<EntryIndex, SparseEntry<SparseIndex>> entries_;
  bool needs_cleaning_ = false;
};

template <typename SparseIndex>
inline StrongSparseVectorOfDoubles<SparseIndex>::StrongSparseVectorOfDoubles(
    StrongSparseVectorOfDoubles&&) noexcept = default;

}  // namespace minimip
#endif  // SRC_LP_INTERFACE_STRONG_SPARSE_VECTOR_H_
