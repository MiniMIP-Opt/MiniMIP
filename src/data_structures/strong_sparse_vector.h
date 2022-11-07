// Copyright 2022 the MiniMIP Project
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
//
// -----------------------------------------------------------------------------
//
// Implementation of `SparseRow` and `SparseCol` by templatizing a lazy, strong,
// sparse vector containing doubles.
//
// "Sparse" means that (in the "clean" state) we only store non-zeros. To this
// end, we store pairs of (sparse index, value), which are sorted by sparse
// index, there are no duplicates, and value is a non-zero.
//
// "Lazy" means that the entries are actually appended in O(1), thus - in the
// interim - we may store store zeros, duplicated and/or unsorted entries.
// A client calls `CleanUpIfNeeded()` to sort and remove zeros and duplicates.
//
// "Strong" means we use strongly typed ints to distinguish between col, row and
// entry at compile time (i.e., provide strong type safety guarantees).
//
// Example:
//   StrongRow row({{1, 10.}, {5, 50.}});
//   row.AddEntry(RowIndex(2), 20.);
//   row.AddEntry(RowIndex(1), 0.);
//   row.AddEntry(RowIndex(5), 25.0);
//   row.CleanUpIfNeeded();
//   for (const auto& e : row.entries()) {
//     LOG(INFO) << e;
//   }
//
// This will print (in this order):
//   index=2, value=20.
//   index=5, value=25.

#ifndef SRC_DATA_STRUCTURES_STRONG_SPARSE_VECTOR_H_
#define SRC_DATA_STRUCTURES_STRONG_SPARSE_VECTOR_H_

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
using ColEntry = SparseEntry<RowIndex>;
using RowEntry = SparseEntry<ColIndex>;
using SparseRow = StrongSparseVectorOfDoubles<ColIndex>;
using SparseCol = StrongSparseVectorOfDoubles<RowIndex>;

template <typename SparseIndex>
constexpr SparseIndex kInvalidSparseIndex(-1);

[[maybe_unused]] constexpr RowIndex kInvalidRow = kInvalidSparseIndex<RowIndex>;
[[maybe_unused]] constexpr ColIndex kInvalidCol = kInvalidSparseIndex<ColIndex>;

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
  // The class is copyable.
  // TODO(lpawel): Re-consider this. Now, it's easy to use a copy in place of a
  // reference by mistake.
  StrongSparseVectorOfDoubles(const StrongSparseVectorOfDoubles&) = default;
  StrongSparseVectorOfDoubles<SparseIndex>& operator=(
      const StrongSparseVectorOfDoubles<SparseIndex>&) = default;

  // The class is (no-throw) movable.
  StrongSparseVectorOfDoubles(StrongSparseVectorOfDoubles&&) noexcept = default;
  StrongSparseVectorOfDoubles<SparseIndex>& operator=(
      StrongSparseVectorOfDoubles<SparseIndex>&&) noexcept = default;

  StrongSparseVectorOfDoubles() = default;

  explicit StrongSparseVectorOfDoubles(
      absl::StrongVector<EntryIndex, SparseEntry<SparseIndex>> entries)
      : entries_(std::move(entries)), may_need_cleaning_(true) {
    CleanUpIfNeeded();
    DCHECK(IsClean());  // This also verifies whether indices are >= 0.
  }

  const absl::StrongVector<EntryIndex, SparseEntry<SparseIndex>>& entries()
      const {
    DCHECK(!MayNeedCleaning())
        << "The vector is not clean! Call `CleanUpIfNeeded()`.";
    return entries_;
  }

  // Use with care! Use this to in-place modify values and indices of entries.
  // To add new entries use `AddEntry()`. Regardless of what modifications are
  // introduced, the vector will be marked for cleaning.
  absl::StrongVector<EntryIndex, SparseEntry<SparseIndex>>& mutable_entries() {
    DCHECK(!MayNeedCleaning())
        << "The vector is not clean! Call `CleanUpIfNeeded()`.";
    // To be on the safe side, we unconditionally mark the entries for clean up.
    // TODO(lpawel): Reconsider this, if it turns out inefficient.
    may_need_cleaning_ = true;
    return entries_;
  }

  std::vector<SparseIndex> indices() const {
    DCHECK(!MayNeedCleaning())
        << "The vector is not clean! Call `CleanUpIfNeeded()`.";
    std::vector<SparseIndex> indices;
    indices.reserve(entries_.size());
    for (const auto& entry : entries_) {
      indices.push_back(entry.index);
    }
    return indices;
  }

  std::vector<double> values() const {
    DCHECK(!MayNeedCleaning())
        << "The vector is not clean! Call `CleanUpIfNeeded()`.";
    std::vector<double> values;
    values.reserve(entries_.size());
    for (const auto& entry : entries_) {
      values.push_back(entry.value);
    }
    return values;
  }

  // Adds new entry and marks the vector for cleaning (if needed).
  void AddEntry(SparseIndex index, double value) {
    if (!may_need_cleaning_) {
      const SparseIndex last_index = entries_.empty()
                                         ? kInvalidSparseIndex<SparseIndex>
                                         : entries_.back().index;
      if (index <= last_index) {
        // The entries will be unsorted or duplicated and will require cleaning.
        may_need_cleaning_ = true;
      } else if (value == 0.0) {
        // The new index doesn't have a value, so no need to set it to 0.
        return;
      }
    }
    // STL resizes the underlying vector exponentially on growth beyond current
    // capacity (to minimize the number of memory reallocations and moves).
    entries_.emplace_back(index, value);
  }

  // Gets the value for a specific value of sparse index. This assumes the
  // vector is cleaned up and runs in O(log(entries().size()).
  double value(SparseIndex index) const {
    DCHECK(!may_need_cleaning_)
        << "The vector is not clean! Call `CleanUpIfNeeded()`.";
    const auto& ge_entry = std::lower_bound(
        entries_.begin(), entries_.end(), index,
        [](const SparseEntry<SparseIndex>& a, const SparseIndex i) -> bool {
          return a.index < i;
        });
    return (ge_entry == entries_.end() || ge_entry->index != index)
               ? 0.0
               : ge_entry->value;
  }

  double operator[](SparseIndex index) const { return value(index); }

  // Removes all entries, but does not release memory of the underlying storage.
  void Clear() {
    entries_.clear();
    may_need_cleaning_ = false;
  }

  // Removes duplicates (the last entry takes precedence), removes zero entries,
  // and sorts entries. Runs in O(n log(n)), where n = entries().size()).
  void CleanUpIfNeeded() {
    if (!MayNeedCleaning()) return;
    std::stable_sort(entries_.begin(), entries_.end(),
                     [](const SparseEntry<SparseIndex>& lhs,
                        const SparseEntry<SparseIndex>& rhs) -> bool {
                       return lhs.index < rhs.index;
                     });

    EntryIndex new_pos(0);
    for (EntryIndex pos(0); pos < EntryIndex(entries_.size()); ++pos) {
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
    may_need_cleaning_ = false;
  }

  bool MayNeedCleaning() const {
    // Note, `!IsClean() => may_need_cleaning_`. In particular,
    // `may_need_cleaning_` may be true even if `IsClean()` is false (because
    // we unconditionally mark the entries for cleaning when returning
    // `mutable_entries()`).
    DCHECK(IsClean() || may_need_cleaning_);
    return may_need_cleaning_;
  }

  SparseIndex last_index() const {
    DCHECK(!may_need_cleaning_);
    return entries_.empty() ? kInvalidSparseIndex<SparseIndex>
                            : entries_.back().index;
  }

  SparseIndex first_index() const {
    DCHECK(!may_need_cleaning_);
    return entries_.empty() ? kInvalidSparseIndex<SparseIndex>
                            : entries_.front().index;
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

  StrongSparseVectorOfDoubles operator-() const {
    DCHECK(!may_need_cleaning_);
    StrongSparseVectorOfDoubles copy = *this;
    for (auto& [unused, value] : copy.entries_) value = -value;
    return copy;
  }

  StrongSparseVectorOfDoubles& operator*=(const double scalar) {
    DCHECK(!may_need_cleaning_);
    if (scalar == 0.0) {
      entries_.clear();
      return *this;
    }
    for (auto& [unused, value] : entries_) value *= scalar;
    return *this;
  }

  StrongSparseVectorOfDoubles operator*(const double scalar) const {
    StrongSparseVectorOfDoubles copy = *this;
    copy *= scalar;
    return copy;
  }

  StrongSparseVectorOfDoubles& operator/=(const double scalar) {
    DCHECK(!may_need_cleaning_);
    DCHECK_NE(scalar, 0.0);
    for (auto& [unused, value] : entries_) value /= scalar;
    return *this;
  }

  StrongSparseVectorOfDoubles operator/(const double scalar) const {
    StrongSparseVectorOfDoubles copy = *this;
    copy /= scalar;
    return copy;
  }

  StrongSparseVectorOfDoubles& operator+=(
      const StrongSparseVectorOfDoubles& other) {
    return AddMultipleOfVector(1.0, other);
  }

  StrongSparseVectorOfDoubles operator+(
      const StrongSparseVectorOfDoubles& other) const {
    StrongSparseVectorOfDoubles copy = *this;
    copy += other;
    return copy;
  }

  StrongSparseVectorOfDoubles& operator-=(
      const StrongSparseVectorOfDoubles& other) {
    return AddMultipleOfVector(-1.0, other);
  }

  StrongSparseVectorOfDoubles operator-(
      const StrongSparseVectorOfDoubles& other) const {
    StrongSparseVectorOfDoubles copy = *this;
    copy -= other;
    return copy;
  }

  // Performs *this = *this + multiple * vector. Time complexity is linear in
  // the total number of entries in both vectors.
  StrongSparseVectorOfDoubles& AddMultipleOfVector(
      double factor, const StrongSparseVectorOfDoubles& other) {
    DCHECK(!may_need_cleaning_);
    DCHECK(!other.may_need_cleaning_);

    // The algorithm we use is:
    // 1. Merge the added entries into our current entry vector. O(n), since
    //    both vectors are sorted.
    // 2. Iterate over entries and for each duplicated index, add the values
    //    into a single entry and set the duplicate values to 0. O(n)
    // 3. Remove all entries that have value 0. O(n)

    const int original_size = entries_.size();
    entries_.resize(entries_.size() + other.entries_.size());
    // Note: resize may invalidate iterators, so we must create it here.
    const auto new_entries_it = entries_.begin() + original_size;
    std::transform(other.entries_.begin(), other.entries_.end(), new_entries_it,
                   [factor](const SparseEntry<SparseIndex>& e) {
                     return SparseEntry<SparseIndex>(e.index, e.value * factor);
                   });
    std::inplace_merge(
        entries_.begin(), new_entries_it, entries_.end(),
        [](const SparseEntry<SparseIndex>& e1,
           const SparseEntry<SparseIndex>& e2) { return e1.index < e2.index; });

    SparseIndex current_unique_index(kInvalidSparseIndex<SparseIndex>);
    double* value_at_current_unique_index = nullptr;
    for (auto& [index, value] : entries_) {
      if (index != current_unique_index) {
        // Note: this also happens on the first iteration, so the pointer is
        // initialized correctly.
        current_unique_index = index;
        value_at_current_unique_index = &value;
      } else {
        *value_at_current_unique_index += value;
        value = 0.0;
      }
    }

    entries_.erase(std::remove_if(entries_.begin(), entries_.end(),
                                  [](const SparseEntry<SparseIndex>& e) {
                                    return e.value == 0.0;
                                  }),
                   entries_.end());
    return *this;
  }

  template <typename OtherSparseIndex>
  double DotProduct(
      const StrongSparseVectorOfDoubles<OtherSparseIndex>& other) const {
    DCHECK(!may_need_cleaning_);
    DCHECK(!other.may_need_cleaning_);
    double product = 0.0;
    auto it1 = entries_.begin();
    auto it2 = other.entries_.begin();
    while (it1 != entries_.end() && it2 != other.entries_.end()) {
      if (it1->index.value() == it2->index.value()) {
        product += it1->value * it2->value;
        ++it1;
        ++it2;
      } else if (it1->index.value() < it2->index.value()) {
        ++it1;
      } else {
        ++it2;
      }
    }
    return product;
  }

  // We allow all instantiations to access each other's private members.
  template <typename OtherSparseIndex>
  friend class StrongSparseVectorOfDoubles;

 private:
  // With this, the index and the corresponding value are kept next to each
  // other in memory. Thus, likely they end up in the same page while
  // accessing (which is good!). However, the compiler will most likely pad
  // this to 16 bytes (see, https://godbolt.org/z/jnvfv3rdn). So, we're
  // over-using memory by 25%.
  // TODO(lpawel): Investigate if this becomes an issue.
  absl::StrongVector<EntryIndex, SparseEntry<SparseIndex>> entries_;
  bool may_need_cleaning_ = false;
};

// Allow the scalar to be the first operand.
template <typename SparseIndex>
StrongSparseVectorOfDoubles<SparseIndex> operator*(
    double scalar, const StrongSparseVectorOfDoubles<SparseIndex>& vector) {
  return vector * scalar;
}

// Note, we need no-throw movable so that `std::vector<SparseRow>::resize()`
// uses move semantics.
static_assert(std::is_move_assignable<SparseRow>(),
              "SparseRow is not move assignable.");
static_assert(std::is_nothrow_move_assignable<SparseRow>(),
              "SparseRow is not nothrow move assignable.");

}  // namespace minimip
#endif  // SRC_DATA_STRUCTURES_STRONG_SPARSE_VECTOR_H_
