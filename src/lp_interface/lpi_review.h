#ifndef MINIMIP_SRC_LP_INTERFACE_LPI_H
#define MINIMIP_SRC_LP_INTERFACE_LPI_H

// NOTE: This doesn't compile (yet, e.g., there are missing includes for abseil),
// and it's here mainly for the purpuse of the review of LP interface.

#include <string>
#include <vector>

namespace minimip {

// LPInterface abstract class. This is used my MiniMIP to communicate with an
// underlying LP solver.
class LpInterface {
 public:
  virtual ~LPInterface() = default;
  
  // ==========================================================================
  // LP model setters.
  // ==========================================================================

  struct LpColumnData {
    // This column has constraint coefficients equal to `coeffs[i]` at
    // constraint `rows[i]`. Both vectors must be of equal length. All row
    // indices must be numbers between 0 and `num_rows()`. All coeffs must be
    // finite, non-zero numbers.
    // TODO(lpawel): Think whether we should accept coeffs equal to 0.0.
    std::vector<int> rows;
    std::vector<double> coeffs;

    // Lower and upper bound of the column. `lower_bound` and upper_bound may
    // be equal to -infinity and +infinity, respectively, to indicate that no
    // corresponding bound is imposed.
    double lower_bound;
    double upper_bound;

    // Objective coefficient of the column. Must be finite.
    double objective_coefficient;

    // Optional name of the column (i.e., may be empty).
    std::string name;
  };

  struct LpRowData {
    // This row has constraint coefficients equal to `coeffs[i]` for variable
    // `cols[i]`. Both vectors must be of equal length. All row indices must be
    // numbers between 0 and `num_rows()`. All coeffs must be finite, non-zero
    // numbers.
    // TODO(lpawel): Think whether we should accept coeffs equal to 0.0.
    std::vector<int> cols;
    std::vector<double> coeffs;

    // Lower and upper bound of the row. `lower_bound` and `upper_bound` may
    // be equal to -infinity and +infinity, respectively, to indicate that no
    // corresponding bound is imposed.
    double lower_bound;
    double upper_bound;

    // Optional name of the row (i.e., may be empty).
    std::string name;
  }

  // All row indices inside elements from `cols` must be between 0 and
  // `rows.size() - 1`. All col indices inside elements from `rows` must be
  // between 0 and `cols.size() - 1`. The coefficient matrix is populated from
  // `cols` (i.e., for all SparseLpRow from `rows`, their `cols` and `coeffs`
  // must be empty).
  // TODO(lpawel): Consider allowing to populate the LP from rows too.
  virtual absl::Status PopulateFromSparseColumns(
      bool is_minimization,
      const std::vector<LpColumnData>& cols,
      const std::vector<LpColumnData>& rows) = 0;

  // Appends a single column to the LP from the right side of the coefficient
  // matrix.
  virtual absl::Status AddColumn(const LpColumnData& col) = 0;

  // Deletes the column at index `col`, which must be a number between 0 and
  // `num_columns()`.
  virtual absl::Status DeleteColumn(int col) = 0;

  // Appends a single row to the LP from the bottom side of the coefficient
  // matrix.
  virtual absl::Status AddRow(const LpRowData& row) = 0;

  // Deletes the row at index `row`, which must be a number between 0 and
  // `num_rows()`.
  virtual absl::Status DeleteRow(int row) = 0;

  // Sets new lower and upper bounds for a column and row.
  // NOTE: We do not provide a bulk version, because both Glop and SoPLEX
  // would need to iterate over the updated columns anyway.
  virtual absl::Status SetColumnBounds(int col, double lower_bound,
                                       double upper_bound) = 0;
  virtual absl::Status SetRowBounds(int row, double lower_bound,
                                    double upper_bound) = 0;

  // Sets the objective sense (i.e., direction of optimization).
  // NOTE: For the purpose of LP solving, MiniMIP assumes the objective
  // offset is 0.
  virtual void SetMaximization(bool is_maximization) = 0;
    
  // Sets new objective coefficient for a column.
  // NOTE: We do not provide a bulk version, because both Glop and SoPLEX
  // would need to iterate over the updated columns anyway.  
  virtual absl::Status SetObjectiveCoefficent(int col,
                                              double objective_coefficient) = 0;

  // Removes all cols and rows (and associated coeffiecients, bounds, and
  // names).
  virtual void ClearLpData() = 0;

  // Clears all information (e.g., basis) except for cols and rows (and the
  // associated coefficients, bounds, and names).
  virtual void ClearLpState() = 0;                                              

  // ==========================================================================
  // LP model getters.
  // ==========================================================================
  
  // Simple getters.
  // TODO(lpawel): Style guide allows for simple getters to be in the lower
  // case format. For consistency though, always use camel format?
  virtual int num_rows() const = 0;
  virtual int num_columns() const = 0;
  virtual int num_non_zeros() const = 0;
  virtual bool is_maximization() const = 0;

  // Gets the objective coefficient for `col` (which must be between 0 and
  // `num_cols()-1`).
  virtual double GetObjectiveCoefficient(int col) const = 0;

  // Gets the sparse column at `col` index (which must be between 0 and
  // `num_cols()`).
  // NOTE: We do not provide a bulk version for getting a row / column, because
  // both Glop and SoPLEX would need to iterate over the rows and cols anyway.    
  virtual LpColumnData GetColumnData(int col) const = 0;

  // Gets the sparse row at `row` index (which must be between 0 and
  // `num_rows()-1`).
  // NOTE: We do not provide a bulk version for getting a row / column, because
  // both Glop and SoPLEX would need to iterate over the rows and cols anyway.      
  virtual LpRowData GetRowData(int row) const = 0;

  // Gets the lower and upper bound for `col` (which must be between 0 and
  // `num_cols()-1`).
  struct ColBounds {
    double lower_bound;
    double upper_bound;
  }
  virtual ColBounds GetColBounds(int col) const = 0;

  // Gets the lower and upper bound for `row` (which must be between 0 and
  // `num_rows()-1`).
  struct RowBounds {
    double lower_bound;
    double upper_bound;
  };
  virtual RowBounds GetRowBounds(int row) const = 0;    

  // Gets a coefficient from the constaint matrix for `col` and `row` (`col` and
  // `row` must be between 0 and `num_cols()-1` and `num_rows()-1`,
  // respectively).
  virtual double GetMatrixCoefficient(int col, int row) const = 0;

  // ==========================================================================
  // Solving methods.
  // ==========================================================================

  // Solve the current LP model with primal or dual simplex. Note, these
  // methods call the underlying LP solver "as is" (i.e., they do not clear or
  // set its state, and the solve may continue from the currently set basis, if
  // the underlying solver supports that). It is the responsability of the
  // caller to reset the LP solver with a specific basis.
  virtual absl::Status SolveLpWithPrimalSimplex() = 0;
  virtual absl::Status SolveLpWithDualSimplex() = 0;

  // These two functions should always be called before and after the underlying
  // LP solver is used to perform strong branching solves.
  virtual absl::Status StartStrongBranchingingMode() = 0;
  virtual absl::Status EndStrongBranchingingMode() = 0;
  
  // Executes strong branching probing of the up and down branch for the `col`
  // around `primal_value`. Note, `primal_value` may be integer in which case,
  // "integer strong branching" is executed (as of 2021/09/29, this is supported
  // by SoPlex, but not by Glop).
  struct StrongBranchResult {
    bool down_branch_is_valid;
    bool up_branch_is_valid;
    double dual_bound_on_down_branch;
    double dual_bound_on_up_branch;
    int num_simplex_iterations;
  };
  virtual absl::StatusOr<StronBranchResult> StrongBranchOnSingleColumn(
    int col, double primal_value, int iteration_limit) = 0;

  // ==========================================================================
  // Solution information getters.
  // ==========================================================================

  // Whether LP was solved.
  virtual bool IsSolved() const = 0;

  // Whether LP solution is stable (aka, "precise" or "reliable").
  //
  // This function should return true if the solution is reliable, i.e.,
  // feasible and optimal (or proven infeasible/unbounded) wrt the original
  // problem (within assumed numerical tolerance). E.g., the solution may be
  // optimal wrt a scaled version of the problem, but not feasible wrt 
  // original, unscaled problem -- in this case, the method should return false.
  virtual bool IsStable() const = 0;

  // Getters of the status of the computed solution.
  virtual bool IsOptimal() const = 0;
  virtual bool IsPrimalFeasible() const = 0;
  virtual bool IsPrimalInfeasible() const = 0;
  virtual bool IsPrimalUnbounded() const = 0;
  virtual bool IsDualFeasible() const = 0;
  virtual bool IsDualInfeasible() const = 0;
  virtual bool IsDualUnbounded() const = 0;

  // Whether LP was proven to have a primal unbounded ray (but not necessary
  // a primal feasible point). This does not mean that the solver knows and can
  // return the primal ray
  virtual bool ExistsPrimalRay() const = 0;

  // Whether LP was proven to have a primal unbounded ray (but not necessary
  // a primal feasible point). The solver knows and can return the primal ray.
  virtual bool HasPrimalRay() const = 0;

  // Whether LP was proven to have a dual unbounded ray (but not necessary
  // a dual feasible point). This does not mean that the solver knows and can
  // return the dual ray.
  virtual bool ExistsDualRay() const = 0;

  // Whether LP was proven to have a dual unbounded ray (but not necessary
  // a dual feasible point). The solver knows and can return the dual ray.
  virtual bool HasDualRay() const = 0;

  // Whether the solver hit various limits.
  // TODO(lpawel): Add support for a "deterministic effort" counter and limit.    
  virtual bool ObjectiveLimitIsExceeded() const = 0;
  virtual bool IterationLimitIsExceeded() const = 0;
  virtual bool TimeLimitIsExceeded() const = 0;

  // The number of simplex iterations of the call of primal or dual solve.
  virtual int num_iterations_of_last_solve() const = 0;

  // Returns the objective value of the solution. Can only be called if
  // `IsPrimalFeasible()` is true.
  virtual double GetObjectiveValue() const = 0;

  // Gets primal and dual solution vectors. Can only be called if
  // `IsOptimal()` is true.
  // TODO(lpawel): This should not fail (why would it?). But, if the
  // underlying solver requires it, we could change the output type to
  // `absl::StatusOr<DenseLpSolution>`.
  struct DenseLpSolution {
    std::vector<double> primal_values,
    std::vector<double> dual_values,
    std::vector<double> row_activities,
    std::vector<double> reduced_costs
  };
  virtual DenseLpSolution GetSolution() const = 0;
  
  // Returns the primal ray for unbounded LPs. Can only be called if
  // `HasPrimalRay()` is true.
  virtual absl::StatusOr<std::vector<double>> GetPrimalRay() const = 0;
    
  // Returns dual Farkas row multipliers from the proof of infeasibility.
  // Can only be called if `HasDualRay()` is true.
  virtual absl::StatusOr<std::vector<double>> GetDualFarkasMultipliers()
      const = 0;

  // ==========================================================================
  // Getters and setters of the basis.
  // ==========================================================================

  // Basis statuses of the model variables (if used for cols) or of the slack
  // variables (if used for rows).
  enum class BasisStatus {
    kBasic = 0,
    kAtLowerBound = 1,
    kAtUpperBound = 2,
    kFixed = 3,
    kFree = 4,
  };

  // Gets and sets the LP basis status for all cols and rows.
  struct LpBasis {
    std::vector<BasisStatus> col_basis_status;
    std::vector<BasisStatus> row_basis_status;
  };
  virtual absl::StatusOr<LpBasis> GetLpBasis() const = 0;
  virtual absl::Status SetBasis(const LpBasis& basis) = 0;

  // Returns the indices of the basic columns and rows. Basic column n is
  // indicated by value n. Basic row m is indicated by value -1-m.
  // The returned vector is dense and has size equal to `num_rows()`.
  virtual absl::StatusOr<std::vector<int>> GetBasisIndices() const = 0;

  // ==========================================================================
  // Getters of vectors in the inverted basis matrix.
  //
  // NOTE, LPInterface assumes that all slack variables were added by the
  // underlying LP Solver with coefficient equal to +1. If internally the
  // LP solver uses a -1 coefficient for a slack variable, then the row
  // associated with that slack variable must be negated here.  
  // ==========================================================================

  // Gets a dense col (or row) of the inverse basis matrix B^-1, which
  // corresponds to `GetBasisIndices()[basis_index]` model or slack variable.
  //
  // TODO(lpawel): It's good practive to mark methods as `const` if they should
  // not affect internal state. However -- I'm not sure whether for these
  // getters it's a reasonable expectation (there may be some things being
  // lazily computed and cached).
  advanced
  virtual absl::StaturOr<std::vector<double>> GetDenseRowOfBInverted(
      int basis_index) const = 0;
  virtual absl::StaturOr<std::vector<double>> GetDenseColumnOfBInverted(
      int basis_index) const = 0;

  // Get a dense col (or row) of the inverse basis matrix multiplied by
  // constraint matrix B^-1 * A, which corresponds to
  // `GetBasisIndices()[basis_index]` model or slack variable.
  virtual absl::StatusOr<std::vector<double>> GetDenseColumnOfBInvertedTimesA(
    int basis_index) const = 0;
  virtual absl::StatusOr<std::vector<double>> GetDenseRowOfBInvertedTimesA(
    int basis_index) const = 0;
  
  // Same as above but returns a sparse vector.
  // TODO(lpawel): Think whether we really need both dense and sparse getters.
  struct SparseVector {
    std::vector<int> indices;
    std::vector<double> values;
  };
  virtual absl::StaturOr<SparseVector> GetSparseRowOfBInverted(
      int basis_index) const = 0;
  virtual absl::StaturOr<SparseVector> GetSparseColumnOfBInverted(
      int basis_index) const = 0;
  virtual absl::StatusOr<SparseVector> GetSparseColumnOfBInvertedTimesA(
      int basis_index) const = 0;
  virtual absl::StatusOr<SparseVector> GetSparseRowOfBInvertedTimesA(
      int basis_index) const = 0;

  // ==========================================================================
  // Getters and setters of the parameters.
  // ==========================================================================

  // LP pricing strategy. The underlying LP solvers may not implement all those
  // pricing strategies, in which case trying to use it will result in an error
  // status.
  enum class PricingMethod {
    // The underlying LP solver should use its preferred strategy.
    kAuto = 0,

    // Full pricing.
    kFull = 1,        

    // Partial pricing.
    kPartial = 2,

    // Steepest edge pricing.
    kSteepestEdge = 3,

    // Steepest edge pricing without initial dual norms.
    kSteepEdgeQStart = 4, 

    // Devex pricing.
    kDevex = 5,
  };

  // LP scaling method to apply before solving.
  enum class ScalingMethod {
    // Do not use scaling.
    kOff = 0,

    // Use the preferred scaling 
    kAuto = 1,

    // TODO(lpawel): Add other scaling methods as needed.
  };

  // All of LP solver knobs. The underlying LP solvers may not implement all of
  // those (e.g., Glop is single threaded only) in which case an error status
  // is returned if an unsupported value is being set.
  struct Options {
    bool solve_from_scratch;
    ScalingMethod scaling_method;
    bool use_presolving;
    bool use_polishing;
    PricingMethod pricing_method;
    double primal_feasibility_tolerance;
    double dual_feasibility_tolerance;
    double markovitz_tolerance;
    int refactorization_period;
    double objective_limit;
    int simplex_iteration_limit;
    double time_limit_in_seconds;
    bool use_cpu_user_time;  // Otherwise, uses wallclock time.
    int num_threads;
    int random_seed;
    bool print_logs;
  };

  // `Options` are very small, thus it's perfectly fine to pass them by copy.
  virtual Options GetDefaultOptions() const = 0;
  virtual Options GetCurrentOptions() const = 0;
  virtual absl::Status SetOptions(const Options& options) = 0;

  // ==========================================================================
  // Various extra methods.
  // ==========================================================================

  // Returns the value treated as infinity by the underlying LP solver.
  virtual double Infinity() const = 0;

  // Whether `value` is plus or minus infinity.
  virtual bool IsInfinity(double value) const = 0;
};

} // namespace minimip
#endif // MINIMIP_SRC_LP_INTERFACE_LPI_H

