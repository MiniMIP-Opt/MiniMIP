#ifndef SRC_LP_INTERFACE_LPI_H_
#define SRC_LP_INTERFACE_LPI_H_

#include <string>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "src/lp_interface/lp_types.h"
#include "src/messagehandler/message_handler.h"

namespace minimip {

// LPInterface abstract class. This is used by MiniMIP to communicate with an
// underlying LP solver.
class LPInterface : private messagehandler {
 public:
  virtual ~LPInterface() = default;

  // ==========================================================================
  // LP model setters.
  // ==========================================================================

  virtual absl::Status LoadSparseColumnLP(
      LPObjectiveSense obj_sense,  // objective sense
      const std::vector<double>&
          objective_values,  // objective function values of columns
      const std::vector<double>& lower_bounds,      // lower bounds of columns
      const std::vector<double>& upper_bounds,      // upper bounds of columns
      std::vector<std::string>& col_names,          // column names
      const std::vector<double>& left_hand_sides,   // left hand sides of rows
      const std::vector<double>& right_hand_sides,  // right hand sides of rows
      std::vector<std::string>& row_names,          // row names
      const std::vector<SparseVector>& cols         // sparse columns
      ) = 0;

  // add column to the LP
  virtual absl::Status AddColumn(
      const SparseVector& col,       // column to be added
      const double lower_bound,      // lower bound of new column
      const double upper_bound,      // upper bound of new column
      const double objective_value,  // objective function value of new column
      std::string col_name           // column name
      ) = 0;

  // add columns to the LP
  virtual absl::Status AddColumns(
      const std::vector<SparseVector>& cols,    // columns to be added
      const std::vector<double>& lower_bounds,  // lower bounds of new columns
      const std::vector<double>& upper_bounds,  // upper bounds of new columns
      const std::vector<double>&
          objective_values,  // objective function values of new columns
      std::vector<std::string>& col_names  // column names
      ) = 0;

  // deletes all columns in the given range from LP
  virtual absl::Status DeleteColumns(
      int first_col,  // first column to be deleted
      int last_col    // last column to be deleted
      ) = 0;

  // add row to the LP
  virtual absl::Status AddRow(
      SparseVector row,        // row to be added
      double left_hand_side,   // left hand side of new row
      double right_hand_side,  // right hand side of new row
      std::string row_name     // row name
      ) = 0;

  // adds rows to the LP
  virtual absl::Status AddRows(
      const std::vector<SparseVector>& rows,  // number of rows to be added
      const std::vector<double>&
          left_hand_sides,  // left hand sides of new rows
      const std::vector<double>&
          right_hand_sides,                // right hand sides of new rows
      std::vector<std::string>& row_names  // row names
      ) = 0;

  // deletes all rows in the given range from LP
  virtual absl::Status DeleteRows(int first_row,  // first row to be deleted
                                  int last_row    // last row to be deleted
                                  ) = 0;

  // deletes rows from LP; the new position of a row must not be greater that
  // its old position
  virtual absl::Status DeleteRowSet(
      std::vector<bool>& deletion_status  // deletion status of rows
      ) = 0;

  // clears the whole LP
  virtual absl::Status Clear() = 0;

  // clears current LPInterface state (like basis information) of the solver
  virtual absl::Status ClearState() = 0;

  // change lower bound and upper bound of column
  virtual absl::Status SetColumnBounds(int col, double lower_bound,
                                       double upper_bound) = 0;

  // change left- and right-hand side of row
  virtual absl::Status SetRowSides(int row, double left_hand_side,
                                   double right_hand_side) = 0;

  // changes the objective sense
  virtual absl::Status SetObjectiveSense(
      LPObjectiveSense obj_sense  // new objective sense
      ) = 0;

  // changes a single objective value of a column in the LP
  virtual absl::Status SetObjectiveCoefficient(
      int col,  // column index to change objective value for
      double objective_coefficient  // new objective value
      ) = 0;

  // changes multiple objective values of columns in the LP
  virtual absl::Status SetObjectiveCoefficients(
      const std::vector<int>&
          indices,  // column indices to change objective value for
      const std::vector<double>&
          objective_coefficients  // new objective values for columns
  );
  // ==========================================================================
  // LP model getters.
  // ==========================================================================

  // gets the number of rows in the LP
  virtual int GetNumberOfRows() const = 0;

  // gets the number of columns in the LP
  virtual int GetNumberOfColumns() const = 0;

  // gets the number of non-zero elements in the LP constraint matrix
  virtual int GetNumberOfNonZeros() const = 0;

  // gets the objective sense of the LP
  virtual LPObjectiveSense GetObjectiveSense() const = 0;

  // gets columns from LP problem object
  virtual SparseVector GetSparseColumnCoefficients(int col) const = 0;

  // gets rows from LP problem object
  virtual SparseVector GetSparseRowCoefficients(int row) const = 0;

  // gets objective coefficient of column from LP problem object
  virtual double GetObjectiveCoefficient(int col) const = 0;

  // gets current lower bound of column from LP problem object
  virtual double GetLowerBound(int col) const = 0;

  // gets current upper bound of column from LP problem object
  virtual double GetUpperBound(int col) const = 0;

  // gets current left hand sides of row from LP problem object
  virtual double GetLeftHandSide(int row) const = 0;

  // gets current right hand sides of row from LP problem object
  virtual double GetRightHandSide(int row) const = 0;

  // gets the matrix coefficient of column and row from LP problem object
  virtual double GetMatrixCoefficient(int col,  // column number of coefficient
                                      int row   // row number of coefficient
  ) const = 0;

  // ==========================================================================
  // Solving methods.
  // ==========================================================================

  // calls primal simplex to solve the LP
  virtual absl::Status SolveLPWithPrimalSimplex() = 0;

  // calls dual simplex to solve the LP
  virtual absl::Status SolveLPWithDualSimplex() = 0;

  // start strong branching - call before any strong branching
  virtual absl::Status StartStrongBranching() = 0;

  // end strong branching - call after any strong branching
  virtual absl::Status EndStrongBranching() = 0;

  // performs strong branching iterations on one branching candidate
  struct StrongBranchResult {
    double dual_bound_down_branch;  // stores dual bound after branching
                                    // column down
    double dual_bound_up_branch;    // stores dual bound after branching
                                    // column up
    bool down_valid;  // whether the returned down value is a valid dual
                      // bound; otherwise, it can only be used as an
                      // estimate value
    bool up_valid;    // whether the returned up value is a valid dual bound;
                      // otherwise, it can only be used as an estimate value
    int iterations;   // stores total number of strong branching iterations
  };

  virtual absl::StatusOr<StrongBranchResult> StrongBranchValue(
      int col,             // column to apply strong branching on
      double primal_sol,   // current primal solution value of column
      int iteration_limit  // iteration limit for strong branchings
      ) = 0;

  // ==========================================================================
  // Solution information getters.
  // ==========================================================================

  // returns whether a solve method was called after the last modification of
  // the LP
  virtual bool IsSolved() const = 0;

  // returns true if current LP solution is stable
  //
  // This function should return true if the solution is reliable, i.e.,
  // feasible and optimal (or proven infeasible/unbounded) with respect to the
  // original problem. The optimality status might be with respect to a scaled
  // version of the problem, but the solution might not be feasible to the
  // unscaled original problem; in this case, minimip::LPInterface.IsStable()
  // should return false.
  virtual bool IsStable() const = 0;

  // returns true if LP was solved to optimality
  virtual bool IsOptimal() const = 0;

  // returns true if LP is proven to be primal feasible
  virtual bool IsPrimalFeasible() const = 0;

  // returns true if LP is proven to be primal infeasible
  virtual bool IsPrimalInfeasible() const = 0;

  // returns true if LP is proven to be primal unbounded
  virtual bool IsPrimalUnbounded() const = 0;

  // returns true if LP is proven to be dual feasible
  virtual bool IsDualFeasible() const = 0;

  // returns true if LP is proven to be dual infeasible
  virtual bool IsDualInfeasible() const = 0;

  // returns true if LP is proven to be dual unbounded
  virtual bool IsDualUnbounded() const = 0;

  // returns true if LP is proven to have a primal unbounded ray (but not
  // necessary a primal feasible point);
  //  this does not necessarily mean that the solver knows and can return the
  //  primal ray
  virtual bool ExistsPrimalRay() const = 0;

  // returns true if LP is proven to have a primal unbounded ray (but not
  // necessary a primal feasible point),
  //  and the solver knows and can return the primal ray
  virtual bool HasPrimalRay() const = 0;

  // returns true if LP is proven to have a dual unbounded ray (but not
  // necessary a dual feasible point); this does not necessarily mean that the
  // solver knows and can return the dual ray
  virtual bool ExistsDualRay() const = 0;

  // returns true if LP is proven to have a dual unbounded ray (but not
  // necessary a dual feasible point), and the solver knows and can return the
  // dual ray
  virtual bool HasDualRay() const = 0;

  // returns true if the objective limit was reached
  virtual bool ObjectiveLimitIsExceeded() const = 0;

  // returns true if the iteration limit was reached
  virtual bool IterationLimitIsExceeded() const = 0;

  // returns true if the time limit was reached
  virtual bool TimeLimitIsExceeded() const = 0;

  // gets objective value of solution
  virtual double GetObjectiveValue() = 0;

  // gets primal and dual solution vectors for feasible LPs
  //
  // Before calling these functions, the caller must ensure that the LP has been
  // solved to optimality, i.e., that minimip::LPInterface.IsOptimal() returns
  // true.

  // gets primal solution vector
  virtual absl::StatusOr<std::vector<double>> GetPrimalSolution() const = 0;

  // gets row activity vector
  virtual absl::StatusOr<std::vector<double>> GetRowActivity() const = 0;

  // gets dual solution vector
  virtual absl::StatusOr<std::vector<double>> GetDualSolution() const = 0;

  // gets reduced cost vector
  virtual absl::StatusOr<std::vector<double>> GetReducedCost() const = 0;

  // gets primal ray for unbounded LPs
  virtual absl::StatusOr<std::vector<double>> GetPrimalRay() const = 0;

  // gets dual Farkas proof for infeasibility
  virtual absl::StatusOr<std::vector<double>> GetDualFarkasMultiplier()
      const = 0;

  // gets the number of LP iterations of the last solve call
  virtual int GetIterations() const = 0;

  // ==========================================================================
  // Getters and setters of the basis.
  // ==========================================================================

  // gets current basis status for columns and rows
  virtual absl::StatusOr<std::vector<LPBasisStatus>> GetColumnBasisStatus()
      const = 0;
  virtual absl::StatusOr<std::vector<LPBasisStatus>> GetRowBasisStatus()
      const = 0;

  // sets current basis status for columns and rows
  virtual absl::Status SetBasisStatus(
      const std::vector<LPBasisStatus>& column_basis_status,
      const std::vector<LPBasisStatus>& row_basis_status) = 0;

  // returns the indices of the basic columns and rows; basic column n gives
  // value n, basic row m gives value -1-m
  virtual std::vector<int> GetBasisIndices() const = 0;

  // ==========================================================================
  // Getters of vectors in the inverted basis matrix.
  // ==========================================================================

  // get row of inverse basis matrix B^-1
  //
  // NOTE: The LP interface defines slack variables to have coefficient +1. This
  // means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables
  //       whose coefficient is -1, should be negated; see also the explanation
  //       in lpi.h.

  // make only dense available
  // virtual absl::StaturOr<const std::vector<double>&> GetSparseRowOfBInverted(

  virtual absl::StatusOr<SparseVector> GetSparseRowOfBInverted(
      int row_number) const = 0;

  // get column of inverse basis matrix B^-1
  //
  // NOTE: The LP interface defines slack variables to have coefficient +1. This
  // means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables
  //       whose coefficient is -1, should be negated
  //
  // column number of B^-1; this is NOT the number of the
  // column in the LP; you have to call
  // minimip::LPInterface.GetBasisIndices() to get the
  // array which links the B^-1 column numbers to the row
  // and column numbers of the LP! c must be between 0 and
  // num_rows-1, since the basis has the size num_rows *
  // num_rows
  virtual absl::StatusOr<SparseVector> GetSparseColumnOfBInverted(
      int col_number) const = 0;

  // get row of inverse basis matrix times constraint matrix B^-1 * A
  //
  // NOTE: The LP interface defines slack variables to have coefficient +1. This
  // means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables
  //       whose coefficient is -1, should be negated; see also the explanation
  //       in lpi.h.
  virtual absl::StatusOr<SparseVector> GetSparseRowOfBInvertedTimesA(
      int row_number) const = 0;

  // get column of inverse basis matrix times constraint matrix B^-1 * A
  //
  // NOTE: The LP interface defines slack variables to have coefficient +1. This
  // means that if, internally, the LP solver
  //       uses a -1 coefficient, then rows associated with slacks variables
  //       whose coefficient is -1, should be negated; see also the explanation
  //       in lpi.h.
  virtual absl::StatusOr<SparseVector> GetSparseColumnOfBInvertedTimesA(
      int col_number) const = 0;

  // ==========================================================================
  // Getters and setters of the parameters.
  // ==========================================================================

  // gets integer parameter of LP
  virtual absl::StatusOr<int> GetIntegerParameter(
      LPParameter type  // parameter number
  ) const = 0;

  // sets integer parameter of LP
  virtual absl::Status SetIntegerParameter(
      LPParameter type,  // parameter number
      int param_val      // parameter value
      ) = 0;

  // gets floating point parameter of LP
  virtual absl::StatusOr<double> GetRealParameter(
      LPParameter type  // parameter number
  ) const = 0;

  // sets floating point parameter of LP
  virtual absl::Status SetRealParameter(LPParameter type,  // parameter number
                                        double param_val   // parameter value
                                        ) = 0;

  // ==========================================================================
  // Numerical methods.
  // ==========================================================================

  // returns value treated as infinity in the LP solver
  virtual double Infinity() const = 0;

  // checks if given value is treated as infinity in the LP solver
  virtual bool IsInfinity(double value  // value to be checked for infinity
  ) const = 0;

  // ==========================================================================
  // File interface methods.
  // ==========================================================================

  // reads LP from a file
  virtual absl::Status ReadLP(const char* file_name  // file name
                              ) = 0;

  // writes LP to a file
  virtual absl::Status WriteLP(const char* file_name  // file name
  ) const = 0;
};

}  // namespace minimip
#endif  // SRC_LP_INTERFACE_LPI_H_
