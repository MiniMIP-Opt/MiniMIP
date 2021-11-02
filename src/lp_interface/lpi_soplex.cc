
#include "src/lp_interface/lpi_soplex.h"

// check the return value of setParam methods
#define CHECK_SOPLEX_PARAM(x)                                            \
  if (!x) {                                                              \
    MiniMIPmessagePrintWarning("SoPlex: unsupported parameter value\n"); \
  }

// remember the original value of the MiniMIP_DEBUG define and undefine it
#ifdef MINIMIP_DEBUG
#define ___DEBUG
#undef MINIMIP_DEBUG
#endif

// reset the MiniMIP_DEBUG define to its original MiniMIP value
#undef MINIMIP_DEBUG
#ifdef ___DEBUG
#define MINIMIP_DEBUG
#undef ___DEBUG
#endif

#define SOPLEX_VERBLEVEL 0  // verbosity level for LPINFO

// Macro for a single SoPlex call for which exceptions have to be catched -
// return an LP error. We make no distinction between different exception types,
// e.g., between memory allocation and other exceptions.
#ifndef NDEBUG
#define SOPLEX_TRY(x)                                                          \
  do {                                                                         \
    try {                                                                      \
      (x);                                                                     \
    } catch (const SPxMemoryException& E) {                                    \
      std::string s = E.what();                                                \
      MiniMIPerrorMessage("SoPlex threw a memory exception: %s\n", s.c_str()); \
      return absl::Status(absl::StatusCode::kInternal, "Error");               \
    } catch (const SPxException& E) {                                          \
      std::string s = E.what();                                                \
      MiniMIPerrorMessage("SoPlex threw an exception: %s\n", s.c_str());       \
      return absl::Status(absl::StatusCode::kInternal, "LP Error");            \
    }                                                                          \
  } while (false)

#else
#define SOPLEX_TRY(x)                                                          \
  do {                                                                         \
    try {                                                                      \
      (x);                                                                     \
    } catch (const SPxMemoryException& E) {                                    \
      std::string s = E.what();                                                \
      MiniMIPerrorMessage("SoPlex threw a memory exception: %s\n", s.c_str()); \
      return absl::Status(absl::StatusCode::kInternal, "Error");               \
    } catch (const SPxException&) {                                            \
      return absl::Status(absl::StatusCode::kInternal, "LP Error");            \
    }                                                                          \
  } while (FALSE)
#endif

// Macro for a single SoPlex call for which exceptions have to be catched -
// abort if they arise. MINIMIP_ABORT() is not accessible here.
#define SOPLEX_TRY_ABORT(x)                                              \
  do {                                                                   \
    try {                                                                \
      (x);                                                               \
    } catch (const SPxException& E) {                                    \
      std::string s = E.what();                                          \
      MiniMIPerrorMessage("SoPlex threw an exception: %s\n", s.c_str()); \
      abort();                                                           \
    }                                                                    \
  } while (FALSE)

namespace minimip {

LPSoplexInterface::LPSoplexInterface()
    : spx_(new soplex::SoPlex),
      pricing_(LPPricing::kDefault),
      solved_(false),
      lp_info_(false),
      from_scratch_(false),
      col_basis_status_(0),
      row_basis_status_(0) {}

LPSoplexInterface::~LPSoplexInterface() { delete spx_; }

// provides access for temporary storage of basis status of rows
DataArray<SPxSolver::VarStatus>& LPSoplexInterface::RowsBasisStatus() {
  return row_basis_status_;
}
// provides access for temporary storage of basis status of columns
DataArray<SPxSolver::VarStatus>& LPSoplexInterface::ColumnsBasisStatus() {
  return col_basis_status_;
}

void LPSoplexInterface::SetFromScratch(bool from_scratch) {
  from_scratch_ = from_scratch;
}

void LPSoplexInterface::SetLPInfo(bool lp_info) { lp_info_ = lp_info; }

// get objective limit according to objective sense
double LPSoplexInterface::GetObjectiveLimit() const {
  return (spx_->intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE)
             ? spx_->realParam(SoPlex::OBJLIMIT_UPPER)
             : spx_->realParam(SoPlex::OBJLIMIT_LOWER);
}

double LPSoplexInterface::FeasibilityTolerance() const {
  return spx_->realParam(SoPlex::FEASTOL);
}

// return optimality tolerance
double LPSoplexInterface::OptimalityTolerance() const {
  return spx_->realParam(SoPlex::OPTTOL);
}

// set feasibility tolerance and store value in case SoPlex only accepts a
// larger tolerance
void LPSoplexInterface::SetFeasibilityTolerance(const Real d) {
  assert(spx_->setRealParam(SoPlex::FEASTOL, d));
}

// set optimality tolerance and store value in case SoPlex only accepts a larger
// tolerance
void LPSoplexInterface::SetOptimalityTolerance(const Real d) {
  assert(spx_->setRealParam(SoPlex::OPTTOL, d));
}

void LPSoplexInterface::TrySolve(bool print_warning = true) {
  Real time_spent;
  Real time_limit;

  try {
    static_cast<void>(spx_->optimize());
  } catch (const SPxException& x) {
    std::string s = x.what();
    if (print_warning) {
      //  MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
      //  %s\n", s.c_str());
    }
    // since it is not clear if the status in SoPlex are set correctly we want
    // to make sure that if an error is thrown the status is not OPTIMAL
    // anymore.
    assert(spx_->status() != SPxSolver::OPTIMAL);
  }
  assert(spx_->intParam(SoPlex::ITERLIMIT) < 0 ||
         spx_->numIterations() <= spx_->intParam(SoPlex::ITERLIMIT));

  // update time limit
  time_spent = spx_->solveTime();

  // get current time limit
  if (time_spent > 0) {
    time_limit = spx_->realParam(SoPlex::TIMELIMIT);
    if (time_limit > time_spent)
      time_limit -= time_spent;
    else
      time_limit = 0;

    // set new time limit
    assert(time_limit >= 0);

    assert(spx_->setRealParam(SoPlex::TIMELIMIT, time_limit));
  }
}

bool LPSoplexInterface::CheckConsistentBounds() const {
  for (int i = 0; i < spx_->numColsReal(); ++i) {
    if (spx_->lowerReal(i) >
        spx_->upperReal(i) + spx_->realParam(SoPlex::EPSILON_ZERO)) {
      MiniMIPerrorMessage(
          "inconsistent bounds on column %d: lower=%.17g, upper=%.17g\n", i,
          spx_->lowerReal(i), spx_->upperReal(i));
      return false;
    }
  }
  return true;
}

bool LPSoplexInterface::CheckConsistentSides() const {
  for (int i = 0; i < spx_->numRowsReal(); ++i) {
    if (spx_->lhsReal(i) >
        spx_->rhsReal(i) + spx_->realParam(SoPlex::EPSILON_ZERO)) {
      MiniMIPerrorMessage(
          "inconsistent sides on row %d: lhs=%.17g, rhs=%.17g\n", i,
          spx_->lhsReal(i), spx_->rhsReal(i));
      return false;
    }
  }
  return true;
}

SPxSolver::Status LPSoplexInterface::doSolve(bool print_warning = true) {
  SPxOut::Verbosity verbosity;
  SPxSolver::Status spx_status;

  // store and set verbosity
  verbosity = spx_->spxout.getVerbosity();
  spx_->spxout.setVerbosity(
      (SPxOut::Verbosity)(GetLPInfo() ? SOPLEX_VERBLEVEL : 0));

  assert(CheckConsistentBounds());
  assert(CheckConsistentSides());

  TrySolve(print_warning);
  spx_status = spx_->status();

  // restore verbosity
  spx_->spxout.setVerbosity(verbosity);

  return spx_status;
}

bool LPSoplexInterface::GetFromScratch() const { return from_scratch_; }

// @todo member variable?
bool LPSoplexInterface::GetLPInfo() const { return lp_info_; }

absl::Status LPSoplexInterface::SoPlexSolve() {
  // store and set verbosity
  SPxOut::Verbosity verbosity = spx_->spxout.getVerbosity();
  spx_->spxout.setVerbosity(
      (SPxOut::Verbosity)(GetLPInfo() ? SOPLEX_VERBLEVEL : 0));

  MiniMIPdebugMessage("calling SoPlex solve(): %d cols, %d rows\n",
                      spx_->numColsReal(), spx_->numRowsReal());

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

  // delete starting basis if solving from scratch
  if (GetFromScratch()) {
    try {
      spx_->clearBasis();
    }
#ifndef NDEBUG
    catch (const SPxException& x) {
      std::string s = x.what();
      // MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
      // %s\n", s.c_str());
#else
    catch (const SPxException&) {
#endif
      assert(spx_->status() != SPxSolver::OPTIMAL);
      return absl::Status(absl::StatusCode::kInternal, "LP Error");
    }
  }
  assert(!GetFromScratch() || spx_->status() == SPxSolver::NO_PROBLEM);

  SPxSolver::Status status = doSolve();
  MiniMIPdebugMessage(" -> SoPlex status: %d, basis status: %d\n",
                      spx_->status(), spx_->basisStatus());
  solved_ = true;

  // restore verbosity
  spx_->spxout.setVerbosity(verbosity);

  switch (status) {
    case SPxSolver::ABORT_TIME:
    case SPxSolver::ABORT_ITER:
    case SPxSolver::ABORT_VALUE:
    case SPxSolver::SINGULAR:
    case SPxSolver::REGULAR:
    case SPxSolver::UNKNOWN:
    case SPxSolver::OPTIMAL:
    case SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS:
    case SPxSolver::UNBOUNDED:
    case SPxSolver::INFEASIBLE:
      return absl::OkStatus();
    default:
      return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }
}

// save the current basis
void LPSoplexInterface::SavePreStrongbranchingBasis() {
  row_basis_status_.reSize(spx_->numRowsReal());
  col_basis_status_.reSize(spx_->numColsReal());

  try {
    spx_->getBasis(row_basis_status_.get_ptr(), col_basis_status_.get_ptr());
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    // MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
    // %s\n", s.c_str());

    // since it is not clear if the status in SoPlex are set correctly we want
    // to make sure that if an error is thrown the status is not OPTIMAL
    // anymore.
    assert(spx_->status() != SPxSolver::OPTIMAL);
  }
#else
  catch (const SPxException&) {
  }
#endif
}

// restore basis
void LPSoplexInterface::RestorePreStrongbranchingBasis() {
  assert(row_basis_status_.size() == spx_->numRowsReal());
  assert(col_basis_status_.size() == spx_->numColsReal());

  try {
    spx_->setBasis(row_basis_status_.get_ptr(), col_basis_status_.get_ptr());
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    // MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
    // %s\n", s.c_str());
#else
  catch (const SPxException&) {
#endif
    // since it is not clear if the status in SoPlex are set correctly we want
    // to make sure that if an error is thrown the status is not OPTIMAL
    // anymore.
    assert(spx_->status() != SPxSolver::OPTIMAL);
  }
}

void LPSoplexInterface::InvalidateSolution() { solved_ = false; }

bool LPSoplexInterface::PreStrongBranchingBasisFreed() const {
  return ((row_basis_status_.size() == 0) && (col_basis_status_.size() == 0));
}

void LPSoplexInterface::FreePreStrongBranchingBasis() {
  row_basis_status_.clear();
  col_basis_status_.clear();
}

absl::Status LPSoplexInterface::StrongBranch(
    int col,                         // column to apply strong branching on
    double primal_sol,               // current primal solution value of column
    int iteration_limit,             // iteration limit for strong branchings
    double& dual_bound_down_branch,  // stores dual bound after branching column
                                     // down
    double&
        dual_bound_up_branch,  // stores dual bound after branching column up
    bool& down_valid,  // stores whether the returned down value is a valid dual
                       // bound; otherwise, it can only be used as an estimate
                       // value
    bool&
        up_valid,  // stores whether the returned up value is a valid dual
                   // bound; otherwise, it can only be used as an estimate value
    int&
        iterations  // stores total number of strong branching iterations, or -1
) {
  SPxSolver::Status status;
  double old_lb;
  double old_ub;
  double new_lb;
  double new_ub;
  bool from_parent_basis;
  bool error;
  int old_iter_limit;
  SPxOut::Verbosity verbosity;

  // store and set verbosity
  verbosity = spx_->spxout.getVerbosity();
  spx_->spxout.setVerbosity(
      (SPxOut::Verbosity)(GetLPInfo() ? SOPLEX_VERBLEVEL : 0));

  MiniMIPdebugMessage("calling StrongBranch() on variable %d (%d iterations)\n",
                      col, iteration_limit);

#ifndef STRONGBRANCH_RESTOREBASIS
  from_parent_basis = false;
#endif
  error          = false;
  old_iter_limit = spx_->intParam(SoPlex::ITERLIMIT);

  // get current bounds of column
  old_lb     = spx_->lowerReal(col);
  old_ub     = spx_->upperReal(col);
  down_valid = false;
  up_valid   = false;
  iterations = 0;

  // set the algorithm type to use dual simplex
  static_cast<void>(
      spx_->setIntParam(SoPlex::ALGORITHM, SoPlex::ALGORITHM_DUAL));

  // down branch
  new_ub = EPSCEIL(primal_sol - 1.0, FeasibilityTolerance());

  if (new_ub >= old_lb - 0.5) {
    MiniMIPdebugMessage(
        "strong branching down on x%d (%g) with %d iterations\n", col,
        primal_sol, iteration_limit);

    spx_->changeUpperReal(col, new_ub);
    assert(spx_->lowerReal(col) <= spx_->upperReal(col));

    static_cast<void>(spx_->setIntParam(SoPlex::ITERLIMIT, iteration_limit));
    do {
#ifndef STRONGBRANCH_RESTOREBASIS
      bool repeat_strong_branching;
#endif
      status = spx_->optimize();

      MiniMIPdebugMessage(" --> Terminate with status %d\n", status);
      switch (status) {
        case SPxSolver::OPTIMAL:
          dual_bound_down_branch = spx_->objValueReal();
          down_valid             = true;
          MiniMIPdebugMessage(" --> Terminate with value %f\n",
                              dual_bound_down_branch);
          break;
        case SPxSolver::ABORT_TIME:  // SoPlex does not return a proven dual
                                     // bound, if it is aborted
        case SPxSolver::ABORT_ITER:
        case SPxSolver::ABORT_CYCLING:
        case SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS:
          dual_bound_down_branch = spx_->objValueReal();
          break;
        case SPxSolver::ABORT_VALUE:
        case SPxSolver::INFEASIBLE:
          dual_bound_down_branch = GetObjectiveLimit();
          down_valid             = true;
          break;
        default:
          error = true;
          break;
      }
      iterations += spx_->numIterations();

#ifdef STRONGBRANCH_RESTOREBASIS
      // we restore the pre-strong-branching basis by default (and don't solve
      // again)
      assert(!PreStrongBranchingBasisFreed());
      RestorePreStrongbranchingBasis();
      from_parent_basis = false;
#else
      // if cycling or singular basis occured and we started not from the
      // pre-strong-branching basis, then we restore the pre-strong-branching
      // basis and try again with reduced iteration limit
      repeat_strong_branching =
          ((status == SPxSolver::ABORT_CYCLING ||
            status == SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS ||
            status == SPxSolver::SINGULAR) &&
           !from_parent_basis && spx_->numIterations() < iteration_limit);

      if (repeat_strong_branching) {
        MiniMIPdebugMessage(
            " --> Repeat strong branching down with %d iterations after "
            "restoring basis\n",
            iteration_limit - spx_->numIterations());
        spx_->setIntParam(SoPlex::ITERLIMIT,
                          iteration_limit - spx_->numIterations());
        RestorePreStrongbranchingBasis();
        from_parent_basis = true;
        error = false;
      } else {
        from_parent_basis = false;
      }
#endif
    } while (from_parent_basis);

    spx_->changeUpperReal(col, old_ub);
    assert(spx_->lowerReal(col) <= spx_->upperReal(col));
  } else {
    dual_bound_down_branch = GetObjectiveLimit();
    down_valid             = true;
  }

  // up branch
  if (!error) {
    new_lb = EPSFLOOR(primal_sol + 1.0, FeasibilityTolerance());
    if (new_lb <= old_ub + 0.5) {
      MiniMIPdebugMessage(
          "strong branching  up  on x%d (%g) with %d iterations\n", col,
          primal_sol, iteration_limit);

      spx_->changeLowerReal(col, new_lb);
      assert(spx_->lowerReal(col) <= spx_->upperReal(col));

      static_cast<void>(spx_->setIntParam(SoPlex::ITERLIMIT, iteration_limit));
      do {
#ifndef STRONGBRANCH_RESTOREBASIS
        bool repeat_strong_branching;
#endif
        status = spx_->optimize();

        MiniMIPdebugMessage(" --> Terminate with status %d\n", status);
        switch (status) {
          case SPxSolver::OPTIMAL:
            dual_bound_up_branch = spx_->objValueReal();
            up_valid             = true;
            MiniMIPdebugMessage(" --> Terminate with value %f\n",
                                spx_->objValueReal());
            break;
          case SPxSolver::ABORT_TIME:  // SoPlex does not return a proven dual
                                       // bound, if it is aborted
          case SPxSolver::ABORT_ITER:
          case SPxSolver::ABORT_CYCLING:
          case SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS:
            dual_bound_up_branch = spx_->objValueReal();
            break;
          case SPxSolver::ABORT_VALUE:
          case SPxSolver::INFEASIBLE:
            dual_bound_up_branch = GetObjectiveLimit();
            up_valid             = true;
            break;
          default:
            error = true;
            break;
        }
        iterations += spx_->numIterations();

#ifdef STRONGBRANCH_RESTOREBASIS
        // we restore the pre-strong-branching basis by default (and don't solve
        // again)
        assert(!PreStrongBranchingBasisFreed());
        RestorePreStrongbranchingBasis();
        from_parent_basis = false;
#else
        // if cycling or singular basis occured and we started not from the
        // pre-strong-branching basis, then we restore the pre-strong-branching
        // basis and try again with reduced iteration limit
        repeat_strong_branching =
            ((status == SPxSolver::ABORT_CYCLING ||
              status == SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS ||
              status == SPxSolver::SINGULAR) &&
             !from_parent_basis && spx_->numIterations() < iteration_limit);

        if (repeat_strong_branching) {
          MiniMIPdebugMessage(
              " --> Repeat strong branching  up  with %d iterations after "
              "restoring basis\n",
              iteration_limit - spx_->numIterations());
          RestorePreStrongbranchingBasis();
          spx_->setIntParam(SoPlex::ITERLIMIT,
                            iteration_limit - spx_->numIterations());
          error = false;
          from_parent_basis = true;
        } else {
          from_parent_basis = false;
        }
#endif
      } while (from_parent_basis);

      spx_->changeLowerReal(col, old_lb);
      assert(spx_->lowerReal(col) <= spx_->upperReal(col));
    } else {
      dual_bound_up_branch = GetObjectiveLimit();
      up_valid             = true;
    }
  }
  // reset old iteration limit
  static_cast<void>(spx_->setIntParam(SoPlex::ITERLIMIT, old_iter_limit));

  // restore verbosity
  spx_->spxout.setVerbosity(verbosity);

  if (error) {
    MiniMIPdebugMessage("StrongBranch() returned SoPlex status %d\n",
                        static_cast<int>(status));
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }
  return absl::OkStatus();
}

// ==========================================================================
// LP model setters.
// ==========================================================================

absl::Status LPSoplexInterface::LoadColumnLP(
    LPObjectiveSense obj_sense,  // objective sense
    int num_cols,                // number of columns
    const std::vector<double>&
        objective_values,  // objective function values of columns
    const std::vector<double>& lower_bounds,      // lower bounds of columns
    const std::vector<double>& upper_bounds,      // upper bounds of columns
    std::vector<std::string>& col_names,          // column names
    int num_rows,                                 // number of rows
    const std::vector<double>& left_hand_sides,   // left hand sides of rows
    const std::vector<double>& right_hand_sides,  // right hand sides of rows
    std::vector<std::string>& row_names,          // row names
    int num_non_zeros,  // number of non-zero elements in the constraint matrix
    const std::vector<int>& begin_cols,  // start index of each column in
                                         // row_indices- and vals-array
    const std::vector<int>&
        row_indices,                 // row indices of constraint matrix entries
    const std::vector<double>& vals  // values of constraint matrix entries
) {
#ifndef NDEBUG
  for (int j = 0; j < num_non_zeros; j++) assert(vals[j] != 0);
#endif

  MiniMIPdebugMessage("calling LoadColumnLP()\n");

  assert(PreStrongBranchingBasisFreed());

  try {
    LPRowSet rows(num_rows);
    DSVector empty_vector(0);
    int i;

    spx_->clearLPReal();

    // set objective sense
    static_cast<void>(spx_->setIntParam(
        SoPlex::OBJSENSE, (obj_sense == LPObjectiveSense::kMinimization
                               ? SoPlex::OBJSENSE_MINIMIZE
                               : SoPlex::OBJSENSE_MAXIMIZE)));

    // create empty rows with given sides
    for (i = 0; i < num_rows; ++i)
      rows.add(left_hand_sides[i], empty_vector, right_hand_sides[i]);
    spx_->addRowsReal(rows);

    // create column vectors with coefficients and bounds
    MINIMIP_CALL(AddColumns(num_cols, objective_values, lower_bounds,
                            upper_bounds, col_names, num_non_zeros, begin_cols,
                            row_indices, vals));
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    // MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
    // %s\n", s.c_str());
#else
  catch (const SPxException&) {
#endif
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }
  return absl::OkStatus();
}

// adds columns to the LP
//
// NOTE: The indices array is not checked for duplicates, problems may appear if
// indices are added more than once.
absl::Status LPSoplexInterface::AddColumns(
    int num_cols,  // number of columns to be added
    const std::vector<double>&
        objective_values,  // objective function values of new columns
    const std::vector<double>& lower_bounds,  // lower bounds of new columns
    const std::vector<double>& upper_bounds,  // upper bounds of new columns
    std::vector<std::string>& col_names,      // column names
    int num_non_zeros,  // number of non-zero elements to be added to the
                        // constraint matrix
    const std::vector<int>&
        begin_cols,  // start index of each column in indices- and vals-array
    const std::vector<int>&
        indices,                     // row indices of constraint matrix entries
    const std::vector<double>& vals  // values of constraint matrix entries
) {
  MiniMIPdebugMessage("calling AddColumns()\n");

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

#ifndef NDEBUG
  if (num_non_zeros > 0) {
    // perform check that no new rows are added - this is likely to be a mistake
    int num_rows = spx_->numRowsReal();
    for (int j = 0; j < num_non_zeros; ++j) {
      assert(0 <= indices[j] && indices[j] < num_rows);
      assert(vals[j] != 0.0);
    }
  }
#endif

  try {
    std::vector<int> integer_indices(indices.begin(), indices.end());
    LPColSet cols(num_cols);
    DSVector col_Vector(num_cols);
    int start;
    int last;
    int i;

    // create column vectors with coefficients and bounds
    for (i = 0; i < num_cols; ++i) {
      col_Vector.clear();
      if (num_non_zeros > 0) {
        start = begin_cols[i];
        last  = (i == num_cols - 1 ? num_non_zeros : begin_cols[i + 1]);
        col_Vector.add(last - start, &integer_indices[start], &vals[start]);
      }
      cols.add(objective_values[i], lower_bounds[i], col_Vector,
               upper_bounds[i]);
    }
    spx_->addColsReal(cols);
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    // MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
    // %s\n", s.c_str());
#else
  catch (const SPxException&) {
#endif
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }

  return absl::OkStatus();
}

// deletes all columns in the given range from LP
absl::Status LPSoplexInterface::DeleteColumns(
    int first_col,  // first column to be deleted
    int last_col    // last column to be deleted
) {
  MiniMIPdebugMessage("calling DeleteColumns()\n");

  assert(0 <= first_col && first_col <= last_col &&
         last_col < spx_->numColsReal());

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

  SOPLEX_TRY(spx_->removeColRangeReal(first_col, last_col));

  return absl::OkStatus();
}

// adds rows to the LP
//
// NOTE: The indices array is not checked for duplicates, problems may appear if
// indices are added more than once.
absl::Status LPSoplexInterface::AddRows(
    int num_rows,                                // number of rows to be added
    const std::vector<double>& left_hand_sides,  // left hand sides of new rows
    const std::vector<double>&
        right_hand_sides,                 // right hand sides of new rows
    std::vector<std::string>& row_names,  // row names
    int num_non_zeros,  // number of non-zero elements to be added to the
                        // constraint matrix
    const std::vector<int>&
        begin_rows,  // start index of each row in indices- and vals-array
    const std::vector<int>&
        indices,  // column indices of constraint matrix entries
    const std::vector<double>& vals  // values of constraint matrix entries
) {
  MiniMIPdebugMessage("calling AddRows()\n");

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

#ifndef NDEBUG
  if (num_non_zeros > 0) {
    // perform check that no new columns are added - this is likely to be a
    // mistake
    int num_cols = spx_->numColsReal();
    for (int j = 0; j < num_non_zeros; ++j) {
      assert(vals[j] != 0.0);
      assert(0 <= indices[j] && indices[j] < num_cols);
    }
  }
#endif

  try {
    LPRowSet rows(num_rows);
    DSVector row_Vector;
    int start;
    int last;
    int i;
    std::vector<int> integer_indices(indices.begin(), indices.end());

    // create row vectors with given sides
    for (i = 0; i < num_rows; ++i) {
      row_Vector.clear();
      if (num_non_zeros > 0) {
        start = begin_rows[i];
        last  = (i == num_rows - 1 ? num_non_zeros : begin_rows[i + 1]);
        row_Vector.add(last - start, &integer_indices[start], &vals[start]);
      }
      rows.add(left_hand_sides[i], row_Vector, right_hand_sides[i]);
    }
    spx_->addRowsReal(rows);
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    // MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
    // %s\n", s.c_str());
#else
  catch (const SPxException&) {
#endif
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }
  return absl::OkStatus();
}

// deletes all rows in the given range from LP
absl::Status LPSoplexInterface::DeleteRows(
    int first_row,  // first row to be deleted
    int last_row    // last row to be deleted
) {
  MiniMIPdebugMessage("calling DeleteRows()\n");

  assert(0 <= first_row && first_row <= last_row &&
         last_row < spx_->numRowsReal());

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

  SOPLEX_TRY(spx_->removeRowRangeReal(first_row, last_row));

  return absl::OkStatus();
}

// deletes rows from LP; the new position of a row must not be greater that its
// old position
absl::Status LPSoplexInterface::DeleteRowSet(
    std::vector<bool>& deletion_status  // deletion status of rows
) {
  int num_rows;
  int i;

  MiniMIPdebugMessage("calling DeleteRowSet()\n");

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

  std::vector<int> int_deletion_status(deletion_status.begin(),
                                       deletion_status.end());

  num_rows = spx_->numRowsReal();

  // SoPlex removeRows() method deletes the rows with deletion_status[i] < 0, so
  // we have to negate the values
  for (i = 0; i < num_rows; ++i) int_deletion_status[i] *= -1;

  SOPLEX_TRY(spx_->removeRowsReal(int_deletion_status.data()));

  return absl::OkStatus();
}

// clears the whole LP
absl::Status LPSoplexInterface::Clear() {
  MiniMIPdebugMessage("calling Clear()\n");

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());
  SOPLEX_TRY(spx_->clearLPReal());

  return absl::OkStatus();
}

// clears current LP Interface state (like basis information) of the solver
absl::Status LPSoplexInterface::ClearState() {
  MiniMIPdebugMessage("calling ClearState()\n");

  try {
    spx_->clearBasis();
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    //  MiniMIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception:
    //  %s\n", s.c_str());
#else
  catch (const SPxException&) {
#endif
    assert(spx_->status() != SPxSolver::OPTIMAL);
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }
  return absl::OkStatus();
}

// changes lower and upper bounds of columns
absl::Status LPSoplexInterface::SetColumnBounds(
    int num_cols,                     // number of columns to change bounds for
    const std::vector<int>& indices,  // column indices
    const std::vector<double>& lower_bounds,  // values for the new lower bounds
    const std::vector<double>& upper_bounds   // values for the new upper bounds
) {
  int i;

  MiniMIPdebugMessage("calling SetColumnBounds()\n");

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

  try {
    for (i = 0; i < num_cols; ++i) {
      assert(0 <= indices[i] && indices[i] < spx_->numColsReal());

      if (IsInfinity(lower_bounds[i])) {
        MiniMIPerrorMessage(
            "LP Error: fixing lower bound for variable %d to infinity.\n",
            indices[i]);
        return absl::Status(absl::StatusCode::kInternal, "LP Error");
      }
      if (IsInfinity(-upper_bounds[i])) {
        MiniMIPerrorMessage(
            "LP Error: fixing upper bound for variable %d to -infinity.\n",
            indices[i]);
        return absl::Status(absl::StatusCode::kInternal, "LP Error");
      }
      spx_->changeBoundsReal(indices[i], lower_bounds[i], upper_bounds[i]);
      assert(spx_->lowerReal(indices[i]) <=
             spx_->upperReal(indices[i]) +
                 spx_->realParam(SoPlex::EPSILON_ZERO));
    }
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    // MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
    // %s\n", s.c_str());
#else
  catch (const SPxException&) {
#endif
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }
  return absl::OkStatus();
}

// changes left and right hand sides of rows
absl::Status LPSoplexInterface::SetRowSides(
    int num_rows,                     // number of rows to change sides for
    const std::vector<int>& indices,  // row indices
    const std::vector<double>&
        left_hand_sides,  // new values for left hand sides
    const std::vector<double>&
        right_hand_sides  // new values for right hand sides
) {
  int i;

  MiniMIPdebugMessage("calling SetRowSides()\n");

  if (num_rows <= 0) return absl::OkStatus();

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

  try {
    for (i = 0; i < num_rows; ++i) {
      assert(0 <= indices[i] && indices[i] < spx_->numRowsReal());
      spx_->changeRangeReal(indices[i], left_hand_sides[i],
                            right_hand_sides[i]);
      assert(spx_->lhsReal(indices[i]) <=
             spx_->rhsReal(indices[i]) + spx_->realParam(SoPlex::EPSILON_ZERO));
    }
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    // MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
    // %s\n", s.c_str());
#else
  catch (const SPxException&) {
#endif
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }
  return absl::OkStatus();
}

// changes the objective sense
absl::Status LPSoplexInterface::SetObjectiveSense(
    LPObjectiveSense obj_sense  // new objective sense
) {
  MiniMIPdebugMessage("calling SetObjectiveSense()\n");

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

  SOPLEX_TRY(static_cast<void>(spx_->setIntParam(
      SoPlex::OBJSENSE, obj_sense == LPObjectiveSense::kMinimization
                            ? SoPlex::OBJSENSE_MINIMIZE
                            : SoPlex::OBJSENSE_MAXIMIZE)));

  return absl::OkStatus();
}

// changes objective values of columns in the LP
absl::Status LPSoplexInterface::SetObjectiveCoefficient(
    const int col, // column index to change objective value for
    const double objective_coefficient // new objective values for column
) {
  spx_->changeObjReal(col, objective_coefficient);
  return absl::OkStatus();
}

// changes objective values of columns in the LP
absl::Status LPSoplexInterface::SetObjectiveCoefficients(
    const std::vector<int>&
        indices,  // column indices to change objective value for
    const std::vector<double>&
        objective_coefficients  // new objective values for columns
) {
  size_t i;

  MiniMIPdebugMessage("calling SetObjectiveCoefficients()\n");

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

  try {
    for (i = 0; i < indices.size(); ++i) {
      assert(0 <= indices[i] && indices[i] < spx_->numColsReal());
      spx_->changeObjReal(indices[i], objective_coefficients[i]);
    }
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    // MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
    // %s\n", s.c_str());
#else
  catch (const SPxException&) {
#endif
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }
  return absl::OkStatus();
}

// ==========================================================================
// LP model getters.
// ==========================================================================

// gets the number of rows in the LP
int LPSoplexInterface::GetNumberOfRows() const {
  MiniMIPdebugMessage("calling GetNumberOfRows()\n");

  return spx_->numRowsReal();
}

// gets the number of columns in the LP
int LPSoplexInterface::GetNumberOfColumns() const {
  MiniMIPdebugMessage("calling GetNumberOfColumns()\n");

  return spx_->numColsReal();
}

// gets the number of nonzero elements in the LP constraint matrix
int LPSoplexInterface::GetNumberOfNonZeros() const {
  int i;
  // SoPlex has no direct method to return the number of nonzeros, so we have to
  // count them manually
  int num_non_zeros = 0;

  MiniMIPdebugMessage("calling GetNumberOfNonZeros()\n");

  if (spx_->numRowsReal() < spx_->numColsReal()) {
    for (i = 0; i < spx_->numRowsReal(); ++i)
      num_non_zeros += spx_->rowVectorRealInternal(i).size();
  } else {
    for (i = 0; i < spx_->numColsReal(); ++i)
      num_non_zeros += spx_->colVectorRealInternal(i).size();
  }
  std::cout << " num_non_zeros " << num_non_zeros << std::endl;

  return num_non_zeros;
}

// gets the objective sense of the LP
LPObjectiveSense LPSoplexInterface::GetObjectiveSense() const {
  MiniMIPdebugMessage("calling GetObjectiveSense()\n");

  return (spx_->intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE)
             ? LPObjectiveSense::kMinimization
             : LPObjectiveSense::kMaximization;
}

// gets columns from LP problem object
//
// Either both, lower_bound and upper_bound, have to be 0, or both have to be
// non-0, either n_non_zeroes, begin_cols, indices, and obj_coeffs have to be 0,
// or all of them have to be non-0.
LPInterface::SparseVector LPSoplexInterface::GetSparseColumnCoefficients(
    int col) const {
  MiniMIPdebugMessage("calling GetSparseColumnCoefficients()\n");

  assert(0 <= col && col < spx_->numColsReal());

  LPInterface::SparseVector sparse_column;

  if (spx_->boolParam(SoPlex::PERSISTENTSCALING)) {
    DSVector cvec;
    spx_->getColVectorReal(col, cvec);
    for (int j = 0; j < cvec.size(); ++j) {
      sparse_column.indices.push_back(cvec.index(j));
      sparse_column.values.push_back(cvec.value(j));
    }
  } else {
    const SVector& cvec = spx_->colVectorRealInternal(col);
    for (int j = 0; j < cvec.size(); ++j) {
      sparse_column.indices.push_back(cvec.index(j));
      sparse_column.values.push_back(cvec.value(j));
    }
  }
  return sparse_column;
}

// gets rows from LP problem object
//
// Either both, left_hand_side and right_hand_side, have to be 0, or both have
// to be non-0, either n_non_zeroes, begin_cols, indices, and obj_coeffs have to
// be 0, or all of them have to be non-0.
LPInterface::SparseVector LPSoplexInterface::GetSparseRowCoefficients(
    int row) const {
  MiniMIPdebugMessage("calling GetSparseRowCoefficients()\n");

  assert(0 <= row && row < spx_->numRowsReal());

  LPInterface::SparseVector sparse_row;

  if (spx_->boolParam(SoPlex::PERSISTENTSCALING)) {
    DSVector rvec;
    spx_->getRowVectorReal(row, rvec);
    for (int j = 0; j < rvec.size(); ++j) {
      sparse_row.indices.push_back(rvec.index(j));
      sparse_row.values.push_back(rvec.value(j));
    }
  } else {
    const SVector& rvec = spx_->rowVectorRealInternal(row);
    for (int j = 0; j < rvec.size(); ++j) {
      sparse_row.indices.push_back(rvec.index(j));
      sparse_row.values.push_back(rvec.value(j));
    }
  }
  return sparse_row;
}

// gets objective coefficients from LP problem object
double LPSoplexInterface::GetObjectiveCoefficient(int col) const {
  MiniMIPdebugMessage("calling GetObjectiveCoefficient()\n");

  assert(0 <= col && col < spx_->numColsReal());

  return spx_->objReal(col);
}

// gets current lower bound of column from LP problem object
double LPSoplexInterface::GetLowerBound(int col) const {
  MiniMIPdebugMessage("calling GetLowerBound()\n");

  assert(0 <= col && col < spx_->numColsReal());

  return spx_->lowerReal(col);
}

// gets current upper bound of column from LP problem object
double LPSoplexInterface::GetUpperBound(int col) const {
  MiniMIPdebugMessage("calling GetUpperBound()\n");

  assert(0 <= col && col < spx_->numColsReal());

  return spx_->upperReal(col);
}

// gets current left hand sides of row from LP problem object
double LPSoplexInterface::GetLeftHandSide(int row) const {
  MiniMIPdebugMessage("calling GetLeftHandSide()\n");

  assert(0 <= row && row < spx_->numRowsReal());

  return spx_->lhsReal(row);
}

// gets current right hand sides of row from LP problem object
double LPSoplexInterface::GetRightHandSide(int row) const {
  MiniMIPdebugMessage("calling GetRightHandSide()\n");

  assert(0 <= row && row < spx_->numRowsReal());

  return spx_->rhsReal(row);
}

// gets the matrix coefficient of column and row from LP problem object
double LPSoplexInterface::GetMatrixCoefficient(
    int col,  // column number of coefficient
    int row   // row number of coefficient
) const {
  MiniMIPdebugMessage("calling GetMatrixCoefficient()\n");

  assert(0 <= col && col < spx_->numColsReal());
  assert(0 <= row && row < spx_->numRowsReal());

  return spx_->coefReal(row, col);
}

// ==========================================================================
// Solving methods.
// ==========================================================================

// calls primal simplex to solve the LP
absl::Status LPSoplexInterface::SolveLPWithPrimalSimplex() {
  MiniMIPdebugMessage("calling SolveLPWithPrimalSimplex()\n");

  static_cast<void>(
      spx_->setIntParam(SoPlex::ALGORITHM, SoPlex::ALGORITHM_PRIMAL));
  return SoPlexSolve();
}

// calls dual simplex to solve the LP
absl::Status LPSoplexInterface::SolveLPWithDualSimplex() {
  MiniMIPdebugMessage("calling SolveLPWithDualSimplex()\n");

  static_cast<void>(
      spx_->setIntParam(SoPlex::ALGORITHM, SoPlex::ALGORITHM_DUAL));
  return SoPlexSolve();
}

// start strong branching - call before any strong branching
absl::Status LPSoplexInterface::StartStrongBranching() {
  assert(PreStrongBranchingBasisFreed());
  SavePreStrongbranchingBasis();

  return absl::OkStatus();
}

// end strong branching - call after any strong branching
absl::Status LPSoplexInterface::EndStrongBranching() {
  assert(!PreStrongBranchingBasisFreed());
  RestorePreStrongbranchingBasis();
  FreePreStrongBranchingBasis();

  return absl::OkStatus();
}

// performs strong branching iterations on one branching candidate
absl::StatusOr<LPInterface::StrongBranchResult>
LPSoplexInterface::StrongBranchValue(
    int col,             // column to apply strong branching on
    double primal_sol,   // current primal solution value of column
    int iteration_limit  // iteration limit for strong branchings
) {
  StrongBranchResult strong_branch_result;
  // pass call on to StrongBranch()
  absl::Status absl_status_code = StrongBranch(
      col, primal_sol, iteration_limit,
      strong_branch_result.dual_bound_down_branch,
      strong_branch_result.dual_bound_up_branch,
      strong_branch_result.down_valid, strong_branch_result.up_valid,
      strong_branch_result.iterations);

  // pass absl::Status(absl::StatusCode::kInternal, "LP Error") to MiniMIP
  // without a back trace
  //@TODO: look into possible statuscode returns.
  if (absl_status_code != absl::OkStatus())
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  else
    return strong_branch_result;
}

// ============================================================================
// Solution information getters.
// ============================================================================

// returns whether a solve method was called after the last modification of the
// LP
bool LPSoplexInterface::IsSolved() const { return solved_; }

// returns true iff LP is proven to have a primal unbounded ray (but not
// necessary a primal feasible point); this does not necessarily meanthat the
// solver knows and can return the primal ray
bool LPSoplexInterface::ExistsPrimalRay() const {
  MiniMIPdebugMessage("calling ExistsPrimalRay()\n");

  return (spx_->status() == SPxSolver::UNBOUNDED);
}

// returns true iff LP is proven to have a primal unbounded ray (but not
// necessary a primal feasible point), and the solver knows and can return the
// primal ray
bool LPSoplexInterface::HasPrimalRay() const {
  MiniMIPdebugMessage("calling HasPrimalRay()\n");

  return spx_->hasPrimalRay();
}

// returns true iff LP is proven to be primal unbounded
bool LPSoplexInterface::IsPrimalUnbounded() const {
  MiniMIPdebugMessage("calling IsPrimalUnbounded()\n");

  assert(spx_->status() != SPxSolver::UNBOUNDED ||
         spx_->basisStatus() == SPxBasis::UNBOUNDED);

  // if SoPlex returns unbounded, this may only mean that an unbounded ray is
  // available, not necessarily a primal
  //* feasible point; hence we have to check the perturbation
  return spx_->status() == SPxSolver::UNBOUNDED;
}

// returns true iff LP is proven to be primal infeasible
bool LPSoplexInterface::IsPrimalInfeasible() const {
  MiniMIPdebugMessage("calling IsPrimalInfeasible()\n");

  return (spx_->status() == SPxSolver::INFEASIBLE);
}

// returns true iff LP is proven to be primal feasible
bool LPSoplexInterface::IsPrimalFeasible() const {
  MiniMIPdebugMessage("calling IsPrimalFeasible()\n");

  return spx_->basisStatus() == SPxBasis::OPTIMAL ||
         spx_->basisStatus() == SPxBasis::PRIMAL;
}

// returns true iff LP is proven to have a dual unbounded ray (but not necessary
// a dual feasible point);
//*  this does not necessarily meanthat the solver knows and can return the dual
// ray
bool LPSoplexInterface::ExistsDualRay() const {
  MiniMIPdebugMessage("calling ExistsDualRay()\n");

  return (spx_->status() == SPxSolver::INFEASIBLE);
}

// returns true iff LP is proven to have a dual unbounded ray (but not necessary
// a dual feasible point),
//*  and the solver knows and can return the dual ray
bool LPSoplexInterface::HasDualRay() const {
  MiniMIPdebugMessage("calling HasDualRay()\n");

  return spx_->hasDualFarkas();
}

// returns true iff LP is proven to be dual unbounded
bool LPSoplexInterface::IsDualUnbounded() const {
  MiniMIPdebugMessage("calling IsDualUnbounded()\n");

  return spx_->status() == SPxSolver::INFEASIBLE &&
         spx_->basisStatus() == SPxBasis::DUAL;
}

// returns true iff LP is proven to be dual infeasible
bool LPSoplexInterface::IsDualInfeasible() const {
  MiniMIPdebugMessage("calling IsDualInfeasible()\n");

  return (spx_->status() == SPxSolver::UNBOUNDED);
}

// returns true iff LP is proven to be dual feasible
bool LPSoplexInterface::IsDualFeasible() const {
  MiniMIPdebugMessage("calling IsDualFeasible()\n");

  return (spx_->basisStatus() == SPxBasis::OPTIMAL) ||
         spx_->basisStatus() == SPxBasis::DUAL;
}

// returns true iff LP was solved to optimality
bool LPSoplexInterface::IsOptimal() const {
  MiniMIPdebugMessage("calling IsOptimal()\n");

  assert((spx_->basisStatus() == SPxBasis::OPTIMAL) ==
         (IsPrimalFeasible() && IsDualFeasible()));

  return (spx_->status() == SPxSolver::OPTIMAL);
}

// returns true iff current LP solution is stable
//*
//*  This function should return true if the solution is reliable, i.e.,
// feasible and optimal (or proven
//*  infeasible/unbounded) with respect to the original problem. The optimality
// status might be with respect to a scaled
//*  version of the problem, but the solution might not be feasible to the
// unscaled original problem; in this case,
//*  minimip::LPInterface.IsStable() should return false.
bool LPSoplexInterface::IsStable() const {
  MiniMIPdebugMessage("calling IsStable()\n");

  if (spx_->status() == SPxSolver::ERROR ||
      spx_->status() == SPxSolver::SINGULAR)
    return false;
  if (spx_->status() == SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS) return false;
  return true;
}

// returns true iff the objective limit was reached
bool LPSoplexInterface::ObjectiveLimitIsExceeded() const {
  MiniMIPdebugMessage("calling ObjectiveLimitIsExceeded()\n");

  return (spx_->status() == SPxSolver::ABORT_VALUE);
}

// returns true iff the iteration limit was reached
bool LPSoplexInterface::IterationLimitIsExceeded() const {
  MiniMIPdebugMessage("calling IterationLimitIsExceeded()\n");

  return (spx_->status() == SPxSolver::ABORT_ITER);
}

// returns true iff the time limit was reached
bool LPSoplexInterface::TimeLimitIsExceeded() const {
  MiniMIPdebugMessage("calling TimeLimitIsExceeded()\n");

  return (spx_->status() == SPxSolver::ABORT_TIME);
}

// gets objective value of solution
double LPSoplexInterface::GetObjectiveValue() {
  MiniMIPdebugMessage("calling GetObjectiveValue()\n");

  return spx_->objValueReal();
}

// gets primal and dual solution vectors for feasible LPs
//*
//*  Before calling this function, the caller must ensure that the LP has been
// solved to optimality, i.e., that
//*  minimip::LPInterface.IsOptimal() returns true.

absl::StatusOr<std::vector<double>> LPSoplexInterface::GetPrimalSolution()
    const {
  MiniMIPdebugMessage("calling GetPrimalSolution()\n");
  std::vector<double> primal_sol;

  try {
    static_cast<void>(
        spx_->getPrimalReal(primal_sol.data(), spx_->numColsReal()));
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    // MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
    // %s\n", s.c_str());
#else
  catch (const SPxException&) {
#endif
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }

  return primal_sol;
}

absl::StatusOr<std::vector<double>> LPSoplexInterface::GetDualSolution() const {
  MiniMIPdebugMessage("calling GetDualSolution()\n");
  std::vector<double> dual_sol;

  try {
    static_cast<void>(spx_->getDualReal(dual_sol.data(), spx_->numRowsReal()));
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    // MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
    // %s\n", s.c_str());
#else
  catch (const SPxException&) {
#endif
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }

  return dual_sol;
}

absl::StatusOr<std::vector<double>> LPSoplexInterface::GetRowActivity() const {
  MiniMIPdebugMessage("calling GetRowActivity()\n");
  std::vector<double> activity;

  try {
    static_cast<void>(spx_->getSlacksReal(
        activity.data(),
        spx_->numRowsReal()));  // in SoPlex, the activities are called "slacks"
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    // MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
    // %s\n", s.c_str());
#else
  catch (const SPxException&) {
#endif
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }

  return activity;
}

absl::StatusOr<std::vector<double>> LPSoplexInterface::GetReducedCost() const {
  MiniMIPdebugMessage("calling GetReducedCost()\n");
  std::vector<double> reduced_cost;
  try {
    static_cast<void>(
        spx_->getRedCostReal(reduced_cost.data(), spx_->numColsReal()));
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    // MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
    // %s\n", s.c_str());
#else
  catch (const SPxException&) {
#endif
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }

  return reduced_cost;
}

// gets primal ray for unbounded LPs
absl::StatusOr<std::vector<double>> LPSoplexInterface::GetPrimalRay() const {
  MiniMIPdebugMessage("calling GetPrimalRay()\n");
  std::vector<double> primal_ray;

  assert(spx_->hasPrimalRay());

  try {
    static_cast<void>(
        spx_->getPrimalRayReal(primal_ray.data(), spx_->numColsReal()));
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    // MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
    // %s\n", s.c_str());
#else
  catch (const SPxException&) {
#endif
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }
  return primal_ray;
}

// gets dual Farkas proof for infeasibility
absl::StatusOr<std::vector<double>> LPSoplexInterface::GetDualFarkasMultiplier()
    const {
  MiniMIPdebugMessage("calling GetDualFarkasMultiplier()\n");
  std::vector<double> dual_farkas_multiplier;

  assert(spx_->hasDualFarkas());

  try {
    static_cast<void>(spx_->getDualFarkasReal(dual_farkas_multiplier.data(),
                                              spx_->numRowsReal()));
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    // MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
    // %s\n", s.c_str());
#else
  catch (const SPxException&) {
#endif
    return absl::Status(absl::StatusCode::kInternal, "LP Error");
  }
  return dual_farkas_multiplier;
}

// gets the number of LP iterations of the last solve call
int LPSoplexInterface::GetIterations() const {
  MiniMIPdebugMessage("calling GetIterations()\n");

  return spx_->numIterations();
}

// ==========================================================================
// Getters and setters of the basis.
// ==========================================================================

// gets current basis status for columns and rows
absl::StatusOr<std::vector<LPBasisStatus>>
LPSoplexInterface::GetColumnBasisStatus() const {
  std::vector<LPBasisStatus> column_basis_status;
  MiniMIPdebugMessage("calling GetColumnBasisStatus()\n");
  assert(PreStrongBranchingBasisFreed());

  if (column_basis_status.size() != 0) {
    for (int i = 0; i < spx_->numColsReal(); ++i) {
      //         double obj_coeffs = 0.0;
      switch (spx_->basisColStatus(i)) {
        case SPxSolver::BASIC:
          column_basis_status.push_back(LPBasisStatus::kBasic);
          break;
        case SPxSolver::FIXED:
          // Get reduced cost estimation. If the estimation is not correct this
          // should not hurt: If the basis is loaded into SoPlex again, the
          // status is converted to FIXED again; in this case there is no
          // problem at all. If the basis is saved and/or used in some other
          // solver, it usually is very cheap to perform the pivots necessary to
          // get an optimal basis.
          // @todo implement getRedCostEst()
          //
          // MINIMIP_CALL( getRedCostEst(spx_, i, &obj_coeffs) );
          //             if( obj_coeffs < 0.0 )  // reduced costs < 0 => UPPER
          //             else => LOWER
          //                column_basis_status.push_back(BASESTAT_UPPER;
          //             else
          column_basis_status.push_back(LPBasisStatus::kAtLowerBound);
          break;
        case SPxSolver::ON_LOWER:
          column_basis_status.push_back(LPBasisStatus::kAtLowerBound);
          break;
        case SPxSolver::ON_UPPER:
          column_basis_status.push_back(LPBasisStatus::kAtUpperBound);
          break;
        case SPxSolver::ZERO:
          column_basis_status.push_back(LPBasisStatus::kFree);
          break;
        case SPxSolver::UNDEFINED:
        default:
          MiniMIPerrorMessage("invalid basis status\n");
          // MINIMIP_ABORT();
          assert(false);
          return absl::Status(absl::StatusCode::kInvalidArgument,
                              "Invalid Data");
      }
    }
  }
  return column_basis_status;
}

absl::StatusOr<std::vector<LPBasisStatus>>
LPSoplexInterface::GetRowBasisStatus() const {
  std::vector<LPBasisStatus> row_basis_status;
  MiniMIPdebugMessage("calling GetRowBasisStatus()\n");
  assert(PreStrongBranchingBasisFreed());

  if (row_basis_status.size() != 0) {
    for (int i = 0; i < spx_->numRowsReal(); ++i) {
      switch (spx_->basisRowStatus(i)) {
        case SPxSolver::BASIC:
          row_basis_status.push_back(LPBasisStatus::kBasic);
          break;
        case SPxSolver::FIXED:
          row_basis_status.push_back(LPBasisStatus::kFixed);
        case SPxSolver::ON_LOWER:
          row_basis_status.push_back(LPBasisStatus::kAtLowerBound);
          break;
        case SPxSolver::ON_UPPER:
          row_basis_status.push_back(LPBasisStatus::kAtUpperBound);
          break;
        case SPxSolver::ZERO:
          MiniMIPerrorMessage(
              "slack variable has basis status ZERO (should not occur)\n");
          return absl::Status(absl::StatusCode::kInternal, "LP Error");
        case SPxSolver::UNDEFINED:
        default:
          MiniMIPerrorMessage("invalid basis status\n");
          // std::abort();
          return absl::Status(absl::StatusCode::kInvalidArgument,
                              "Invalid Data");
      }
    }
  }
  return row_basis_status;
}

// sets current basis status for columns and rows
absl::Status LPSoplexInterface::SetBase(
    const std::vector<LPBasisStatus>&
        column_basis_status,  // array with column basis status
    const std::vector<LPBasisStatus>&
        row_basis_status  // array with row basis status
) {
  MiniMIPdebugMessage("calling SetBase()\n");

  int num_cols = GetNumberOfColumns();
  int num_rows = GetNumberOfRows();

  assert(PreStrongBranchingBasisFreed());
  InvalidateSolution();

  DataArray<SPxSolver::VarStatus>& _colstat = ColumnsBasisStatus();
  DataArray<SPxSolver::VarStatus>& _rowstat = RowsBasisStatus();

  _colstat.reSize(num_cols);
  _rowstat.reSize(num_rows);

  for (int i = 0; i < num_rows; ++i) {
    switch (row_basis_status[i]) {
      case LPBasisStatus::kBasic:
        _rowstat[i] = SPxSolver::BASIC;
        break;
      case LPBasisStatus::kAtLowerBound:
        _rowstat[i] = SPxSolver::ON_LOWER;
        break;
      case LPBasisStatus::kAtUpperBound:
        _rowstat[i] = SPxSolver::ON_UPPER;
        break;
      case LPBasisStatus::kFixed:
        _rowstat[i] = SPxSolver::FIXED;
        break;
      case LPBasisStatus::kFree:
        MiniMIPerrorMessage(
            "slack variable has basis status ZERO (should not occur)\n");
        return absl::Status(absl::StatusCode::kInternal, "LP Error");
      default:
        MiniMIPerrorMessage("invalid basis status\n");

        // MINIMIP_ABORT();
        assert(false);

        return absl::Status(absl::StatusCode::kInvalidArgument, "Invalid Data");
    }
  }

  for (int i = 0; i < num_cols; ++i) {
    switch (column_basis_status[i]) {
      case LPBasisStatus::kBasic:
        _colstat[i] = SPxSolver::BASIC;
        break;
      case LPBasisStatus::kAtLowerBound:
        _colstat[i] = SPxSolver::ON_LOWER;
        break;
      case LPBasisStatus::kAtUpperBound:
        _colstat[i] = SPxSolver::ON_UPPER;
        break;
      case LPBasisStatus::kFixed:
        _colstat[i] = SPxSolver::FIXED;
        break;
      case LPBasisStatus::kFree:
        _colstat[i] = SPxSolver::ZERO;
        break;
      default:
        MiniMIPerrorMessage("invalid basis status\n");
        // MINIMIP_ABORT();
        assert(false);
        return absl::Status(absl::StatusCode::kInvalidArgument, "Invalid Data");
    }
  }

  SOPLEX_TRY(spx_->setBasis(_rowstat.get_ptr(), _colstat.get_ptr()));
  FreePreStrongBranchingBasis();

  return absl::OkStatus();
}

// returns the indices of the basic columns and rows; basic column n gives value
// n, basic row m gives value -1-m
std::vector<int> LPSoplexInterface::GetBasisIndices() const {
  std::vector<int> basis_indices;
  MiniMIPdebugMessage("calling GetBasisInd()\n");

  assert(PreStrongBranchingBasisFreed());

  spx_->getBasisInd(basis_indices.data());

  return basis_indices;
}

// ==========================================================================
// Getters of vectors in the inverted basis matrix.
// ==========================================================================

// get row of inverse basis matrix B^-1
//
// NOTE: The LP interface defines slack variables to have coefficient +1. This
// means that if, internally, the LP solver
//       uses a -1 coefficient, then rows associated with slacks variables whose
//       coefficient is -1, should be negated; see also the explanation in
//       lpi.h.
absl::StatusOr<LPInterface::SparseVector>
LPSoplexInterface::GetSparseRowOfBInverted(int row_number) const {
  LPInterface::SparseVector sparse_row;
  int num_indices;
  MiniMIPdebugMessage("calling GetSparseRowOfBInverted()\n");

  assert(PreStrongBranchingBasisFreed());
  assert(row_number >= 0);
  assert(row_number < spx_->numRowsReal());

  std::vector<int> integer_indices(sparse_row.indices.begin(),
                                   sparse_row.indices.end());

  if (!spx_->getBasisInverseRowReal(row_number, sparse_row.values.data(),
                                    integer_indices.data(), &num_indices))
    return absl::Status(absl::StatusCode::kInternal, "LP Error");

  sparse_row.indices.assign(integer_indices.begin(), integer_indices.end());

  return sparse_row;
}

// get column of inverse basis matrix B^-1
//
// NOTE: The LP interface defines slack variables to have coefficient +1. This
// means that if, internally, the LP solver
//       uses a -1 coefficient, then rows associated with slacks variables whose
//       coefficient is -1, should be negated; see also the explanation in
//       lpi.h.
//
// column number of B^-1; this is NOT the number of the
// column in the LP; you have to call
// minimip::LPInterface.GetBasisIndices() to get the array
// which links the B^-1 column numbers to the row and
// column numbers of the LP! c must be between 0 and
// num_rows-1, since the basis has the size num_rows *
// num_rows

absl::StatusOr<LPInterface::SparseVector>
LPSoplexInterface::GetSparseColumnOfBInverted(int col_number) const {
  MiniMIPdebugMessage("calling GetColumnOfBInverted()\n");

  LPInterface::SparseVector sparse_column;
  int num_indices;

  assert(PreStrongBranchingBasisFreed());

  std::vector<int> integer_indices;

  if (!spx_->getBasisInverseColReal(col_number, sparse_column.values.data(),
                                    integer_indices.data(), &num_indices))
    return absl::Status(absl::StatusCode::kInternal, "LP Error");

  sparse_column.indices.assign(integer_indices.begin(), integer_indices.end());

  return sparse_column;
}

// get row of inverse basis matrix times constraint matrix B^-1 * A
//
// NOTE: The LP interface defines slack variables to have coefficient +1. This
// means that if, internally, the LP solver
//       uses a -1 coefficient, then rows associated with slacks variables whose
//       coefficient is -1, should be negated; see also the explanation in
//       lpi.h.
absl::StatusOr<LPInterface::SparseVector>
LPSoplexInterface::GetSparseRowOfBInvertedTimesA(int row_number) const {
  LPInterface::SparseVector sparse_row;
  std::vector<double> buf;
  std::vector<double> binv;
  size_t num_rows;
  int num_cols;
  int c;

  MiniMIPdebugMessage("calling GetSparseRowOfBInvertedTimesA()\n");

  assert(PreStrongBranchingBasisFreed());

  num_rows = spx_->numRowsReal();
  num_cols = spx_->numColsReal();

  buf.resize(num_rows);

  // calculate the row in B^-1
  absl::StatusOr<LPInterface::SparseVector> absl_tmp_row =
      GetSparseRowOfBInverted(row_number);
  if (absl_tmp_row.ok()) {
    buf                = absl_tmp_row->values;
    sparse_row.indices = absl_tmp_row->indices;
  } else {
    std::cout << absl_tmp_row.status();
  }

  binv.assign(buf.begin(), buf.end());

  assert(binv.size() == num_rows);

  // @todo exploit sparsity in binv by looping over num_rows
  // calculate the scalar product of the row in B^-1 and A
  Vector binv_vec(num_rows, binv.data());

  // temporary unscaled column of A
  DSVector acol;

  for (c = 0; c < num_cols; ++c) {
    spx_->getColVectorReal(c, acol);
    sparse_row.values.push_back(binv_vec * acol);  // scalar product
  }

  return sparse_row;
}

// get column of inverse basis matrix times constraint matrix B^-1 * A
//
// NOTE: The LP interface defines slack variables to have coefficient +1. This
// means that if, internally, the LP solver
//       uses a -1 coefficient, then rows associated with slacks variables whose
//       coefficient is -1, should be negated; see also the explanation in
//       lpi.h.
absl::StatusOr<LPInterface::SparseVector>
LPSoplexInterface::GetSparseColumnOfBInvertedTimesA(int col_number) const {
  LPInterface::SparseVector sparse_row;
  // create a new uninitialized full vector
  DVector col(spx_->numRowsReal());

  // temporary sparse vector used for unscaling (memory is automatically
  // enlarged)
  DSVector colsparse;

  MiniMIPdebugMessage("calling GetColumnOfBInvertedTimesA()\n");

  assert(PreStrongBranchingBasisFreed());

  // extract column col_number of A
  assert(col_number >= 0);
  assert(col_number < spx_->numColsReal());

  // @todo implement this with sparse vectors

  // col needs to be cleared because copying colVectorReal only regards nonzeros
  col.clear();

  spx_->getColVectorReal(col_number, colsparse);
  // the copy is necessary to transform the sparse column into a dense vector
  col = colsparse;

  // solve
  if (!spx_->getBasisInverseTimesVecReal(col.get_ptr(),
                                         sparse_row.values.data()))
    return absl::Status(absl::StatusCode::kInternal, "LP Error");

  return sparse_row;
}

// ==========================================================================
// Getters and setters of the parameters.
// ==========================================================================

// gets integer parameter of LP
absl::StatusOr<int> LPSoplexInterface::GetIntegerParameter(
    LPParameter type  // parameter number
) const {
  int scale_param;
  int param_val;

  MiniMIPdebugMessage("calling GetIntegerParameter()\n");

  switch (type) {
    case LPParameter::kFromScratch:
      param_val = GetFromScratch();
      break;
    case LPParameter::kLPInfo:
      param_val = GetLPInfo();
      break;
    case LPParameter::kLPIterationLimit:
      if (spx_->intParam(SoPlex::ITERLIMIT) == -1) param_val = INT_MAX;
      break;
      param_val = spx_->intParam(SoPlex::ITERLIMIT);
    case LPParameter::kPresolving:
      param_val = spx_->intParam(SoPlex::SIMPLIFIER) == SoPlex::SIMPLIFIER_AUTO;
      break;
    case LPParameter::kPricing:
      param_val = (int)pricing_;
      break;
    case LPParameter::kScaling:
      scale_param = spx_->intParam(SoPlex::SCALER);

      if (scale_param == SoPlex::SCALER_OFF)
        param_val = 0;
      else if (scale_param == SoPlex::SCALER_BIEQUI)
        param_val = 1;
      else {
        assert(scale_param == SoPlex::SCALER_LEASTSQ);
        param_val = 2;
      }
      break;
    case LPParameter::kTiming:
      param_val = static_cast<int>(spx_->intParam(SoPlex::TIMER));
      break;
    case LPParameter::kRandomSeed:
      param_val = static_cast<int>(spx_->randomSeed());
      break;
    case LPParameter::kRefactor:
      param_val = static_cast<int>(spx_->intParam(SoPlex::FACTOR_UPDATE_MAX));
      break;
    default:
      return absl::Status(absl::StatusCode::kInvalidArgument,
                          "Parameter Unknown");
  }

  return param_val;
}

// sets integer parameter of LP
absl::Status LPSoplexInterface::SetIntegerParameter(
    LPParameter type,  // parameter number
    int param_val      // parameter value
) {
  MiniMIPdebugMessage("calling SetIntegerParameter()\n");

  switch (type) {
    case LPParameter::kFromScratch:
      assert(param_val == true || param_val == false);
      SetFromScratch(static_cast<bool>(param_val));
      break;
    case LPParameter::kLPInfo:
      assert(param_val == true || param_val == false);
      SetLPInfo(static_cast<bool>(param_val));
      break;
    case LPParameter::kLPIterationLimit:
      assert(param_val >= 0);
      // -1 <= param_val, -1 meaning no time limit, 0 stopping immediately
      if (param_val >= INT_MAX)
        static_cast<void>(spx_->setIntParam(SoPlex::ITERLIMIT, -1));
      break;
    case LPParameter::kPresolving:
      assert(param_val == true || param_val == false);
      static_cast<void>(spx_->setIntParam(
          SoPlex::SIMPLIFIER,
          (param_val ? SoPlex::SIMPLIFIER_AUTO : SoPlex::SIMPLIFIER_OFF)));
      break;
    case LPParameter::kPricing:
      pricing_ = (LPPricing)param_val;
      switch (pricing_) {
        case LPPricing::kDefault:
        case LPPricing::kAuto:
          static_cast<void>(
              spx_->setIntParam(SoPlex::PRICER, SoPlex::PRICER_AUTO));
          break;
        case LPPricing::kFull:
          static_cast<void>(
              spx_->setIntParam(SoPlex::PRICER, SoPlex::PRICER_STEEP));
          break;
        case LPPricing::kPartial:
          static_cast<void>(
              spx_->setIntParam(SoPlex::PRICER, SoPlex::PRICER_PARMULT));
          break;
        case LPPricing::kSteep:
          static_cast<void>(
              spx_->setIntParam(SoPlex::PRICER, SoPlex::PRICER_STEEP));
          break;
        case LPPricing::kSteepQStart:
          static_cast<void>(
              spx_->setIntParam(SoPlex::PRICER, SoPlex::PRICER_QUICKSTEEP));
          break;
        case LPPricing::kDevex:
          static_cast<void>(
              spx_->setIntParam(SoPlex::PRICER, SoPlex::PRICER_DEVEX));
          break;
        default:
          return absl::Status(absl::StatusCode::kInternal, "LP Error");
      }
      break;
    case LPParameter::kScaling:
      assert(param_val >= 0 && param_val <= 2);
      if (param_val == 0)
        static_cast<void>(
            spx_->setIntParam(SoPlex::SCALER, SoPlex::SCALER_OFF));
      else if (param_val == 1)
        static_cast<void>(
            spx_->setIntParam(SoPlex::SCALER, SoPlex::SCALER_BIEQUI));
      else
        static_cast<void>(
            spx_->setIntParam(SoPlex::SCALER, SoPlex::SCALER_LEASTSQ));
      break;
    case LPParameter::kTiming:
      assert(param_val >= 0 && param_val < 3);
      static_cast<void>(spx_->setIntParam(SoPlex::TIMER, param_val));
      break;
    case LPParameter::kRandomSeed:
      spx_->setRandomSeed(
          static_cast<unsigned long>(static_cast<long>(param_val)));
      break;
    case LPParameter::kRefactor:
      assert(param_val >= 0);
      static_cast<void>(
          spx_->setIntParam(SoPlex::FACTOR_UPDATE_MAX, param_val));
      break;
    default:
      return absl::Status(absl::StatusCode::kInvalidArgument,
                          "Parameter Unknown");
  }

  return absl::OkStatus();
}

// gets floating point parameter of LP
absl::StatusOr<double> LPSoplexInterface::GetRealParameter(
    LPParameter type  // parameter number
) const {
  double param_val;
  MiniMIPdebugMessage("calling GetRealParameter()\n");

  switch (type) {
    case LPParameter::kFeasibilityTolerance:
      param_val = FeasibilityTolerance();
      break;
    case LPParameter::kDualFeasibilityTolerance:
      param_val = OptimalityTolerance();
      break;
    case LPParameter::kObjectiveLimit:
      if (spx_->intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE)
        param_val = spx_->realParam(SoPlex::OBJLIMIT_UPPER);
      else
        param_val = spx_->realParam(SoPlex::OBJLIMIT_LOWER);
      break;
    case LPParameter::kLPTimeLimit:
      param_val = spx_->realParam(SoPlex::TIMELIMIT);
      break;
    case LPParameter::kMarkowitz:
      param_val = spx_->realParam(SoPlex::MIN_MARKOWITZ);
      break;
    default:
      return absl::Status(absl::StatusCode::kInvalidArgument,
                          "Parameter Unknown");
  }

  return param_val;
}

// sets floating point parameter of LP
absl::Status LPSoplexInterface::SetRealParameter(
    LPParameter type,  // parameter number
    double param_val   // parameter value
) {
  MiniMIPdebugMessage("calling SetRealParameter()\n");

  switch (type) {
    case LPParameter::kFeasibilityTolerance:
      // 0 < param_val
      assert(param_val > 0.0);
      SetFeasibilityTolerance(param_val);
      break;
    case LPParameter::kDualFeasibilityTolerance:
      // 0 < param_val
      assert(param_val > 0.0);
      SetOptimalityTolerance(param_val);
      break;
    case LPParameter::kObjectiveLimit:
      // no restrictions on param_val
      if (spx_->intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE)
        static_cast<void>(
            spx_->setRealParam(SoPlex::OBJLIMIT_UPPER, param_val));
      else
        static_cast<void>(
            spx_->setRealParam(SoPlex::OBJLIMIT_LOWER, param_val));
      break;
    case LPParameter::kLPTimeLimit:
      assert(param_val > 0.0);
      // soplex requires 0 < param_val < DEFAULT_INFINITY (= 1e100), -1 means
      // unlimited
      static_cast<void>(spx_->setRealParam(SoPlex::TIMELIMIT, param_val));
      break;
    case LPParameter::kMarkowitz:
      // 1e-4 <= param_val <= 0.999
      if (param_val < 1e-4)
        param_val = 1e-4;
      else if (param_val > 0.9999)
        param_val = 0.9999;

      static_cast<void>(spx_->setRealParam(SoPlex::MIN_MARKOWITZ, param_val));
      break;
    default:
      return absl::Status(absl::StatusCode::kInvalidArgument,
                          "Parameter Unknown");
  }

  return absl::OkStatus();
}

// ==========================================================================
// Numerical methods.
// ==========================================================================

// returns value treated as infinity in the LP solver
double LPSoplexInterface::Infinity() const {
  MiniMIPdebugMessage("calling Infinity()\n");

  return spx_->realParam(SoPlex::INFTY);
}

// checks if given value is treated as infinity in the LP solver
bool LPSoplexInterface::IsInfinity(
    double value  // value to be checked for infinity
) const {
  MiniMIPdebugMessage("calling IsInfinity()\n");

  return (value >= spx_->realParam(SoPlex::INFTY));
}

// returns, whether the given file exists
bool LPSoplexInterface::FileExists(const char* file_name  // file name
) const {
  FILE* file;

  file = fopen(file_name, "r");
  if (file == NULL) return false;

  fclose(file);

  return true;
}

// ==========================================================================
// File interface methods.
// ==========================================================================

// reads LP from a file
absl::Status LPSoplexInterface::ReadLP(const char* file_name  // file name
) {
  MiniMIPdebugMessage("calling ReadLP()\n");

  assert(file_name != NULL);

  assert(PreStrongBranchingBasisFreed());

  if (!FileExists(file_name))
    return absl::Status(absl::StatusCode::kInternal, "Read Errror");

  try {
    assert(spx_->intParam(SoPlex::READMODE) == SoPlex::READMODE_REAL);
    if (!spx_->readFile((char*)(file_name)))
      return absl::Status(absl::StatusCode::kInternal, "Read Errror");
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    // MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
    // %s\n", s.c_str());
#else
  catch (const SPxException&) {
#endif
    return absl::Status(absl::StatusCode::kInternal, "Read Errror");
  }

  return absl::OkStatus();
}

// writes LP to a file
absl::Status LPSoplexInterface::WriteLP(const char* file_name  // file name
) const {
  MiniMIPdebugMessage("calling WriteLP()\n");

  assert(file_name != NULL);

  try {
    static_cast<void>(spx_->writeFileReal(file_name));
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    // MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception:
    // %s\n", s.c_str());
#else
  catch (const SPxException&) {
#endif
    return absl::Status(absl::StatusCode::kInternal, "Write Errror");
  }
  return absl::OkStatus();
}

}  // namespace minimip
