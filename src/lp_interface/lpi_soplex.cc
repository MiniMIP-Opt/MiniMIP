
#include "src/lp_interface/lpi_soplex.h"

/* check the return value of setParam methods */
#define CHECK_SOPLEX_PARAM(x)                                            \
  if (!x) {                                                              \
    MiniMIPmessagePrintWarning("SoPlex: unsupported parameter value\n"); \
  }

/* remember the original value of the MiniMIP_DEBUG define and undefine it */
#ifdef MINIMIP_DEBUG
#define ___DEBUG
#undef MINIMIP_DEBUG
#endif

/* reset the MiniMIP_DEBUG define to its original MiniMIP value */
#undef MINIMIP_DEBUG
#ifdef ___DEBUG
#define MINIMIP_DEBUG
#undef ___DEBUG
#endif

#define SOPLEX_VERBLEVEL 0 /**< verbosity level for LPINFO */

/** Macro for a single SoPlex call for which exceptions have to be catched - return an LP error. We
 *  make no distinction between different exception types, e.g., between memory allocation and other
 *  exceptions.
 */
#ifndef NDEBUG
#define SOPLEX_TRY(x)                                                          \
  do {                                                                         \
    try {                                                                      \
      (x);                                                                     \
    } catch (const SPxMemoryException& E) {                                    \
      std::string s = E.what();                                                \
      MiniMIPerrorMessage("SoPlex threw a memory exception: %s\n", s.c_str()); \
      return RetCode::kError;                                                   \
    } catch (const SPxException& E) {                                          \
      std::string s = E.what();                                                \
      MiniMIPerrorMessage("SoPlex threw an exception: %s\n", s.c_str());       \
      return RetCode::kLPError;                                                \
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
      return ERROR;                                                            \
    } catch (const SPxException&) {                                            \
      return LP_ERROR;                                                         \
    }                                                                          \
  } while (FALSE)
#endif

/* Macro for a single SoPlex call for which exceptions have to be catched - abort if they
 * arise. MINIMIP_ABORT() is not accessible here.
 */
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

/** constructor */
LPSoplexInterface::LPSoplexInterface() : spx_(new soplex::SoPlex),
                                         pricing_(LPPricing::kDefault),
                                         solved_(false),
                                         lp_info_(false),
                                         from_scratch_(false),
                                         col_basis_status_(0),
                                         row_basis_status_(0) {}

/** destructor */
LPSoplexInterface::~LPSoplexInterface() {
  delete spx_;
}

/** provides access for temporary storage of basis status of rows */
DataArray<SPxSolver::VarStatus>& LPSoplexInterface::RowsBasisStatus() {
  return row_basis_status_; /*lint !e1536*/
}
/** provides access for temporary storage of basis status of columns */
DataArray<SPxSolver::VarStatus>& LPSoplexInterface::ColumnsBasisStatus() {
  return col_basis_status_; /*lint !e1536*/
}

void LPSoplexInterface::SetFromScratch(bool from_scratch) {
  from_scratch_ = from_scratch;
}

void LPSoplexInterface::SetLPInfo(bool lp_info) {
  lp_info_ = lp_info;
}

/** get objective limit according to objective sense */
LPValue LPSoplexInterface::GetObjectiveLimit() const {
  return (spx_->intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE)
    ? spx_->realParam(SoPlex::OBJLIMIT_UPPER)
    : spx_->realParam(SoPlex::OBJLIMIT_LOWER);
}

LPValue LPSoplexInterface::FeasibilityTolerance() const {
  return spx_->realParam(SoPlex::FEASTOL);
}

/** return optimality tolerance */
LPValue LPSoplexInterface::OptimalityTolerance() const {
  return spx_->realParam(SoPlex::OPTTOL);
}

/** set feasibility tolerance and store value in case SoPlex only accepts a larger tolerance */
void LPSoplexInterface::SetFeasibilityTolerance(
  const Real d) {
  assert(spx_->setRealParam(SoPlex::FEASTOL, d));
}

/** set optimality tolerance and store value in case SoPlex only accepts a larger tolerance */
void LPSoplexInterface::SetOptimalityTolerance(
  const Real d) {
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
      /*  MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception: %s\n", s.c_str()); */
    }
    /* since it is not clear if the status in SoPlex are set correctly we want to make sure that if an error is thrown
     * the status is not OPTIMAL anymore. */
    assert(spx_->status() != SPxSolver::OPTIMAL);
  }
  assert(spx_->intParam(SoPlex::ITERLIMIT) < 0 || spx_->numIterations() <= spx_->intParam(SoPlex::ITERLIMIT));

  /* update time limit */
  time_spent = spx_->solveTime();

  /* get current time limit */
  if (time_spent > 0) {
    time_limit = spx_->realParam(SoPlex::TIMELIMIT);
    if (time_limit > time_spent)
      time_limit -= time_spent;
    else
      time_limit = 0;

    /* set new time limit */
    assert(time_limit >= 0);

    assert(spx_->setRealParam(SoPlex::TIMELIMIT, time_limit));
  }
}

bool LPSoplexInterface::CheckConsistentBounds() const {
  for (int i = 0; i < spx_->numColsReal(); ++i) {
    if (spx_->lowerReal(i) > spx_->upperReal(i) + spx_->realParam(SoPlex::EPSILON_ZERO)) {
      MiniMIPerrorMessage("inconsistent bounds on column %d: lower=%.17g, upper=%.17g\n",
                          i, spx_->lowerReal(i), spx_->upperReal(i));
      return false;
    }
  }
  return true;
}

bool LPSoplexInterface::CheckConsistentSides() const {
  for (int i = 0; i < spx_->numRowsReal(); ++i) {
    if (spx_->lhsReal(i) > spx_->rhsReal(i) + spx_->realParam(SoPlex::EPSILON_ZERO)) {
      MiniMIPerrorMessage("inconsistent sides on row %d: lhs=%.17g, rhs=%.17g\n",
                          i, spx_->lhsReal(i), spx_->rhsReal(i));
      return false;
    }
  }
  return true;
}

SPxSolver::Status LPSoplexInterface::doSolve(bool print_warning = true) {
  SPxOut::Verbosity verbosity;
  SPxSolver::Status spx_status;

  /* store and set verbosity */
  verbosity = spx_->spxout.getVerbosity();
  spx_->spxout.setVerbosity((SPxOut::Verbosity)(GetLPInfo() ? SOPLEX_VERBLEVEL : 0));

  assert(CheckConsistentBounds());
  assert(CheckConsistentSides());

  TrySolve(print_warning);
  spx_status = spx_->status();

  /* restore verbosity */
  spx_->spxout.setVerbosity(verbosity);

  return spx_status;
}

bool LPSoplexInterface::GetFromScratch() const {
  return from_scratch_;
}

/* @todo member variable? */
bool LPSoplexInterface::GetLPInfo() const {
  return lp_info_;
}

RetCode LPSoplexInterface::SoPlexSolve() {
  /* store and set verbosity */
  SPxOut::Verbosity verbosity = spx_->spxout.getVerbosity();
  spx_->spxout.setVerbosity((SPxOut::Verbosity)(GetLPInfo() ? SOPLEX_VERBLEVEL : 0));

  MiniMIPdebugMessage("calling SoPlex solve(): %d cols, %d rows\n", spx_->numColsReal(), spx_->numRowsReal());

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

  /* delete starting basis if solving from scratch */
  if (GetFromScratch()) {
    try {
      spx_->clearBasis();
    }
#ifndef NDEBUG
    catch (const SPxException& x) {
      std::string s = x.what();
      /* MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception: %s\n", s.c_str()); */
#else
    catch (const SPxException&) {
#endif
      assert(spx_->status() != SPxSolver::OPTIMAL);
      return RetCode::kLPError;
    }
  }
  assert(!GetFromScratch() || spx_->status() == SPxSolver::NO_PROBLEM);

  SPxSolver::Status status = doSolve();
  MiniMIPdebugMessage(" -> SoPlex status: %d, basis status: %d\n", spx_->status(), spx_->basisStatus());
  solved_ = true;

  /* restore verbosity */
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
      return RetCode::kOkay;
    default:
      return RetCode::kLPError;
  } /*lint !e788*/
}

/** save the current basis */
void LPSoplexInterface::SavePreStrongbranchingBasis() {
  row_basis_status_.reSize(spx_->numRowsReal());
  col_basis_status_.reSize(spx_->numColsReal());

  try {
    spx_->getBasis(row_basis_status_.get_ptr(), col_basis_status_.get_ptr());
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    /* MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception: %s\n", s.c_str()); */

    /* since it is not clear if the status in SoPlex are set correctly we want to make sure that if an error is thrown
     * the status is not OPTIMAL anymore. */
    assert(spx_->status() != SPxSolver::OPTIMAL);
  }
#else
  catch (const SPxException&) {
  }
#endif
}

/** restore basis */
void LPSoplexInterface::RestorePreStrongbranchingBasis() {
  assert(row_basis_status_.size() == spx_->numRowsReal());
  assert(col_basis_status_.size() == spx_->numColsReal());

  try {
    spx_->setBasis(row_basis_status_.get_ptr(), col_basis_status_.get_ptr());
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    /* MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception: %s\n", s.c_str()); */
#else
  catch (const SPxException&) {
#endif
    /* since it is not clear if the status in SoPlex are set correctly we want to make sure that if an error is thrown
     * the status is not OPTIMAL anymore. */
    assert(spx_->status() != SPxSolver::OPTIMAL);
  }
}

void LPSoplexInterface::InvalidateSolution() {
  solved_ = false;
}

bool LPSoplexInterface::PreStrongBranchingBasisFreed() const {
  return ((row_basis_status_.size() == 0) && (col_basis_status_.size() == 0));
}

void LPSoplexInterface::FreePreStrongBranchingBasis() {
  row_basis_status_.clear();
  col_basis_status_.clear();
}

RetCode LPSoplexInterface::StrongBranch(
  LPIndex col,                       /**< column to apply strong branching on */
  LPValue primal_sol,              /**< current primal solution value of column */
  LPNum iteration_limit,           /**< iteration limit for strong branchings */
  LPValue& dual_bound_down_branch, /**< stores dual bound after branching column down */
  LPValue& dual_bound_up_branch,   /**< stores dual bound after branching column up */
  bool& down_valid,                /**< stores whether the returned down value is a valid dual bound;
                                    *   otherwise, it can only be used as an estimate value */
  bool& up_valid,                  /**< stores whether the returned up value is a valid dual bound;
                                    *   otherwise, it can only be used as an estimate value */
  LPNum& iterations                /**< stores total number of strong branching iterations, or -1 */
) {
  SPxSolver::Status status;
  LPValue old_lb;
  LPValue old_ub;
  LPValue new_lb;
  LPValue new_ub;
  bool from_parent_basis;
  bool error;
  int old_iter_limit;
  SPxOut::Verbosity verbosity;

  /* store and set verbosity */
  verbosity = spx_->spxout.getVerbosity();
  spx_->spxout.setVerbosity((SPxOut::Verbosity)(GetLPInfo() ? SOPLEX_VERBLEVEL : 0));

  MiniMIPdebugMessage("calling StrongBranch() on variable %d (%d iterations)\n", col, iteration_limit);

#ifndef STRONGBRANCH_RESTOREBASIS
  from_parent_basis = false;
#endif
  error = false;
  old_iter_limit = spx_->intParam(SoPlex::ITERLIMIT);

  /* get current bounds of column */
  old_lb = spx_->lowerReal(col);
  old_ub = spx_->upperReal(col);
  down_valid = false;
  up_valid = false;
  iterations = 0;

  /* set the algorithm type to use dual simplex */
  static_cast<void>(spx_->setIntParam(SoPlex::ALGORITHM, SoPlex::ALGORITHM_DUAL));

  /* down branch */
  new_ub = EPSCEIL(primal_sol - 1.0, FeasibilityTolerance());

  if (new_ub >= old_lb - 0.5) {
    MiniMIPdebugMessage("strong branching down on x%d (%g) with %d iterations\n", col, primal_sol, iteration_limit);

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
          down_valid = true;
          MiniMIPdebugMessage(" --> Terminate with value %f\n", dual_bound_down_branch);
          break;
        case SPxSolver::ABORT_TIME: /* SoPlex does not return a proven dual bound, if it is aborted */
        case SPxSolver::ABORT_ITER:
        case SPxSolver::ABORT_CYCLING:
        case SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS:
          dual_bound_down_branch = spx_->objValueReal();
          break;
        case SPxSolver::ABORT_VALUE:
        case SPxSolver::INFEASIBLE:
          dual_bound_down_branch = GetObjectiveLimit();
          down_valid = true;
          break;
        default:
          error = true;
          break;
      } /*lint !e788*/
      iterations += spx_->numIterations();

#ifdef STRONGBRANCH_RESTOREBASIS
      /* we restore the pre-strong-branching basis by default (and don't solve again) */
      assert(!PreStrongBranchingBasisFreed());
      RestorePreStrongbranchingBasis();
      from_parent_basis = false;
#else
      /* if cycling or singular basis occured and we started not from the pre-strong-branching basis, then we restore the
       * pre-strong-branching basis and try again with reduced iteration limit */
      repeat_strong_branching = ((status == SPxSolver::ABORT_CYCLING || status == SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS || status == SPxSolver::SINGULAR) && !from_parent_basis && spx_->numIterations() < static_cast<int>(iteration_limit));

      if (repeat_strong_branching) {
        MiniMIPdebugMessage(" --> Repeat strong branching down with %d iterations after restoring basis\n",
           iteration_limit - spx_->numIterations());
        spx_->setIntParam(SoPlex::ITERLIMIT, iteration_limit - spx_->numIterations());
        RestorePreStrongbranchingBasis();
        from_parent_basis = true;
        error = false;
      }
      /* otherwise don't solve again */
      else
        from_parent_basis = false;
#endif
    } while (from_parent_basis);

    spx_->changeUpperReal(col, old_ub);
    assert(spx_->lowerReal(col) <= spx_->upperReal(col));
  } else {
    dual_bound_down_branch = GetObjectiveLimit();
    down_valid = true;
  }

  /* up branch */
  if (!error) {
    new_lb = EPSFLOOR(primal_sol + 1.0, FeasibilityTolerance());
    if (new_lb <= old_ub + 0.5) {
      MiniMIPdebugMessage("strong branching  up  on x%d (%g) with %d iterations\n", col, primal_sol, iteration_limit);

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
            up_valid = true;
            MiniMIPdebugMessage(" --> Terminate with value %f\n", spx_->objValueReal());
            break;
          case SPxSolver::ABORT_TIME: /* SoPlex does not return a proven dual bound, if it is aborted */
          case SPxSolver::ABORT_ITER:
          case SPxSolver::ABORT_CYCLING:
          case SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS:
            dual_bound_up_branch = spx_->objValueReal();
            break;
          case SPxSolver::ABORT_VALUE:
          case SPxSolver::INFEASIBLE:
            dual_bound_up_branch = GetObjectiveLimit();
            up_valid = true;
            break;
          default:
            error = true;
            break;
        } /*lint !e788*/
        iterations += spx_->numIterations();

#ifdef STRONGBRANCH_RESTOREBASIS
        /* we restore the pre-strong-branching basis by default (and don't solve again) */
        assert(!PreStrongBranchingBasisFreed());
        RestorePreStrongbranchingBasis();
        from_parent_basis = false;
#else
        /* if cycling or singular basis occured and we started not from the pre-strong-branching basis, then we restore the
         * pre-strong-branching basis and try again with reduced iteration limit */
        repeat_strong_branching = ((status == SPxSolver::ABORT_CYCLING || status == SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS || status == SPxSolver::SINGULAR) && !from_parent_basis && spx_->numIterations() < static_cast<int>(iteration_limit));

        if (repeat_strong_branching) {
          MiniMIPdebugMessage(" --> Repeat strong branching  up  with %d iterations after restoring basis\n", iteration_limit - spx_->numIterations());
          RestorePreStrongbranchingBasis();
          spx_->setIntParam(SoPlex::ITERLIMIT, iteration_limit - spx_->numIterations());
          error = false;
          from_parent_basis = true;
        }
        /* otherwise don't solve again */
        else
          from_parent_basis = false;
#endif
      } while (from_parent_basis);

      spx_->changeLowerReal(col, old_lb);
      assert(spx_->lowerReal(col) <= spx_->upperReal(col));
    } else {
      dual_bound_up_branch = GetObjectiveLimit();
      up_valid = true;
    }
  }
  /* reset old iteration limit */
  static_cast<void>(spx_->setIntParam(SoPlex::ITERLIMIT, old_iter_limit));

  /* restore verbosity */
  spx_->spxout.setVerbosity(verbosity);

  if (error) {
    MiniMIPdebugMessage("StrongBranch() returned SoPlex status %d\n", int(status)); /*lint !e644*/
    return RetCode::kLPError;
  }
  return RetCode::kOkay;
}

RetCode LPSoplexInterface::LoadColumnLP(
  LPObjectiveSense obj_sense,           /**< objective sense */
  LPNum num_cols,                       /**< number of columns */
  const LPValueArray& objective_values, /**< objective function values of columns */
  const LPValueArray& lower_bounds,     /**< lower bounds of columns */
  const LPValueArray& upper_bounds,     /**< upper bounds of columns */
  StringArray& col_names,               /**< column names */
  LPNum num_rows,                       /**< number of rows */
  const LPValueArray& left_hand_sides,  /**< left hand sides of rows */
  const LPValueArray& right_hand_sides, /**< right hand sides of rows */
  StringArray& row_names,               /**< row names */
  LPNum num_non_zeros,                  /**< number of non-zero elements in the constraint matrix */
  const LPIndexArray& begin_cols,       /**< start index of each column in row_indices- and vals-array */
  const LPIndexArray& row_indices,      /**< row indices of constraint matrix entries */
  const LPValueArray& vals              /**< values of constraint matrix entries */
) {
#ifndef NDEBUG
   for (int j = 0; j < static_cast<int>(num_non_zeros); j++)
      assert(vals[j] != 0);
#endif

  MiniMIPdebugMessage("calling LoadColumnLP()\n");

  assert(PreStrongBranchingBasisFreed());

  try {
    LPRowSet rows(num_rows);
    DSVector empty_vector(0);
    int i;

    spx_->clearLPReal();

    /* set objective sense */
    static_cast<void>(spx_->setIntParam(SoPlex::OBJSENSE, (obj_sense == LPObjectiveSense::kMinimize ? SoPlex::OBJSENSE_MINIMIZE : SoPlex::OBJSENSE_MAXIMIZE)));

    /* create empty rows with given sides */
    for (i = 0; i < static_cast<int>(num_rows); ++i)
      rows.add(left_hand_sides[i], empty_vector, right_hand_sides[i]);
    spx_->addRowsReal(rows);

    /* create column vectors with coefficients and bounds */
    MINIMIP_CALL(AddColumns(num_cols, objective_values, lower_bounds, upper_bounds, col_names, num_non_zeros, begin_cols, row_indices, vals));
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    /* MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception: %s\n", s.c_str()); */
#else
  catch (const SPxException&) {
#endif
    return RetCode::kLPError;
  }
  return RetCode::kOkay;
}

/** adds columns to the LP
 *
 *  @note The indices array is not checked for duplicates, problems may appear if indices are added more than once.
 */
RetCode LPSoplexInterface::AddColumns(
  LPNum num_cols,                       /**< number of columns to be added */
  const LPValueArray& objective_values, /**< objective function values of new columns */
  const LPValueArray& lower_bounds,     /**< lower bounds of new columns */
  const LPValueArray& upper_bounds,     /**< upper bounds of new columns */
  StringArray& col_names,               /**< column names */
  LPNum num_non_zeros,                  /**< number of non-zero elements to be added to the constraint matrix */
  const LPIndexArray& begin_cols,       /**< start index of each column in indices- and vals-array */
  const LPIndexArray& indices,          /**< row indices of constraint matrix entries */
  const LPValueArray& vals              /**< values of constraint matrix entries */
) {
  MiniMIPdebugMessage("calling AddColumns()\n");

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

#ifndef NDEBUG
  if (num_non_zeros > 0) {
    /* perform check that no new rows are added - this is likely to be a mistake */
    int num_rows = spx_->numRowsReal();
    for (unsigned int j = 0; j < num_non_zeros; ++j) {
      assert(0 <= indices[j] && indices[j] < static_cast<unsigned int>(num_rows));
      assert(vals[j] != 0.0);
    }
  }
#endif

  try {
    IntArray integer_indices(indices.begin(), indices.end());
    LPColSet cols(num_cols);
    DSVector col_Vector(num_cols);
    int start;
    int last;
    int i;

    /* create column vectors with coefficients and bounds */
    for (i = 0; i < static_cast<int>(num_cols); ++i) {
      col_Vector.clear();
      if (num_non_zeros > 0) {
        start = begin_cols[i];
        last = (i == static_cast<int>(num_cols) - 1 ? num_non_zeros : begin_cols[i + 1]);
        col_Vector.add(last - start, &integer_indices[start], &vals[start]);
      }
      cols.add(objective_values[i], lower_bounds[i], col_Vector, upper_bounds[i]);
    }
    spx_->addColsReal(cols);
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    /* MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception: %s\n", s.c_str()); */
#else
  catch (const SPxException&) {
#endif
    return RetCode::kLPError;
  }

  return RetCode::kOkay;
}

/** deletes all columns in the given range from LP */
RetCode LPSoplexInterface::DeleteColumns(
  LPIndex first_col, /**< first column to be deleted */
  LPIndex last_col   /**< last column to be deleted */
) {
  MiniMIPdebugMessage("calling DeleteColumns()\n");

  assert(0 <= first_col && first_col <= last_col && static_cast<int>(last_col) < spx_->numColsReal());

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

  SOPLEX_TRY(spx_->removeColRangeReal(first_col, last_col));

  return RetCode::kOkay;
}

/** deletes columns from LP; the new position of a column must not be greater than its old position */
RetCode LPSoplexInterface::DeleteColumnSet(
  BoolArray& deletion_status /**< deletion status of columns */
) {
  int num_cols;
  int i;

  MiniMIPdebugMessage("calling DeleteColumnSet()\n");

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

  IntArray int_deletion_status(deletion_status.begin(), deletion_status.end());
  num_cols = spx_->numColsReal();

  /* SoPlex removeCols() method deletes the columns with deletion_status[i] < 0, so we have to negate the values */
  for (i = 0; i < num_cols; ++i)
    int_deletion_status[i] *= -1;

  SOPLEX_TRY(spx_->removeColsReal(int_deletion_status.data()));

  return RetCode::kOkay;
}

/** adds rows to the LP
 *
 *  @note The indices array is not checked for duplicates, problems may appear if indices are added more than once.
 */
RetCode LPSoplexInterface::AddRows(
  LPNum num_rows,                       /**< number of rows to be added */
  const LPValueArray& left_hand_sides,  /**< left hand sides of new rows */
  const LPValueArray& right_hand_sides, /**< right hand sides of new rows */
  StringArray& row_names,               /**< row names */
  LPNum num_non_zeros,                  /**< number of non-zero elements to be added to the constraint matrix */
  const LPIndexArray& begin_rows,       /**< start index of each row in indices- and vals-array */
  const LPIndexArray& indices,          /**< column indices of constraint matrix entries */
  const LPValueArray& vals              /**< values of constraint matrix entries */
) {
  MiniMIPdebugMessage("calling AddRows()\n");

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

#ifndef NDEBUG
  if (num_non_zeros > 0) {
    /* perform check that no new columns are added - this is likely to be a mistake */
    int num_cols = spx_->numColsReal();
    for (int j = 0; j < static_cast<int>(num_non_zeros); ++j) {
      assert(vals[j] != 0.0);
      assert(0 <= indices[j] && static_cast<int>(indices[j]) < num_cols);
    }
  }
#endif

  try {
    LPRowSet rows(num_rows);
    DSVector row_Vector;
    unsigned int start;
    unsigned int last;
    unsigned int i;
    IntArray integer_indices(indices.begin(), indices.end());

    /* create row vectors with given sides */
    for (i = 0; i < num_rows; ++i) {
      row_Vector.clear();
      if (num_non_zeros > 0) {
        start = begin_rows[i];
        last = (i == num_rows - 1 ? num_non_zeros : begin_rows[i + 1]);
        row_Vector.add(last - start, &integer_indices[start], &vals[start]);
      }
      rows.add(left_hand_sides[i], row_Vector, right_hand_sides[i]);
    }
    spx_->addRowsReal(rows);
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    /* MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception: %s\n", s.c_str()); */
#else
  catch (const SPxException&) {
#endif
    return RetCode::kLPError;
  }
  return RetCode::kOkay;
}

/** deletes all rows in the given range from LP */
RetCode LPSoplexInterface::DeleteRows(
  LPIndex first_row, /**< first row to be deleted */
  LPIndex last_row   /**< last row to be deleted */
) {
  MiniMIPdebugMessage("calling DeleteRows()\n");

  assert(0 <= first_row && first_row <= last_row && static_cast<int>(last_row) < spx_->numRowsReal());

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

  SOPLEX_TRY(spx_->removeRowRangeReal(first_row, last_row));

  return RetCode::kOkay;
}

/** deletes rows from LP; the new position of a row must not be greater that its old position */
RetCode LPSoplexInterface::DeleteRowSet(
  BoolArray& deletion_status /**< deletion status of rows */
) {
  int num_rows;
  int i;

  MiniMIPdebugMessage("calling DeleteRowSet()\n");

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

  IntArray int_deletion_status(deletion_status.begin(), deletion_status.end());

  num_rows = spx_->numRowsReal();

  /* SoPlex removeRows() method deletes the rows with deletion_status[i] < 0, so we have to negate the values */
  for (i = 0; i < num_rows; ++i)
    int_deletion_status[i] *= -1;

  SOPLEX_TRY(spx_->removeRowsReal(int_deletion_status.data()));

  return RetCode::kOkay;
}

/** clears the whole LP */
RetCode LPSoplexInterface::Clear() {
  MiniMIPdebugMessage("calling Clear()\n");

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());
  SOPLEX_TRY(spx_->clearLPReal());

  return RetCode::kOkay;
}

/** clears current LP Interface state (like basis information) of the solver */
RetCode LPSoplexInterface::ClearState() {
  /*lint --e{715}*/
  MiniMIPdebugMessage("calling ClearState()\n");

  try {
    spx_->clearBasis();
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    /*  MiniMIPmessagePrintWarning(lpi->messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());*/
#else
  catch (const SPxException&) {
#endif
    assert(spx_->status() != SPxSolver::OPTIMAL);
    return RetCode::kLPError;
  }
  return RetCode::kOkay;
}

/** changes lower and upper bounds of columns */
RetCode LPSoplexInterface::ChangeBounds(
  LPNum num_cols,                   /**< number of columns to change bounds for */
  const LPIndexArray& indices,      /**< column indices */
  const LPValueArray& lower_bounds, /**< values for the new lower bounds */
  const LPValueArray& upper_bounds  /**< values for the new upper bounds */
) {
  unsigned int i;

  MiniMIPdebugMessage("calling ChangeBounds()\n");

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

  try {
    for (i = 0; i < num_cols; ++i) {
      assert(0 <= indices[i] && static_cast<int>(indices[i]) < spx_->numColsReal());

      if (IsInfinity(lower_bounds[i])) {
        MiniMIPerrorMessage("LP Error: fixing lower bound for variable %d to infinity.\n", indices[i]);
        return RetCode::kLPError;
      }
      if (IsInfinity(-upper_bounds[i])) {
        MiniMIPerrorMessage("LP Error: fixing upper bound for variable %d to -infinity.\n", indices[i]);
        return RetCode::kLPError;
      }
      spx_->changeBoundsReal(indices[i], lower_bounds[i], upper_bounds[i]);
      assert(spx_->lowerReal(indices[i]) <= spx_->upperReal(indices[i]) + spx_->realParam(SoPlex::EPSILON_ZERO));
    }
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    /* MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception: %s\n", s.c_str());*/
#else
  catch (const SPxException&) {
#endif
    return RetCode::kLPError;
  }
  return RetCode::kOkay;
}

/** changes left and right hand sides of rows */
RetCode LPSoplexInterface::ChangeSides(
  LPNum num_rows,                      /**< number of rows to change sides for */
  const LPIndexArray& indices,         /**< row indices */
  const LPValueArray& left_hand_sides, /**< new values for left hand sides */
  const LPValueArray& right_hand_sides /**< new values for right hand sides */
) {
  unsigned int i;

  MiniMIPdebugMessage("calling ChangeSides()\n");

  if (num_rows <= 0)
    return RetCode::kOkay;

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

  try {
    for (i = 0; i < num_rows; ++i) {
      assert(0 <= indices[i] && static_cast<int>(indices[i]) < spx_->numRowsReal());
      spx_->changeRangeReal(indices[i], left_hand_sides[i], right_hand_sides[i]);
      assert(spx_->lhsReal(indices[i]) <= spx_->rhsReal(indices[i]) + spx_->realParam(SoPlex::EPSILON_ZERO));
    }
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    /* MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception: %s\n", s.c_str()); */
#else
  catch (const SPxException&) {
#endif
    return RetCode::kLPError;
  }
  return RetCode::kOkay;
}

/** changes the objective sense */
RetCode LPSoplexInterface::ChangeObjectiveSense(
  LPObjectiveSense obj_sense /**< new objective sense */
) {
  MiniMIPdebugMessage("calling ChangeObjectiveSense()\n");

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

  SOPLEX_TRY(static_cast<void>(spx_->setIntParam(SoPlex::OBJSENSE, obj_sense == LPObjectiveSense::kMinimize ? SoPlex::OBJSENSE_MINIMIZE : SoPlex::OBJSENSE_MAXIMIZE)));

  return RetCode::kOkay;
}

/** changes objective values of columns in the LP */
RetCode LPSoplexInterface::ChangeObjective(
  LPNum num_cols,                  /**< number of columns to change objective value for */
  const LPIndexArray& indices,     /**< column indices to change objective value for */
  const LPValueArray& new_obj_vals /**< new objective values for columns */
) {
  unsigned int i;

  MiniMIPdebugMessage("calling ChangeObjective()\n");

  InvalidateSolution();

  assert(PreStrongBranchingBasisFreed());

  try {
    for (i = 0; i < num_cols; ++i) {
      assert(0 <= indices[i] && static_cast<int>(indices[i]) < spx_->numColsReal());
      spx_->changeObjReal(indices[i], new_obj_vals[i]);
    }
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    /* MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception: %s\n", s.c_str()); */
#else
  catch (const SPxException&) {
#endif
    return RetCode::kLPError;
  }
  return RetCode::kOkay;
}

/**@} */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of rows in the LP */
LPNum LPSoplexInterface::GetNumberOfRows() {
  MiniMIPdebugMessage("calling GetNumberOfRows()\n");

  return spx_->numRowsReal();
}

/** gets the number of columns in the LP */
LPNum LPSoplexInterface::GetNumberOfColumns() {
  MiniMIPdebugMessage("calling GetNumberOfColumns()\n");

  return spx_->numColsReal();
}

/** gets the number of nonzero elements in the LP constraint matrix */
LPNum LPSoplexInterface::GetNumberOfNonZeros() {
  int i;
  /* SoPlex has no direct method to return the number of nonzeros, so we have to count them manually */
  LPNum num_non_zeros = 0;

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

/** gets the objective sense of the LP */
LPObjectiveSense LPSoplexInterface::GetObjectiveSense() {
  MiniMIPdebugMessage("calling GetObjectiveSense()\n");

  return (spx_->intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE) ? LPObjectiveSense::kMinimize : LPObjectiveSense::kMaximize;
}

/** gets columns from LP problem object
 *
 *  Either both, lower_bound and upper_bound, have to be 0, or both have to be non-0,
 *  either n_non_zeroes, begin_cols, indices, and obj_coeffs have to be 0, or all of them have to be non-0.
 */
RetCode LPSoplexInterface::GetColumns(
  LPIndex first_col,            /**< first column to get from LP */
  LPIndex last_col,             /**< last column to get from LP */
  LPValueArray& lower_bounds, /**< array to store the lower bound vector */
  LPValueArray& upper_bounds, /**< array to store the upper bound vector */
  LPNum& num_non_zeros,       /**< store the number of non-zero elements */
  LPIndexArray& begin_cols,   /**< array to store start index of each column in indices- and vals-array */
  LPIndexArray& indices,      /**< array to store row indices of constraint matrix entries */
  LPValueArray& vals          /**< array to store values of constraint matrix entries */
) {
  unsigned int i;
  int j;

  MiniMIPdebugMessage("calling GetColumns()\n");

  assert(0 <= first_col && first_col <= last_col && static_cast<int>(last_col) < spx_->numColsReal());

  if (spx_->boolParam(SoPlex::PERSISTENTSCALING)) {
    DVector tmp_lower_bound(spx_->numColsReal());
    DVector tmp_upper_bound(spx_->numColsReal());
    spx_->getLowerReal(tmp_lower_bound);
    spx_->getUpperReal(tmp_upper_bound);
    for (i = first_col; i <= last_col; ++i) {
      lower_bounds[i - first_col] = tmp_lower_bound[i];
      upper_bounds[i - first_col] = tmp_upper_bound[i];
    }
  } else {
    const Vector& tmp_lower_bound = spx_->lowerRealInternal();
    const Vector& tmp_upper_bound = spx_->upperRealInternal();
    for (i = first_col; i <= last_col; ++i) {
      lower_bounds[i - first_col] = tmp_lower_bound[i];
      upper_bounds[i - first_col] = tmp_upper_bound[i];
    }
  }

  num_non_zeros = 0;
  for (i = first_col; i <= last_col; ++i) {
    begin_cols[i - first_col] = num_non_zeros;

    if (spx_->boolParam(SoPlex::PERSISTENTSCALING)) {
      DSVector cvec;
      spx_->getColVectorReal(i, cvec);
      for (j = 0; j < cvec.size(); ++j) {
        indices[num_non_zeros] = cvec.index(j);
        vals[num_non_zeros] = cvec.value(j);
        num_non_zeros++;
      }
    } else {
      const SVector& cvec = spx_->colVectorRealInternal(i);
      for (j = 0; j < cvec.size(); ++j) {
        indices[num_non_zeros] = cvec.index(j);
        vals[num_non_zeros] = cvec.value(j);
        num_non_zeros++;
      }
    }
  }
  return RetCode::kOkay;
}

/** gets rows from LP problem object
 *
 *  Either both, left_hand_side and right_hand_side, have to be 0, or both have to be non-0,
 *  either n_non_zeroes, begin_cols, indices, and obj_coeffs have to be 0, or all of them have to be non-0.
 */
RetCode LPSoplexInterface::GetRows(
  LPIndex first_row,                /**< first row to get from LP */
  LPIndex last_row,                 /**< last row to get from LP */
  LPValueArray& left_hand_sides,  /**< array to store left hand side vector */
  LPValueArray& right_hand_sides, /**< array to store right hand side vector */
  LPNum& num_non_zeros,           /**< store the number of non-zero elements */
  LPIndexArray& begin_cols,       /**< array to store start index of each row in indices- and vals-array */
  LPIndexArray& indices,          /**< array to store column indices of constraint matrix entries */
  LPValueArray& vals              /**< array to store values of constraint matrix entries */
) {
  unsigned int i;
  int j;

  MiniMIPdebugMessage("calling GetRows()\n");

  assert(0 <= first_row && first_row <= last_row && static_cast<int>(last_row) < spx_->numRowsReal());

  if (spx_->boolParam(SoPlex::PERSISTENTSCALING)) {
    DVector lhsvec(spx_->numRowsReal());
    DVector rhsvec(spx_->numRowsReal());
    spx_->getLhsReal(lhsvec);
    spx_->getRhsReal(rhsvec);
    for (i = first_row; i <= last_row; ++i) {
      left_hand_sides[i - first_row] = lhsvec[i];
      right_hand_sides[i - first_row] = rhsvec[i];
    }
  } else {
    const Vector& lhsvec = spx_->lhsRealInternal();
    const Vector& rhsvec = spx_->rhsRealInternal();
    for (i = first_row; i <= last_row; ++i) {
      left_hand_sides[i - first_row] = lhsvec[i];
      right_hand_sides[i - first_row] = rhsvec[i];
    }
  }

  num_non_zeros = 0;
  for (i = first_row; i <= last_row; ++i) {
    begin_cols[i - first_row] = num_non_zeros;

    if (spx_->boolParam(SoPlex::PERSISTENTSCALING)) {
      DSVector rvec;
      spx_->getRowVectorReal(i, rvec);
      for (j = 0; j < rvec.size(); ++j) {
        indices[num_non_zeros] = rvec.index(j);
        vals[num_non_zeros] = rvec.value(j);
        num_non_zeros++;
      }
    } else {
      const SVector& rvec = spx_->rowVectorRealInternal(i);
      for (j = 0; j < rvec.size(); ++j) {
        indices[num_non_zeros] = rvec.index(j);
        vals[num_non_zeros] = rvec.value(j);
        num_non_zeros++;
      }
    }
  }
  return RetCode::kOkay;
}

/** gets objective coefficients from LP problem object */
RetCode LPSoplexInterface::GetObjective(
  LPIndex first_col,         /**< first column to get objective coefficient for */
  LPIndex last_col,          /**< last column to get objective coefficient for */
  LPValueArray& obj_coeffs /**< array to store objective coefficients */
) {
  unsigned int i;

  MiniMIPdebugMessage("calling GetObjective()\n");

  assert(0 <= first_col && first_col <= last_col && static_cast<int>(last_col) < spx_->numColsReal());

  for (i = first_col; i <= last_col; ++i)
    obj_coeffs[i - first_col] = spx_->objReal(i);

  return RetCode::kOkay;
}

/** gets current bounds from LP problem object */
RetCode LPSoplexInterface::GetBounds(
  LPIndex first_col,            /**< first column to get bounds for */
  LPIndex last_col,             /**< last column to get bounds for */
  LPValueArray& lower_bounds, /**< array to store lower bound values */
  LPValueArray& upper_bounds  /**< array to store upper bound values */
) {
  unsigned int i;

  MiniMIPdebugMessage("calling GetBounds()\n");

  assert(0 <= first_col && first_col <= last_col && static_cast<int>(last_col) < spx_->numColsReal());

  for (i = first_col; i <= last_col; ++i) {
    lower_bounds[i - first_col] = spx_->lowerReal(i);
    upper_bounds[i - first_col] = spx_->upperReal(i);
  }

  return RetCode::kOkay;
}

/** gets current row sides from LP problem object */
RetCode LPSoplexInterface::GetSides(
  LPIndex first_row,               /**< first row to get sides for */
  LPIndex last_row,                /**< last row to get sides for */
  LPValueArray& left_hand_sides, /**< array to store left hand side values */
  LPValueArray& right_hand_sides /**< array to store right hand side values */
) {
  unsigned int i;

  MiniMIPdebugMessage("calling GetSides()\n");

  assert(0 <= first_row && first_row <= last_row && static_cast<int>(last_row) < spx_->numRowsReal());

  for (i = first_row; i <= last_row; ++i) {
    left_hand_sides[i - first_row] = spx_->lhsReal(i);
    right_hand_sides[i - first_row] = spx_->rhsReal(i);
  }
  return RetCode::kOkay;
}

/** gets a single coefficient */
RetCode LPSoplexInterface::GetCoefficient(
  LPIndex row,   /**< row number of coefficient */
  LPIndex col,   /**< column number of coefficient */
  LPValue& val /**< array to store the value of the coefficient */
) {
  MiniMIPdebugMessage("calling GetCoefficients()\n");

  assert(0 <= col && static_cast<int>(col) < static_cast<int>(spx_->numColsReal()));
  assert(0 <= row && static_cast<int>(row) < static_cast<int>(spx_->numRowsReal()));

  val = spx_->coefReal(row, col);

  return RetCode::kOkay;
}

/**@} */

/**@name Solving Methods */
/**@{ */

/** calls primal simplex to solve the LP */
RetCode LPSoplexInterface::SolvePrimal() {
  MiniMIPdebugMessage("calling SolvePrimal()\n");

  static_cast<void>(spx_->setIntParam(SoPlex::ALGORITHM, SoPlex::ALGORITHM_PRIMAL));
  return SoPlexSolve();
}

/** calls dual simplex to solve the LP */
RetCode LPSoplexInterface::SolveDual() {
  MiniMIPdebugMessage("calling SolveDual()\n");

  static_cast<void>(spx_->setIntParam(SoPlex::ALGORITHM, SoPlex::ALGORITHM_DUAL));
  return SoPlexSolve();
}

/** start strong branching - call before any strong branching */
RetCode LPSoplexInterface::StartStrongbranch() {

  assert(PreStrongBranchingBasisFreed());
  SavePreStrongbranchingBasis();

  return RetCode::kOkay;
}

/** end strong branching - call after any strong branching */
RetCode LPSoplexInterface::EndStrongbranch() {

  assert(!PreStrongBranchingBasisFreed());
  RestorePreStrongbranchingBasis();
  FreePreStrongBranchingBasis();

  return RetCode::kOkay;
}

/** performs strong branching iterations on one @b fractional candidate */
RetCode LPSoplexInterface::StrongbranchFractionalValue(
  LPIndex col,                       /**< column to apply strong branching on */
  LPValue primal_sol,              /**< fractional current primal solution value of column */
  LPNum iteration_limit,           /**< iteration limit for strong branchings */
  LPValue& dual_bound_down_branch, /**< stores dual bound after branching column down */
  LPValue& dual_bound_up_branch,   /**< stores dual bound after branching column up */
  bool& down_valid,                /**< whether the returned down value is a valid dual bound; otherwise, it can only be used as an estimate value */
  bool& up_valid,                  /**< whether the returned up value is a valid dual bound; otherwise, it can only be used as an estimate value */
  LPNum& iterations                /**< stores total number of strong branching iterations */
) {

  RetCode retcode;

  /* pass call on to StrongBranch() */
  retcode = StrongBranch(col, primal_sol, iteration_limit, dual_bound_down_branch, dual_bound_up_branch, down_valid, up_valid, iterations);

  /* pass RetCode::kLPError to MiniMIP without a back trace */
  if (retcode == RetCode::kLPError)
    return RetCode::kLPError;

  /* evaluate retcode */
  MINIMIP_CALL(retcode);

  return RetCode::kOkay;
}

/** performs strong branching iterations on given @b fractional candidates */
RetCode LPSoplexInterface::StrongbranchFractionalValues(
  LPIndexArray& cols,                       /**< columns to apply strong branching on */
  LPNum num_cols,                         /**< number of columns */
  LPValueArray& primal_sols,              /**< fractional current primal solution values of columns */
  LPNum iteration_limit,                  /**< iteration limit for strong branchings */
  LPValueArray& dual_bound_down_branches, /**< stores dual bounds after branching columns down */
  LPValueArray& dual_bound_up_branches,   /**< stores dual bounds after branching columns up */
  BoolArray& down_valids,                 /**< stores whether the returned down values are valid dual bounds;
                                           *   otherwise, they can only be used as an estimate values */
  BoolArray& up_valids,                   /**< stores whether the returned up values are a valid dual bounds;
                                           *   otherwise, they can only be used as an estimate values */
  LPNum& iterations                       /**< stores total number of strong branching iterations */
) {
  RetCode retcode;
  iterations = 0;
  bool up, down;

  for (int j = 0; j < static_cast<int>(num_cols); ++j) {
    /* pass call on to StrongBranch() */
    retcode = StrongBranch(cols[j], primal_sols[j], iteration_limit, (dual_bound_down_branches[j]), (dual_bound_up_branches[j]), down, up, iterations);
    down_valids[j] = down;
    up_valids[j] = up;
    /* pass RetCode::kLPError to MiniMIP without a back trace */
    if (retcode == RetCode::kLPError)
      return RetCode::kLPError;

    /* evaluate retcode */
    MINIMIP_CALL(retcode);
  }
  return RetCode::kOkay;
}

/** performs strong branching iterations on one candidate with @b integral value */
RetCode LPSoplexInterface::StrongbranchIntegerValue(
  LPIndex col,                       /**< column to apply strong branching on */
  LPValue primal_sol,              /**< current integral primal solution value of column */
  LPNum iteration_limit,           /**< iteration limit for strong branchings */
  LPValue& dual_bound_down_branch, /**< stores dual bound after branching column down */
  LPValue& dual_bound_up_branch,   /**< stores dual bound after branching column up */
  bool& down_valid,                /**< stores whether the returned down value is a valid dual bound;
                                    *   otherwise, it can only be used as an estimate value */
  bool& up_valid,                  /**< stores whether the returned up value is a valid dual bound;
                                    *   otherwise, it can only be used as an estimate value */
  LPNum& iterations                /**< stores total number of strong branching iterations */
) {
  RetCode retcode;

  /* pass call on to StrongBranch() */
  retcode = StrongBranch(col, primal_sol, iteration_limit, dual_bound_down_branch, dual_bound_up_branch, down_valid, up_valid, iterations);

  /* pass RetCode::kLPError to MiniMIP without a back trace */
  if (retcode == RetCode::kLPError)
    return RetCode::kLPError;

  /* evaluate retcode */
  MINIMIP_CALL(retcode);

  return RetCode::kOkay;
}

/** performs strong branching iterations on given candidates with @b integral values */
RetCode LPSoplexInterface::StrongbranchIntegerValues(
  LPIndexArray& cols,                       /**< columns to apply strong branching on */
  LPNum num_cols,                         /**< number of columns */
  LPValueArray& primal_sols,              /**< current integral primal solution values of columns */
  LPNum iteration_limit,                  /**< iteration limit for strong branchings */
  LPValueArray& dual_bound_down_branches, /**< stores dual bounds after branching columns down */
  LPValueArray& dual_bound_up_branches,   /**< stores dual bounds after branching columns up */
  BoolArray& down_valids,                 /**< stores whether the returned down values are valid dual bounds;
                                           *   otherwise, they can only be used as an estimate values */
  BoolArray& up_valids,                   /**< stores whether the returned up values are a valid dual bounds;
                                           *   otherwise, they can only be used as an estimate values */
  LPNum& iterations                       /**< stores total number of strong branching iterations */
) {
  RetCode retcode;
  iterations = 0;
  bool up, down;

  for (int j = 0; j < static_cast<int>(num_cols); ++j) {
    /* pass call on to StrongBranch() */
    retcode = StrongBranch(cols[j], primal_sols[j], iteration_limit, (dual_bound_down_branches[j]), (dual_bound_up_branches[j]), down, up, iterations);
    down_valids[j] = down;
    up_valids[j] = up;

    /* pass RetCode::kLPError to MiniMIP without a back trace */
    if (retcode == RetCode::kLPError)
      return RetCode::kLPError;

    /* evaluate retcode */
    MINIMIP_CALL(retcode);
  }

  return RetCode::kOkay;
}
/**@} */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the LP */
bool LPSoplexInterface::WasSolved() {

  return solved_;
}

/** returns true iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily meanthat the solver knows and can return the primal ray
 */
bool LPSoplexInterface::ExistsPrimalRay() {
  MiniMIPdebugMessage("calling ExistsPrimalRay()\n");

  return (spx_->status() == SPxSolver::UNBOUNDED);
}

/** returns true iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
bool LPSoplexInterface::HasPrimalRay() {
  MiniMIPdebugMessage("calling HasPrimalRay()\n");

  return spx_->hasPrimalRay();
}

/** returns true iff LP is proven to be primal unbounded */
bool LPSoplexInterface::IsPrimalUnbounded() {
  MiniMIPdebugMessage("calling IsPrimalUnbounded()\n");

  assert(spx_->status() != SPxSolver::UNBOUNDED || spx_->basisStatus() == SPxBasis::UNBOUNDED);

  /* if SoPlex returns unbounded, this may only mean that an unbounded ray is available, not necessarily a primal
    * feasible point; hence we have to check the perturbation
    */
  return spx_->status() == SPxSolver::UNBOUNDED;
}

/** returns true iff LP is proven to be primal infeasible */
bool LPSoplexInterface::IsPrimalInfeasible() {
  MiniMIPdebugMessage("calling IsPrimalInfeasible()\n");

  return (spx_->status() == SPxSolver::INFEASIBLE);
}

/** returns true iff LP is proven to be primal feasible */
bool LPSoplexInterface::IsPrimalFeasible() {
  MiniMIPdebugMessage("calling IsPrimalFeasible()\n");

  return spx_->basisStatus() == SPxBasis::OPTIMAL || spx_->basisStatus() == SPxBasis::PRIMAL;
}

/** returns true iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
   *  this does not necessarily meanthat the solver knows and can return the dual ray
   */
bool LPSoplexInterface::ExistsDualRay() {
  MiniMIPdebugMessage("calling ExistsDualRay()\n");

  return (spx_->status() == SPxSolver::INFEASIBLE);
}

/** returns true iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
   *  and the solver knows and can return the dual ray
   */
bool LPSoplexInterface::HasDualRay() {
  MiniMIPdebugMessage("calling HasDualRay()\n");

  return spx_->hasDualFarkas();
}

/** returns true iff LP is proven to be dual unbounded */
bool LPSoplexInterface::IsDualUnbounded() {
  MiniMIPdebugMessage("calling IsDualUnbounded()\n");

  return spx_->status() == SPxSolver::INFEASIBLE && spx_->basisStatus() == SPxBasis::DUAL;
}

/** returns true iff LP is proven to be dual infeasible */
bool LPSoplexInterface::IsDualInfeasible() {
  MiniMIPdebugMessage("calling IsDualInfeasible()\n");

  return (spx_->status() == SPxSolver::UNBOUNDED);
}

/** returns true iff LP is proven to be dual feasible */
bool LPSoplexInterface::IsDualFeasible() {
  MiniMIPdebugMessage("calling IsDualFeasible()\n");

  return (spx_->basisStatus() == SPxBasis::OPTIMAL) || spx_->basisStatus() == SPxBasis::DUAL;
}

/** returns true iff LP was solved to optimality */
bool LPSoplexInterface::IsOptimal() {
  MiniMIPdebugMessage("calling IsOptimal()\n");

  assert((spx_->basisStatus() == SPxBasis::OPTIMAL) == (IsPrimalFeasible() && IsDualFeasible()));

  return (spx_->status() == SPxSolver::OPTIMAL);
}

/** returns true iff current LP solution is stable
   *
   *  This function should return true if the solution is reliable, i.e., feasible and optimal (or proven
   *  infeasible/unbounded) with respect to the original problem. The optimality status might be with respect to a scaled
   *  version of the problem, but the solution might not be feasible to the unscaled original problem; in this case,
   *  MiniMIP::LPInterface.IsStable() should return false.
   */
bool LPSoplexInterface::IsStable() {
  MiniMIPdebugMessage("calling IsStable()\n");

  if (spx_->status() == SPxSolver::ERROR || spx_->status() == SPxSolver::SINGULAR)
    return false;
  if (spx_->status() == SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS)
    return false;
  return true;
}

/** returns true iff the objective limit was reached */
bool LPSoplexInterface::IsObjectiveLimitExceeded() {
  MiniMIPdebugMessage("calling IsObjectiveLimitExceeded()\n");

  return (spx_->status() == SPxSolver::ABORT_VALUE);
}

/** returns true iff the iteration limit was reached */
bool LPSoplexInterface::IsIterationLimitExceeded() {
  MiniMIPdebugMessage("calling IsIterationLimitExceeded()\n");

  return (spx_->status() == SPxSolver::ABORT_ITER);
}

/** returns true iff the time limit was reached */
bool LPSoplexInterface::IsTimeLimitExceeded() {
  MiniMIPdebugMessage("calling IsTimeLimitExceeded()\n");

  return (spx_->status() == SPxSolver::ABORT_TIME);
}

/** gets objective value of solution */
RetCode LPSoplexInterface::GetObjectiveValue(
  LPValue& obj_val /**< the objective value */
) {
  MiniMIPdebugMessage("calling GetObjectiveValue()\n");

  obj_val = spx_->objValueReal();

  return RetCode::kOkay;
}

/** gets primal and dual solution vectors for feasible LPs
   *
   *  Before calling this function, the caller must ensure that the LP has been solved to optimality, i.e., that
   *  MiniMIP::LPInterface.IsOptimal() returns true.
   */
RetCode LPSoplexInterface::GetSolution(
  LPValue& obj_val,          /**< stores the objective value */
  LPValueArray& primal_sol,  /**< primal solution vector */
  LPValueArray& dual_sol,    /**< dual solution vector */
  LPValueArray& activity,    /**< row activity vector */
  LPValueArray& reduced_cost /**< reduced cost vector */
) {
  MiniMIPdebugMessage("calling GetSolution()\n");
  obj_val = spx_->objValueReal();

  try {
    /*    if (primal_sol != 0) { */
    static_cast<void>(spx_->getPrimalReal(primal_sol.data(), spx_->numColsReal()));
    /*    }
    if (dual_sol.size() != 0) {*/
    static_cast<void>(spx_->getDualReal(dual_sol.data(), spx_->numRowsReal()));
    /*    }
    if (activity.size() != 0) {*/
    static_cast<void>(spx_->getSlacksReal(activity.data(), spx_->numRowsReal())); /* in SoPlex, the activities are called "slacks" */
                                                                                  /*    }
    if (reduced_cost.size() != 0) {*/
    static_cast<void>(spx_->getRedCostReal(reduced_cost.data(), spx_->numColsReal()));
    /*    }*/
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    /* MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception: %s\n", s.c_str()); */
#else
  catch (const SPxException&) {
#endif
    return RetCode::kLPError;
  }

  return RetCode::kOkay;
}

/** gets primal ray for unbounded LPs */
RetCode LPSoplexInterface::GetPrimalRay(
  LPValueArray& primal_ray /**< primal ray */
) {
  MiniMIPdebugMessage("calling GetPrimalRay()\n");

  assert(spx_->hasPrimalRay());

  try {
    static_cast<void>(spx_->getPrimalRayReal(primal_ray.data(), spx_->numColsReal()));
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    /* MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception: %s\n", s.c_str()); */
#else
  catch (const SPxException&) {
#endif
    return RetCode::kLPError;
  }
  return RetCode::kOkay;
}

/** gets dual Farkas proof for infeasibility */
RetCode LPSoplexInterface::GetDualFarkasMultiplier(
  LPValueArray& dual_farkas_multiplier /**< dual Farkas row multipliers */
) {
  MiniMIPdebugMessage("calling GetDualFarkasMultiplier()\n");

  assert(spx_->hasDualFarkas());

  try {
    static_cast<void>(spx_->getDualFarkasReal(dual_farkas_multiplier.data(), spx_->numRowsReal()));
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    /* MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception: %s\n", s.c_str()); */
#else
  catch (const SPxException&) {
#endif
    return RetCode::kLPError;
  }
  return RetCode::kOkay;
}

/** gets the number of LP iterations of the last solve call */
RetCode LPSoplexInterface::GetIterations(
  LPNum& iterations /**< the number of iterations of the last solve call */
) {
  MiniMIPdebugMessage("calling GetIterations()\n");

  iterations = spx_->numIterations();

  return RetCode::kOkay;
}

/**@} */

/**@name LP Basis Methods */
/**@{ */

/** gets current basis status for columns and rows */
RetCode LPSoplexInterface::GetBase(
  LPBaseStatArray& column_basis_status, /**< array to store column basis status */
  LPBaseStatArray& row_basis_status     /**< array to store row basis status */
) {
  int i;

  MiniMIPdebugMessage("calling GetBase()\n");
  assert(PreStrongBranchingBasisFreed());

  if (row_basis_status.size() != 0) {
    for (i = 0; i < spx_->numRowsReal(); ++i) {
      switch (spx_->basisRowStatus(i)) {
        case SPxSolver::BASIC:
          row_basis_status[i] = LPBaseStat::kBasic; /*lint !e641*/
          break;
        case SPxSolver::FIXED:
        case SPxSolver::ON_LOWER:
          row_basis_status[i] = LPBaseStat::kLower; /*lint !e641*/
          break;
        case SPxSolver::ON_UPPER:
          row_basis_status[i] = LPBaseStat::kUpper; /*lint !e641*/
          break;
        case SPxSolver::ZERO:
          MiniMIPerrorMessage("slack variable has basis status ZERO (should not occur)\n");
          return RetCode::kLPError;
        case SPxSolver::UNDEFINED:
        default:
          MiniMIPerrorMessage("invalid basis status\n");
          /* std::abort(); */
          return RetCode::kInvalidData; /*lint !e527*/
      }
    }
  }

  if (column_basis_status.size() != 0) {
    for (i = 0; i < spx_->numColsReal(); ++i) {
      /*         LPValue obj_coeffs = 0.0; */
      switch (spx_->basisColStatus(i)) {
        case SPxSolver::BASIC:
          column_basis_status[i] = LPBaseStat::kBasic; /*lint !e641*/
          break;
        case SPxSolver::FIXED:
          /* Get reduced cost estimation. If the estimation is not correct this should not hurt:
           * If the basis is loaded into SoPlex again, the status is converted to FIXED again; in
           * this case there is no problem at all. If the basis is saved and/or used in some other
           * solver, it usually is very cheap to perform the pivots necessary to get an optimal
           * basis.
           * @todo implement getRedCostEst()
           *
           * MINIMIP_CALL( getRedCostEst(spx_, i, &obj_coeffs) );
           *             if( obj_coeffs < 0.0 )  // reduced costs < 0 => UPPER  else => LOWER
           *                column_basis_status[i] = BASESTAT_UPPER;  //lint !e641
           *             else */
          column_basis_status[i] = LPBaseStat::kLower; /*lint !e641*/
          break;
        case SPxSolver::ON_LOWER:
          column_basis_status[i] = LPBaseStat::kLower; /*lint !e641*/
          break;
        case SPxSolver::ON_UPPER:
          column_basis_status[i] = LPBaseStat::kUpper; /*lint !e641*/
          break;
        case SPxSolver::ZERO:
          column_basis_status[i] = LPBaseStat::kZero; /*lint !e641*/
          break;
        case SPxSolver::UNDEFINED:
        default:
          MiniMIPerrorMessage("invalid basis status\n");
          /* MINIMIP_ABORT(); */
          assert(false);
          return RetCode::kInvalidData; /*lint !e527*/
      }
    }
  }
  return RetCode::kOkay;
}

/** sets current basis status for columns and rows */
RetCode LPSoplexInterface::SetBase(
  const LPBaseStatArray& column_basis_status, /**< array with column basis status */
  const LPBaseStatArray& row_basis_status     /**< array with row basis status */
) {
  int i;

  MiniMIPdebugMessage("calling SetBase()\n");

  LPNum num_cols = GetNumberOfColumns();
  LPNum num_rows = GetNumberOfRows();

  assert(PreStrongBranchingBasisFreed());
  InvalidateSolution();

  DataArray<SPxSolver::VarStatus>& _colstat = ColumnsBasisStatus();
  DataArray<SPxSolver::VarStatus>& _rowstat = RowsBasisStatus();

  _colstat.reSize(num_cols);
  _rowstat.reSize(num_rows);

  for (i = 0; i < static_cast<int>(num_rows); ++i) {
    switch (row_basis_status[i]) /*lint !e613*/
    {
      case LPBaseStat::kLower:
        _rowstat[i] = SPxSolver::ON_LOWER;
        break;
      case LPBaseStat::kBasic:
        _rowstat[i] = SPxSolver::BASIC;
        break;
      case LPBaseStat::kUpper:
        _rowstat[i] = SPxSolver::ON_UPPER;
        break;
      case LPBaseStat::kZero:
        MiniMIPerrorMessage("slack variable has basis status ZERO (should not occur)\n");
        return RetCode::kLPError; /*lint !e429*/
      default:
        MiniMIPerrorMessage("invalid basis status\n");

        /* MINIMIP_ABORT(); */
        assert(false);

        return RetCode::kInvalidData; /*lint !e527*/
    }
  }

  for (i = 0; i < static_cast<int>(num_cols); ++i) {
    switch (column_basis_status[i]) /*lint !e613*/
    {
      case LPBaseStat::kLower:
        _colstat[i] = SPxSolver::ON_LOWER;
        break;
      case LPBaseStat::kBasic:
        _colstat[i] = SPxSolver::BASIC;
        break;
      case LPBaseStat::kUpper:
        _colstat[i] = SPxSolver::ON_UPPER;
        break;
      case LPBaseStat::kZero:
        _colstat[i] = SPxSolver::ZERO;
        break;
      default:
        MiniMIPerrorMessage("invalid basis status\n");
        /* MINIMIP_ABORT(); */
        assert(false);
        return RetCode::kInvalidData; /*lint !e527*/
    }
  }

  SOPLEX_TRY(spx_->setBasis(_rowstat.get_ptr(), _colstat.get_ptr()));
  FreePreStrongBranchingBasis();

  return RetCode::kOkay;
}

/** returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m */
RetCode LPSoplexInterface::GetBasisIndices(
  IntArray& basis_indices /**< array to store basis indices ready to keep number of rows entries */
) {

  MiniMIPdebugMessage("calling GetBasisInd()\n");

  assert(PreStrongBranchingBasisFreed());

  spx_->getBasisInd(basis_indices.data());

  return RetCode::kOkay;
}

/** get row of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
RetCode LPSoplexInterface::GetBInvertedRow(
  LPIndex row_number,         /**< row number */
  LPValueArray& row_coeffs, /**< array to store the coefficients of the row */
  LPIndexArray& indices,    /**< array to store the non-zero indices */
  int& num_indices          /**< the number of non-zero indices (-1: if we do not store sparsity information) */
) {
  MiniMIPdebugMessage("calling GetBInvertedRow()\n");

  assert(PreStrongBranchingBasisFreed());
  assert(row_number >= 0);
  assert(row_number < (unsigned) spx_->numRowsReal());

  IntArray integer_indices(indices.begin(), indices.end());

  if (!spx_->getBasisInverseRowReal(row_number, row_coeffs.data(), integer_indices.data(), &num_indices))
    return RetCode::kLPError;

  indices.assign(integer_indices.begin(), integer_indices.end());

  return RetCode::kOkay;
}

/** get column of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
RetCode LPSoplexInterface::GetBInvertedColumn(
  LPIndex col_number,         /**< column number of B^-1; this is NOT the number of the column in the LP;
                             *   you have to call MiniMIP::LPInterface.GetBasisIndices() to get the array which links the
                             *   B^-1 column numbers to the row and column numbers of the LP!
                             *   c must be between 0 and num_rows-1, since the basis has the size
                             *   num_rows * num_rows */
  LPValueArray& col_coeffs, /**< array to store the coefficients of the column */
  LPIndexArray& indices,    /**< array to store the non-zero indices */
  int& num_indices          /**< the number of non-zero indices (-1: if we do not store sparsity information) */
) {
  MiniMIPdebugMessage("calling GetBInvertedColumn()\n");

  assert(PreStrongBranchingBasisFreed());

  IntArray integer_indices(indices.begin(), indices.end());

  if (!spx_->getBasisInverseColReal(col_number, col_coeffs.data(), integer_indices.data(), &num_indices))
    return RetCode::kLPError;

  indices.assign(integer_indices.begin(), integer_indices.end());

  return RetCode::kOkay;
}

/** get row of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
RetCode LPSoplexInterface::GetBInvertedARow(
  LPIndex row_number,                   /**< row number */
  const LPValueArray& b_inverted_row, /**< row in (A_B)^-1 from prior call to MiniMIP::LPInterface.GetBInvRow() */
  LPValueArray& row_coeffs,           /**< array to store coefficients of the row */
  LPIndexArray& indices,              /**< array to store the non-zero indices */
  int& num_indices                    /**< thee number of non-zero indices (-1: if we do not store sparsity information) */
) {

  LPValueArray buf;
  LPValueArray binv;
  int num_rows;
  int num_cols;
  int c;

  MiniMIPdebugMessage("calling GetBInvertedARow()\n");

  assert(PreStrongBranchingBasisFreed());

  num_rows = spx_->numRowsReal();
  num_cols = spx_->numColsReal();

  buf.resize(num_rows);

  /* get (or calculate) the row in B^-1 */
  if (b_inverted_row.size() == 0) {
    MINIMIP_CALL(GetBInvertedRow(row_number, buf, indices, num_indices));
    binv.assign(buf.begin(), buf.end());
  } else
    binv.assign(b_inverted_row.begin(), b_inverted_row.end());

  assert(binv.size() == (unsigned) num_rows);

  /* mark sparsity pattern as invalid */
  num_indices = -1;

  /* @todo exploit sparsity in binv by looping over num_rows */
  /* calculate the scalar product of the row in B^-1 and A */
  Vector binv_vec(num_rows, binv.data());

  /* temporary unscaled column of A */
  DSVector acol;

  for (c = 0; c < num_cols; ++c) {
    spx_->getColVectorReal(c, acol);
    row_coeffs[c] = binv_vec * acol; /* scalar product */ /*lint !e1702*/
  }

  return RetCode::kOkay;
}

/** get column of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
RetCode LPSoplexInterface::GetBInvertedAColumn(
  LPIndex col_number,         /**< column number */
  LPValueArray& col_coeffs, /**< array to store coefficients of the column */
  LPIndexArray& indices,    /**< array to store the non-zero indices */
  int& num_indices          /**< the number of non-zero indices (-1: if we do not store sparsity information) */
) {
  /* create a new uninitialized full vector */
  DVector col(spx_->numRowsReal());

  /* temporary sparse vector used for unscaling (memory is automatically enlarged) */
  DSVector colsparse;

  MiniMIPdebugMessage("calling GetBInvertedAColumn()\n");

  assert(PreStrongBranchingBasisFreed());

  /* extract column col_number of A */
  assert(col_number >= 0);
  assert(col_number < (unsigned) spx_->numColsReal());

  /* @todo implement this with sparse vectors */
  /* mark sparsity pattern as invalid */
  num_indices = -1;

  /* col needs to be cleared because copying colVectorReal only regards nonzeros */
  col.clear();

  spx_->getColVectorReal(col_number, colsparse);
  /* the copy is necessary to transform the sparse column into a dense vector */
  col = colsparse;

  /* solve */
  if (!spx_->getBasisInverseTimesVecReal(col.get_ptr(), col_coeffs.data()))
    return RetCode::kLPError;

  return RetCode::kOkay;
}

/**@} */

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of LP */
RetCode LPSoplexInterface::GetIntegerParameter(
  LPParameter type, /**< parameter number */
  LPNum& param_val  /**< returns the parameter value */
) {
  int scale_param;

  MiniMIPdebugMessage("calling GetIntegerParameter()\n");

  switch (type) {
    case LPParameter::kFromScratch:
      param_val = GetFromScratch();
      break;
    case LPParameter::kLPInfo:
      param_val = GetLPInfo();
      break;
    case LPParameter::kLPIterationLimit:
      if (spx_->intParam(SoPlex::ITERLIMIT) == -1)
        param_val = INT_MAX;
      break;
      param_val = spx_->intParam(SoPlex::ITERLIMIT);
    case LPParameter::kPresolving:
      param_val = spx_->intParam(SoPlex::SIMPLIFIER) == SoPlex::SIMPLIFIER_AUTO;
      break;
    case LPParameter::kPricing:
      param_val = (LPNum) pricing_;
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
      param_val = (LPNum) (spx_->intParam(SoPlex::TIMER));
      break;
    case LPParameter::kRandomSeed:
      param_val = (LPNum) spx_->randomSeed();
      break;
    case LPParameter::kRefactor:
      param_val = (LPNum) spx_->intParam(SoPlex::FACTOR_UPDATE_MAX);
      break;
    default:
      return RetCode::kParameterUnknown;
  } /*lint !e788*/

  return RetCode::kOkay;
}

/** sets integer parameter of LP */
RetCode LPSoplexInterface::SetIntegerParameter(
  LPParameter type, /**< parameter number */
  LPNum param_val   /**< parameter value */
) {
  MiniMIPdebugMessage("calling SetIntegerParameter()\n");

  switch (type) {
    case LPParameter::kFromScratch:
      assert(param_val == true || param_val == false);
      SetFromScratch(bool(param_val));
      break;
    case LPParameter::kLPInfo:
      assert(param_val == true || param_val == false);
      SetLPInfo(bool(param_val));
      break;
    case LPParameter::kLPIterationLimit:
      assert(param_val >= 0);
      /* -1 <= param_val, -1 meaning no time limit, 0 stopping immediately */
      if (param_val >= INT_MAX)
        static_cast<void>(spx_->setIntParam(SoPlex::ITERLIMIT, -1));
      break;
    case LPParameter::kPresolving:
      assert(param_val == true || param_val == false);
      static_cast<void>(spx_->setIntParam(SoPlex::SIMPLIFIER, (param_val ? SoPlex::SIMPLIFIER_AUTO : SoPlex::SIMPLIFIER_OFF)));
      break;
    case LPParameter::kPricing:
      pricing_ = (LPPricing) param_val;
      switch (pricing_) {
        case LPPricing::kDefault:
        case LPPricing::kAuto:
          static_cast<void>(spx_->setIntParam(SoPlex::PRICER, SoPlex::PRICER_AUTO));
          break;
        case LPPricing::kFull:
          static_cast<void>(spx_->setIntParam(SoPlex::PRICER, SoPlex::PRICER_STEEP));
          break;
        case LPPricing::kPartial:
          static_cast<void>(spx_->setIntParam(SoPlex::PRICER, SoPlex::PRICER_PARMULT));
          break;
        case LPPricing::kSteep:
          static_cast<void>(spx_->setIntParam(SoPlex::PRICER, SoPlex::PRICER_STEEP));
          break;
        case LPPricing::kSteepQStart:
          static_cast<void>(spx_->setIntParam(SoPlex::PRICER, SoPlex::PRICER_QUICKSTEEP));
          break;
        case LPPricing::kDevex:
          static_cast<void>(spx_->setIntParam(SoPlex::PRICER, SoPlex::PRICER_DEVEX));
          break;
        default:
          return RetCode::kLPError;
      }
      break;
    case LPParameter::kScaling:
      assert(param_val >= 0 && param_val <= 2);
      if (param_val == 0)
        static_cast<void>(spx_->setIntParam(SoPlex::SCALER, SoPlex::SCALER_OFF));
      else if (param_val == 1)
        static_cast<void>(spx_->setIntParam(SoPlex::SCALER, SoPlex::SCALER_BIEQUI));
      else
        static_cast<void>(spx_->setIntParam(SoPlex::SCALER, SoPlex::SCALER_LEASTSQ));
      break;
    case LPParameter::kTiming:
      assert(param_val >= 0 && param_val < 3);
      static_cast<void>(spx_->setIntParam(SoPlex::TIMER, param_val));
      break;
    case LPParameter::kRandomSeed:
      spx_->setRandomSeed((unsigned long) (long) param_val); /* TODO: add c++ style static_cast */
      break;
    case LPParameter::kRefactor:
      assert(param_val >= 0);
      static_cast<void>(spx_->setIntParam(SoPlex::FACTOR_UPDATE_MAX, param_val));
      break;
    default:
      return RetCode::kParameterUnknown;
  } /*lint !e788*/

  return RetCode::kOkay;
}

/** gets floating point parameter of LP */
RetCode LPSoplexInterface::GetRealParameter(
  LPParameter type,    /**< parameter number */
  LPValue& LPValue_val /**< returns the parameter value */
) {
  MiniMIPdebugMessage("calling GetRealParameter()\n");

  switch (type) {
    case LPParameter::kFeasibilityTolerance:
      LPValue_val = FeasibilityTolerance();
      break;
    case LPParameter::kDualFeasibilityTolerance:
      LPValue_val = OptimalityTolerance();
      break;
    case LPParameter::kObjectiveLimit:
      if (spx_->intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE)
        LPValue_val = spx_->realParam(SoPlex::OBJLIMIT_UPPER);
      else
        LPValue_val = spx_->realParam(SoPlex::OBJLIMIT_LOWER);
      break;
    case LPParameter::kLPTimeLimit:
      LPValue_val = spx_->realParam(SoPlex::TIMELIMIT);
      break;
    case LPParameter::kMarkowitz:
      LPValue_val = spx_->realParam(SoPlex::MIN_MARKOWITZ);
      break;
    default:
      return RetCode::kParameterUnknown;
  } /*lint !e788*/

  return RetCode::kOkay;
}

/** sets floating point parameter of LP */
RetCode LPSoplexInterface::SetRealParameter(
  LPParameter type,   /**< parameter number */
  LPValue LPValue_val /**< parameter value */
) {
  MiniMIPdebugMessage("calling SetRealParameter()\n");

  switch (type) {
    case LPParameter::kFeasibilityTolerance:
      /* 0 < LPValue_val */
      assert(LPValue_val > 0.0);
      SetFeasibilityTolerance(LPValue_val);
      break;
    case LPParameter::kDualFeasibilityTolerance:
      /* 0 < LPValue_val */
      assert(LPValue_val > 0.0);
      SetOptimalityTolerance(LPValue_val);
      break;
    case LPParameter::kObjectiveLimit:
      /* no restrictions on LPValue_val */
      if (spx_->intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE)
        static_cast<void>(spx_->setRealParam(SoPlex::OBJLIMIT_UPPER, LPValue_val));
      else
        static_cast<void>(spx_->setRealParam(SoPlex::OBJLIMIT_LOWER, LPValue_val));
      break;
    case LPParameter::kLPTimeLimit:
      assert(LPValue_val > 0.0);
      /* soplex requires 0 < LPValue_val < DEFAULT_INFINITY (= 1e100), -1 means unlimited */
      static_cast<void>(spx_->setRealParam(SoPlex::TIMELIMIT, LPValue_val));
      break;
    case LPParameter::kMarkowitz:
      /* 1e-4 <= LPValue_val <= 0.999 */
      if (LPValue_val < 1e-4)
        LPValue_val = 1e-4;
      else if (LPValue_val > 0.9999)
        LPValue_val = 0.9999;

      static_cast<void>(spx_->setRealParam(SoPlex::MIN_MARKOWITZ, LPValue_val));
      break;
    default:
      return RetCode::kParameterUnknown;
  } /*lint !e788*/

  return RetCode::kOkay;
}

/**@} */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the LP solver */
LPValue LPSoplexInterface::Infinity() {
  MiniMIPdebugMessage("calling Infinity()\n");

  return spx_->realParam(SoPlex::INFTY);
}

/** checks if given value is treated as infinity in the LP solver */
bool LPSoplexInterface::IsInfinity(
  LPValue val /**< value to be checked for infinity */
) {
  MiniMIPdebugMessage("calling IsInfinity()\n");

  return (val >= spx_->realParam(SoPlex::INFTY));
}

/** returns, whether the given file exists */
bool LPSoplexInterface::FileExists(
  const char* file_name /**< file name */
) {
  FILE* file;

  file = fopen(file_name, "r");
  if (file == NULL)
    return false;

  fclose(file);

  return true;
}

/**@} */

/**@name File Interface Methods */
/**@{ */

/** reads LP from a file */
RetCode LPSoplexInterface::ReadLP(
  const char* file_name /**< file name */
) {
  MiniMIPdebugMessage("calling ReadLP()\n");

  assert(file_name != NULL);

  assert(PreStrongBranchingBasisFreed());

  if (!FileExists(file_name))
    return RetCode::kReadError;

  try {
    assert(spx_->intParam(SoPlex::READMODE) == SoPlex::READMODE_REAL);
    if (!spx_->readFile((char*) (file_name)))
      return RetCode::kReadError;
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    /* MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception: %s\n", s.c_str()); */
#else
  catch (const SPxException&) {
#endif
    return RetCode::kReadError;
  }

  return RetCode::kOkay;
}

/** writes LP to a file */
RetCode LPSoplexInterface::WriteLP(
  const char* file_name /**< file name */
) {
  MiniMIPdebugMessage("calling WriteLP()\n");

  assert(file_name != NULL);

  try {
    static_cast<void>(spx_->writeFileReal(file_name));
  }
#ifndef NDEBUG
  catch (const SPxException& x) {
    std::string s = x.what();
    /* MiniMIPmessagePrintWarning(messagehdlr, "SoPlex threw an exception: %s\n", s.c_str()); */
#else
  catch (const SPxException&) {
#endif
    return RetCode::kWriteError;
  }
  return RetCode::kOkay;
}

} /* namespace minimip*/
