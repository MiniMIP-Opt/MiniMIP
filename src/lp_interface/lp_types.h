#ifndef MINIMIP_SRC_LP_INTERFACE_LP_TYPES_H_
#define MINIMIP_SRC_LP_INTERFACE_LP_TYPES_H_

#include <string>
#include <vector>

namespace minimip {

typedef unsigned int LPIndex; /**< unsigned integer representing a number or count in LP solvers*/
typedef std::vector<LPIndex> LPIndexArray;

typedef unsigned int LPNum; /**<  unsigned integer representing an index in LP solvers*/
typedef std::vector<LPNum> LPNumArray;

typedef long long LPLongInt;

typedef double LPValue;
typedef std::vector<LPValue> LPValueArray;

typedef std::vector<int> IntArray;
typedef std::vector<double> DoubleArray;
typedef std::vector<std::string> StringArray;
typedef std::vector<bool> BoolArray;

/** return codes for MiniMIP methods: non-positive return codes are errors */
enum class RetCode {
  OKAY = +1,                  /**< normal termination */
  ERROR = 0,                  /**< unspecified error */
  NO_MEMORY = -1,             /**< insufficient memory error */
  READ_ERROR = -2,            /**< read error */
  WRITE_ERROR = -3,           /**< write error */
  LP_ERROR = -4,              /**< error in LP solver */
  NO_PROBLEM = -5,            /**< no problem exists */
  INVALID_CALL = -6,          /**< method cannot be called at this time in solution process */
  INVALID_DATA = -7,          /**< error in input data */
  INVALID_RESULT = -8,        /**< method returned an invalid result code */
  PARAMETER_UNKNOWN = -9,     /**< the parameter with the given name was not found */
  PARAMETER_WRONG_TYPE = -10, /**< the parameter is not of the expected type */
  PARAMETER_WRONG_VAL = -11,  /**< the value is invalid for the given parameter */
  NOT_IMPLEMENTED = -12       /**< function not implemented */
};

/** LP objective sense */
enum class LPObjectiveSense {
  OBJ_SENSE_MAXIMIZE = -1, /**< maximize objective function */
  OBJ_SENSE_MINIMIZE = +1  /**< minimize objective function */
};

/** LP solver parameters */
enum class LPParameter {
  FROM_SCRATCH = 0,               /**< Solver should start from scratch at next call? */
  SCALING = 1,                    /**< Should LP solver use scaling? */
  PRESOLVING = 2,                 /**< Should LP solver use presolving? */
  POLISHING = 3,                  /**< set solution polishing (0: disable, 1: enable) */
  PRICING = 4,                    /**< pricing strategy */
  FEASIBLITY_TOLERANCE = 5,       /**< feasibility tolerance for primal variables and slacks, strictly positive */
  DUAL_FEASIBILITY_TOLERANCE = 6, /**< feasibility tolerance for dual variables and reduced costs, strictly positive */
  MARKOWITZ = 7,                  /**< Markowitz tolerance */
  REFACTOR = 8,                   /**< set refactorization interval (0: automatic) */
  OBJECTIVE_LIMIT = 9,            /**< objective limit (stop if objective is known be larger/smaller than limit for min/max-imization) */
  LP_ITERATION_LIMIT = 10,        /**< LP iteration limit (> 0) */
  LP_TIME_LIMIT = 11,             /**< LP time limit, positive */
  THREADS = 12,                   /**< number of threads used to solve the LP (-1: automatic) */
  TIMING = 13,                    /**< type of timer (1: cpu, 2: wallclock, 0: off) */
  RANDOMSEED = 14,                /**< inital random seed, e.g., for perturbations in the simplex (0: LP default) */
  LP_INFO = 15,                   /**< Should LP solver output information to the screen? */
};

/** LP pricing strategy */
enum class LPPricing {
  DEFAULT = 0,     /**< the MiniMIP/LP interface should use its preferred strategy */
  AUTO = 1,        /**< the LP solver should use its preferred strategy */
  FULL = 2,        /**< full pricing */
  PARTIAL = 3,     /**< partial pricing */
  STEEP = 4,       /**< steepest edge pricing */
  STEEPQSTART = 5, /**< steepest edge pricing without initial dual norms */
  DEVEX = 6        /**< devex pricing */
};

/** LP basis status for columns and rows */
enum class LPBaseStat {
  BASESTAT_LOWER = 0, /**< (slack) variable is at its lower bound */
  BASESTAT_BASIC = 1, /**< (slack) variable is basic */
  BASESTAT_UPPER = 2, /**< (slack) variable is at its upper bound */
  BASESTAT_ZERO = 3   /**< free variable is non-basic and set to zero */
};

typedef std::vector<LPBaseStat> LPBaseStatArray;

}/* namespace minimip */
#endif /* MINIMIP_SRC_LP_INTERFACE_LP_TYPES_H_ */
