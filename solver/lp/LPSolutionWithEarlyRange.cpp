#include "LPSolutionWithEarlyRange.h"

#include <utility>

int LPSolutionWithEarlyRange::getEarlyRangeLower() const {
    return earlyRangeLower;
}

int LPSolutionWithEarlyRange::getEarlyRangeUpper() const {
    return earlyRangeUpper;
}

LPSolutionWithEarlyRange::LPSolutionWithEarlyRange(int earlyRangeLower, int earlyRangeUpper) : earlyRangeLower(
        earlyRangeLower), earlyRangeUpper(earlyRangeUpper) {}

LPSolutionWithEarlyRange::LPSolutionWithEarlyRange(double objective, const std::vector<double> &variables,
                                                   int earlyRangeLower, int earlyRangeUpper) : LPSolution(objective,
                                                                                                          variables),
                                                                                               earlyRangeLower(
                                                                                                       earlyRangeLower),
                                                                                               earlyRangeUpper(
                                                                                                       earlyRangeUpper) {}

LPSolutionWithEarlyRange::LPSolutionWithEarlyRange(double objective, const std::vector<double> &variables,
                                                   const std::vector<double> &reducedCosts, int earlyRangeLower,
                                                   int earlyRangeUpper) : LPSolution(objective, variables,
                                                                                     reducedCosts),
                                                                          earlyRangeLower(earlyRangeLower),
                                                                          earlyRangeUpper(earlyRangeUpper) {}

LPSolutionWithEarlyRange::LPSolutionWithEarlyRange(double objective, const std::vector<double> &variables,
                                                   const std::vector<double> &reducedCosts,
                                                   const std::vector<double> &shadowPrices, int earlyRangeLower,
                                                   int earlyRangeUpper) : LPSolution(objective, variables, reducedCosts,
                                                                                     shadowPrices),
                                                                          earlyRangeLower(earlyRangeLower),
                                                                          earlyRangeUpper(earlyRangeUpper) {}
