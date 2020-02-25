#ifndef LPSOLUTIONWITHEARLYRANGE_H
#define LPSOLUTIONWITHEARLYRANGE_H

#include <vector>
#include <cstddef>
#include "LPSolution.h"

class LPSolutionWithEarlyRange : public LPSolution {
private:
    int earlyRangeLower = 0;
    int earlyRangeUpper = 0;

public:
    virtual ~LPSolutionWithEarlyRange() = default;

    LPSolutionWithEarlyRange(int earlyRangeLower, int earlyRangeUpper);

    LPSolutionWithEarlyRange(double objective, const std::vector<double> &variables, int earlyRangeLower,
                             int earlyRangeUpper);

    LPSolutionWithEarlyRange(double objective, const std::vector<double> &variables,
                             const std::vector<double> &reducedCosts, int earlyRangeLower, int earlyRangeUpper);

    LPSolutionWithEarlyRange(double objective, const std::vector<double> &variables,
                             const std::vector<double> &reducedCosts, const std::vector<double> &shadowPrices,
                             int earlyRangeLower, int earlyRangeUpper);

    int getEarlyRangeLower() const;

    int getEarlyRangeUpper() const;
};


#endif //LPSOLUTIONWITHEARLYRANGE_H
