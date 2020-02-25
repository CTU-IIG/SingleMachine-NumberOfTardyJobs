#ifndef LPSOLUTIONWITHBRANCHPENALTY_H
#define LPSOLUTIONWITHBRANCHPENALTY_H


#include "LPSolution.h"

typedef struct {
    double upperChange;
    double lowerChange;
} branch_penalty_t;

class LPSolutionWithBranchPenalty : public LPSolution {
private:
    std::vector<std::pair<int, branch_penalty_t>> branchPenalties;

public:
    virtual ~LPSolutionWithBranchPenalty() = default;

    LPSolutionWithBranchPenalty();

    LPSolutionWithBranchPenalty(double objective, const std::vector<double> &variables,
                                const std::vector<std::pair<int, branch_penalty_t>> &branchPenalties);

    LPSolutionWithBranchPenalty(double objective, const std::vector<double> &variables,
                                const std::vector<double> &reducedCosts,
                                const std::vector<std::pair<int, branch_penalty_t>> &branchPenalties);

    LPSolutionWithBranchPenalty(double objective, const std::vector<double> &variables,
                                const std::vector<double> &reducedCosts, const std::vector<double> &shadowPrices,
                                const std::vector<std::pair<int, branch_penalty_t>> &branchPenalties);

    const std::vector<std::pair<int, branch_penalty_t>> &getBranchPenalties() const;

    std::pair<int, double> getMaxMinPseudoCost() const;
};


#endif //LPSOLUTIONWITHBRANCHPENALTY_H
