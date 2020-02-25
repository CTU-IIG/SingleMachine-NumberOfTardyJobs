#include <limits>
#include <cmath>
#include "LPSolutionWithBranchPenalty.h"

const std::vector<std::pair<int, branch_penalty_t>> &LPSolutionWithBranchPenalty::getBranchPenalties() const {
    return branchPenalties;
}

LPSolutionWithBranchPenalty::LPSolutionWithBranchPenalty(double objective, const std::vector<double> &variables,
                                                         const std::vector<std::pair<int, branch_penalty_t>> &branchPenalties)
        : LPSolution(objective, variables), branchPenalties(branchPenalties) {}

LPSolutionWithBranchPenalty::LPSolutionWithBranchPenalty(double objective, const std::vector<double> &variables,
                                                         const std::vector<double> &reducedCosts,
                                                         const std::vector<std::pair<int, branch_penalty_t>> &branchPenalties)
        : LPSolution(objective, variables, reducedCosts), branchPenalties(branchPenalties) {}

LPSolutionWithBranchPenalty::LPSolutionWithBranchPenalty(double objective, const std::vector<double> &variables,
                                                         const std::vector<double> &reducedCosts,
                                                         const std::vector<double> &shadowPrices,
                                                         const std::vector<std::pair<int, branch_penalty_t>> &branchPenalties)
        : LPSolution(objective, variables, reducedCosts, shadowPrices), branchPenalties(branchPenalties) {}

LPSolutionWithBranchPenalty::LPSolutionWithBranchPenalty() {}

std::pair<int, double> LPSolutionWithBranchPenalty::getMaxMinPseudoCost() const {
    double maxMinPseudoCost = std::numeric_limits<double>::lowest();
    int maxMinPseudoCostIndex = -1;

    for (const std::pair<int, branch_penalty_t> &branchPenaltyPair : this->branchPenalties) {
        double minPseudoCost = std::min(std::abs(branchPenaltyPair.second.lowerChange),
                                        std::abs(branchPenaltyPair.second.upperChange));

        if (minPseudoCost > maxMinPseudoCost) {
            maxMinPseudoCost = minPseudoCost;
            maxMinPseudoCostIndex = branchPenaltyPair.first;
        }
    }

    return std::make_pair(maxMinPseudoCostIndex, maxMinPseudoCost);
}
