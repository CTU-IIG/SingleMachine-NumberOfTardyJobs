#ifndef IMPROVEDUPPERBOUND_H
#define IMPROVEDUPPERBOUND_H


#include "../../solver/lp/LPSolver.h"
#include "../../solver/lp/LPSolutionWithBranchPenalty.h"

class ImprovedUpperBound : public LPSolver {
private:
    ILPSolverImpl ilpSolverImpl = ILPSolverImpl::ORIGINAL;
    LowerBoundImpl lowerBoundImpl = LowerBoundImpl::ORIGINAL;
    UpperBoundImpl upperBoundImpl = UpperBoundImpl::ORIGINAL;
    int startDepth = 1;
public:
    ImprovedUpperBound() = delete;

    explicit ImprovedUpperBound(ILPSolverImpl ilpSolverImpl, LowerBoundImpl lowerBoundImpl, UpperBoundImpl upperBoundImpl, int startDepth);

    virtual ~ImprovedUpperBound() = default;

    LPSolution solve(const SchedulingProblem &schedProblem) const override;

    LPSolutionWithBranchPenalty solveWithBranchPenalty(const SchedulingProblem &schedProblem) const override;

private:
    double recursiveSearch(const SchedulingProblem &schedProblem, const LPSolution &lpSolution, int depth = 1, int64_t objectiveSum = 0) const;

    int selectJobWithMostProcessingTimeJob(const SchedulingProblem &schedProblem, const LPSolution &lpSolution) const;
};


#endif //IMPROVEDUPPERBOUND_H
