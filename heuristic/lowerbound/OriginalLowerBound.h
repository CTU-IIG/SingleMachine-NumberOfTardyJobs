#ifndef ORIGINALLOWERBOUND_H
#define ORIGINALLOWERBOUND_H


#include <vector>
#include "../../solver/ilp/BinaryILPSolution.h"
#include "../../scheduling/SchedulingProblem.h"

class OriginalLowerBound {
private:
    ILPSolverImpl ilpSolverImpl = ILPSolverImpl::ORIGINAL;
    LowerBoundImpl lowerBoundImpl = LowerBoundImpl::ORIGINAL;
    UpperBoundImpl upperBoundImpl = UpperBoundImpl::ORIGINAL;
    int upperBoundMaxDepth = 4;
public:
    explicit OriginalLowerBound(ILPSolverImpl ilpSolverImpl, LowerBoundImpl lowerBoundImpl, UpperBoundImpl upperBoundImpl, int upperBoundMaxDepth);

    OriginalLowerBound() = delete;

    virtual ~OriginalLowerBound() = default;

    virtual BinaryILPSolution solve(const SchedulingProblem &schedProblem, const LPSolution &lpSolution) const;

    virtual std::vector<int>
    computeCoreProblem(const SchedulingProblem &schedProblem, const LPSolution &lpSolution) const;

    virtual BinaryILPSolution
    localSearch(const SchedulingProblem &schedProblem, const BinaryILPSolution &solution) const;

private:

    std::vector<std::pair<int, bool>>
    getReduceJobsFromCoreSet(const SchedulingProblem &schedProblem, const LPSolution &lpSolution,
                             const std::vector<int> &coreSet) const;
};


#endif //ORIGINALLOWERBOUND_H
