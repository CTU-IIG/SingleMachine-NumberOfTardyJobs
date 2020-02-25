#ifndef IMPROVEDLOWERBOUND_H
#define IMPROVEDLOWERBOUND_H


#include "OriginalLowerBound.h"

class ImprovedLowerBound : public OriginalLowerBound {
public:
    explicit ImprovedLowerBound(ILPSolverImpl ilpSolverImpl, LowerBoundImpl lowerBoundImpl, UpperBoundImpl upperBoundImpl, int upperBoundMaxDepth);

    ImprovedLowerBound() = delete;

    virtual ~ImprovedLowerBound() = default;

    virtual std::vector<int>
    computeCoreProblem(const SchedulingProblem &schedProblem, const LPSolution &lpSolution) const override;

    virtual BinaryILPSolution
    localSearch(const SchedulingProblem &schedProblem, const BinaryILPSolution &solution) const override;
};


#endif //IMPROVEDLOWERBOUND_H
