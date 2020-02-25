#ifndef BAPTISTE_H
#define BAPTISTE_H


#include <cstdint>
#include "../solver/ilp/BinaryILPSolution.h"
#include "../scheduling/SchedulingProblem.h"

class Baptiste {
private:
    SchedulingProblem schedulingProblem;
    uint32_t threadCount = 1;
    const int64_t ilpModelNonZeroLimit = 14000000;

    volatile int64_t globalBestObjective = 0;
    volatile int64_t totalNodeVisited = 0;
    volatile int64_t rootLowerBound = 0;
    volatile double rootUpperBound = 0;
    volatile int rootReducedProbSize = -1;

    ILPSolverImpl ilpSolverImpl = ILPSolverImpl::ORIGINAL;
    LowerBoundImpl lowerBoundImpl = LowerBoundImpl::ORIGINAL;
    UpperBoundImpl upperBoundImpl = UpperBoundImpl::ORIGINAL;

    int upperBoundMaxDepth = 4;

    BinaryILPSolution branchAndBounds(const SchedulingProblem &schedProblem, int64_t objectiveAdd);

public:
    explicit Baptiste(const SchedulingProblem &schedulingProblem, uint32_t threadCount, int64_t ilpModelNonZeroLimit, ILPSolverImpl ilpSolverImpl, LowerBoundImpl lowerBoundImpl, UpperBoundImpl upperBoundImpl, int upperBoundMaxDepth);

    Baptiste() = delete;

    virtual ~Baptiste() = default;

    BinaryILPSolution solve();

    volatile int64_t getGlobalBestObjective() const;

    volatile int64_t getTotalNodeVisited() const;

    volatile int64_t getRootLowerBound() const;

    volatile double getRootUpperBound() const;

    volatile int getRootReducedProbSize() const;
};


#endif //BAPTISTE_H
