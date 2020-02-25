#ifndef JOBREDUCE_H
#define JOBREDUCE_H


#include "../scheduling/SchedulingProblem.h"
#include "ReducedSchedulingProblem.h"
#include "../solver/lp/LPSolutionWithBranchPenalty.h"

class JobReduce {
private:
    ILPSolverImpl ilpSolverImpl = ILPSolverImpl::ORIGINAL;
    LowerBoundImpl lowerBoundImpl = LowerBoundImpl::ORIGINAL;
    UpperBoundImpl upperBoundImpl = UpperBoundImpl::ORIGINAL;
    int upperBoundMaxDepth = 4;
public:
    explicit JobReduce(ILPSolverImpl ilpSolverImpl, LowerBoundImpl lowerBoundImpl, UpperBoundImpl upperBoundImpl, int upperBoundMaxDepth);

    JobReduce() = delete;

    SchedulingProblem reduceJob_impl(const SchedulingProblem &schedProblem, int jobIndex, bool jobEarly) const;
    SchedulingProblem reduceJobs_impl(const SchedulingProblem &schedProblem, const std::vector<std::pair<int, bool>> &reduceJobs) const;
    SchedulingProblem reduceJob_impl(const SchedulingProblem &schedProblem, const std::vector<int> &jobIndexes, const std::vector<bool> &jobEarly) const;

    ReducedSchedulingProblem reduceJob(const SchedulingProblem &schedProblem, int jobIndex, bool jobEarly) const;
    ReducedSchedulingProblem reduceJobs(const SchedulingProblem &schedProblem, const std::vector<std::pair<int, bool>> &reduceJobs, bool allowExtendedReduction = false) const;

    SchedulingProblem reduceJob(const SchedulingProblem &schedProblem, int jobIndex, bool jobEarly, bool useDominanceTheorem) const;

    std::vector<std::pair<int, bool>> getForceJobByDominance(const SchedulingProblem &schedProblem, int jobIndex, bool jobEarly) const;

    std::vector<std::pair<int, bool>> getForceJobByDominance(const SchedulingProblem &schedProblem, const std::vector<std::pair<int, bool>> &reduceJobs) const;

    ReducedSchedulingProblem reduceJobsByVarFixing(const SchedulingProblem &schedProblem, const BinaryILPSolution &ilpSolution, const LPSolution &lpSolution) const;

    std::vector<std::pair<int, bool>> reductionPreprocessing(const SchedulingProblem &schedProblem) const;
};


#endif //JOBREDUCE_H
