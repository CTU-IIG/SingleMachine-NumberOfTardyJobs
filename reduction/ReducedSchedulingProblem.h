#ifndef REDUCEDSCHEDULINGPROBLEM_H
#define REDUCEDSCHEDULINGPROBLEM_H


#include <vector>
#include <unordered_map>
#include "../scheduling/SchedulingProblem.h"

class ReducedSchedulingProblem {
private:
    const SchedulingProblem originalSchedProblem;
    const SchedulingProblem reducedSchedProblem;
    std::vector<std::pair<int, bool>> reduceJobs;
    mutable int64_t reducedJobsObjective = -1;
    mutable std::unordered_map<int, int> indexMapping;

    int64_t computeReducedJobsObjective() const;

public:
    explicit ReducedSchedulingProblem(const SchedulingProblem &originalSchedProblem,
                                      const SchedulingProblem &reducedSchedProblem,
                                      const std::vector<std::pair<int, bool>> &reduceJobs);

    ReducedSchedulingProblem() = delete;

    virtual ~ReducedSchedulingProblem() = default;

    const std::vector<job_t> &getJobs() const;

    bool empty() const;

    size_t size() const;

    const job_t &operator[](const int index) const;

    int64_t getReducedJobsObjective() const;

    BinaryILPSolution reconstructSchedulingProblem(const BinaryILPSolution &reducedSolution) const;

    SchedulingProblem getReducedSchedulingProblemCopy() const;

    bool isFeasible() const;

    std::vector<int> getNonReducedJobs() const;

    std::vector<std::pair<int, bool>> getReduceJobs() const;

    int mapReducedIndexToOriginal(int reducedIndex) const;

    int64_t getNonZeroCoefCountInModel() const;

    BinaryILPSolution transformSolutionToReduced(const BinaryILPSolution &originalSolution) const;
};


#endif //DREDUCEDSCHEDULINGPROBLEM_H
