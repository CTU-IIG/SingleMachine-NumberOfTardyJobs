#ifndef SCHEDULINGPROBLEM_H
#define SCHEDULINGPROBLEM_H


#include <vector>
#include <cstdint>
#include <fstream>
#include "../common/types.h"
#include "../solver/ilp/BinaryILPSolution.h"
#include "../solver/lp/LPSolution.h"


class SchedulingProblem {
private:
    std::vector<job_t> jobs;

    mutable std::vector<int64_t> orderedTimePoints;
    mutable int64_t totalProcessingTime = -1;
    mutable int64_t nonZeroCoefCountInModel = -1;

    int64_t computeTotalProcessingTime() const;

    std::vector<int64_t> computeOrderedTimePoints() const;

    bool containDeadlines = true;

public:
    explicit SchedulingProblem(const std::vector<job_t> &jobs, bool hasDeadlines);

    SchedulingProblem(const SchedulingProblem &schedulingProblem);

    SchedulingProblem() = delete;

    virtual ~SchedulingProblem() = default;

    const std::vector<job_t> &getJobs() const;

    const std::vector<int64_t> &getOrderedTimePoints() const;

    int64_t getTotalProcessingTime() const;

    bool empty() const;

    size_t size() const;

    const job_t &operator[](const int index) const;

    void sortByDeadline();

    void setDeadline(int index, int32_t deadline);

    static bool compareByDeadline(const job_t &first, const job_t &second);

    static SchedulingProblem createFromStream(std::istream &input, bool clampDeadline = false);

    BinaryILPSolution constructILPSolutionFromLPSolution(const LPSolution &lpSolution) const;

    BinaryILPSolution findFeasibleSolutionByEDF() const;

    bool isFeasible() const;

    bool isFeasibleSolution(const std::vector<bool> &solution) const;

    bool isFeasibleSolution(const std::vector<bool> &solution, const std::vector<int> &presorted) const;

    bool isFeasibleSolution(const BinaryILPSolution &solution) const;

    std::vector<int> createScheduleFromSolution(const std::vector<bool> &solution) const;

    std::vector<int> createScheduleFromSolution(const std::vector<bool> &solution, const std::vector<int> &presorted) const;

    int64_t getNonZeroCoefCountInModel() const;

    bool hasDeadlines() const;
};


#endif //SCHEDULINGPROBLEM_H
