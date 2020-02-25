#ifndef GUROBILPSOLVER_H
#define GUROBILPSOLVER_H

#include <cstdint>
#include "LPSolution.h"
#include "../../scheduling/SchedulingProblem.h"
#include "LPSolver.h"

class GurobiLPSolver : public LPSolver {
private:
    uint32_t threadCount = 1;
    bool printOutput = false;
    bool usePresolve = false;
public:
    explicit GurobiLPSolver(uint32_t threadCount, bool printOutput, bool usePresolve);

    GurobiLPSolver() = delete;

    virtual ~GurobiLPSolver() = default;

    LPSolution solve(const SchedulingProblem &schedProblem) const override;

    LPSolutionWithBranchPenalty solveWithBranchPenalty(const SchedulingProblem &schedProblem) const override;

private:
    LPSolution solve_impl(const SchedulingProblem &schedProblem) const;
};


#endif //GUROBILPSOLVER_H
