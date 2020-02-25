#ifndef LPSOLVER_H
#define LPSOLVER_H


#include "../../scheduling/SchedulingProblem.h"
#include "LPSolution.h"
#include "LPSolutionWithBranchPenalty.h"

class LPSolver {
public:
    LPSolver() = default;

    virtual ~LPSolver() = default;

    virtual LPSolution solve(const SchedulingProblem &schedProblem) const = 0;

    virtual LPSolutionWithBranchPenalty solveWithBranchPenalty(const SchedulingProblem &schedProblem) const = 0;
};


#endif //LPSOLVER_H
