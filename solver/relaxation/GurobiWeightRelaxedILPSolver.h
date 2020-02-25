#ifndef GUROBIWEIGHTRELAXEDILPSOLVER_H
#define GUROBIWEIGHTRELAXEDILPSOLVER_H


#include "../ilp/ILPSolver.h"

class GurobiWeightRelaxedILPSolver : public ILPSolver {
private:
    uint32_t threadCount = 1;
    bool printOutput = false;
    bool usePresolve = false;
public:
    explicit GurobiWeightRelaxedILPSolver(uint32_t threadCount, bool printOutput, bool usePresolve);

    GurobiWeightRelaxedILPSolver() = delete;

    virtual ~GurobiWeightRelaxedILPSolver() = default;

    BinaryILPSolution solve(const SchedulingProblem &schedProblem) const override;

    BinaryILPSolution solve(const SchedulingProblem &schedProblem, const BinaryILPSolution &lowerBound) const override;

private:
    BinaryILPSolution solve_impl(const SchedulingProblem &schedProblem) const;
};


#endif //GUROBIWEIGHTRELAXEDILPSOLVER_H
