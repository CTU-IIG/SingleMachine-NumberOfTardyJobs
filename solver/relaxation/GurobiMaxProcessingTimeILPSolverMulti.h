#ifndef GUROBIMAXPROCESSINGTIMEILPSOLVERMULTI_H
#define GUROBIMAXPROCESSINGTIMEILPSOLVERMULTI_H


#include "../ilp/BinaryILPSolution.h"
#include "../../scheduling/SchedulingProblem.h"
#include "../ilp/ILPSolver.h"

class GurobiMaxProcessingTimeILPSolverMulti : public ILPSolver {
private:
    uint32_t threadCount = 1;
    bool printOutput = false;
    bool usePresolve = false;
public:
    explicit GurobiMaxProcessingTimeILPSolverMulti(uint32_t threadCount, bool printOutput, bool usePresolve);

    GurobiMaxProcessingTimeILPSolverMulti() = delete;

    virtual ~GurobiMaxProcessingTimeILPSolverMulti() = default;

    BinaryILPSolution solve(const SchedulingProblem &schedProblem) const override;

    BinaryILPSolution solve(const SchedulingProblem &schedProblem, const BinaryILPSolution &lowerBound) const override;

private:
    BinaryILPSolution solve_impl(const SchedulingProblem &schedProblem, const BinaryILPSolution &lowerBound) const;
};


#endif //GUROBIMAXPROCESSINGTIMEILPSOLVERMULTI_H
