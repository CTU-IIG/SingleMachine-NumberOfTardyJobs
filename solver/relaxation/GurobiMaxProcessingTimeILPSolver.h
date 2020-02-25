#ifndef GUROBIMAXPROCESSINGTIMEILPSOLVER_H
#define GUROBIMAXPROCESSINGTIMEILPSOLVER_H


#include "../ilp/ILPSolver.h"

class GurobiMaxProcessingTimeILPSolver : public ILPSolver {
private:
    uint32_t threadCount = 1;
    bool printOutput = false;
    bool usePresolve = false;
public:
    explicit GurobiMaxProcessingTimeILPSolver(uint32_t threadCount, bool printOutput, bool usePresolve);

    GurobiMaxProcessingTimeILPSolver() = delete;

    virtual ~GurobiMaxProcessingTimeILPSolver() = default;

    BinaryILPSolution solve(const SchedulingProblem &schedProblem) const override;

    BinaryILPSolution solve(const SchedulingProblem &schedProblem, const BinaryILPSolution &lowerBound) const override;

private:
    BinaryILPSolution solve_impl(const SchedulingProblem &schedProblem, const BinaryILPSolution &lowerBound) const;
};


#endif //GUROBIMAXPROCESSINGTIMEILPSOLVER_H
