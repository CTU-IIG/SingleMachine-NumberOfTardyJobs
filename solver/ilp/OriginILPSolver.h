#ifndef ORIGINILPSOLVER_H
#define ORIGINILPSOLVER_H


#include "BinaryILPSolution.h"
#include "../../scheduling/SchedulingProblem.h"
#include "ILPSolver.h"

class OriginILPSolver : public ILPSolver {
private:
    uint32_t threadCount = 1;
    bool printOutput = false;
    bool usePresolve = false;
public:
    explicit OriginILPSolver(uint32_t threadCount, bool printOutput, bool usePresolve);

    OriginILPSolver() = delete;

    virtual ~OriginILPSolver() = default;


    BinaryILPSolution solve(const SchedulingProblem &schedProblem) const override;

    BinaryILPSolution solve(const SchedulingProblem &schedProblem, const BinaryILPSolution &lowerBound) const override;

private:
    BinaryILPSolution solve_impl(const SchedulingProblem &schedProblem, const BinaryILPSolution &bestSolution) const;
};


#endif //ORIGINILPSOLVER_H
