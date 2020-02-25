#ifndef ILPSOLVER_H
#define ILPSOLVER_H


#include "BinaryILPSolution.h"
#include "../../scheduling/SchedulingProblem.h"

class ILPSolver {
public:
    ILPSolver() = default;

    virtual ~ILPSolver() = default;

    virtual BinaryILPSolution solve(const SchedulingProblem &schedProblem) const = 0;

    virtual BinaryILPSolution solve(const SchedulingProblem &schedProblem, const BinaryILPSolution &lowerBound) const = 0;
};


#endif //ILPSOLVER_H
