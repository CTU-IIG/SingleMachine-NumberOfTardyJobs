#ifndef ILPGRIDSOLVER_H
#define ILPGRIDSOLVER_H


#include "ILPSolver.h"
#include "../lp/LPSolutionWithEarlyRange.h"

enum class EarlyJobConstrainType {NONE, EQUAL, LOWER_EQUAL, GREATER_EQUAL};

class ILPGridSolver : public ILPSolver {
private:
    uint32_t threadCount = 1;
    bool printOutput = false;
    bool usePresolve = false;
public:
    explicit ILPGridSolver(uint32_t threadCount, bool printOutput, bool usePresolve);

    ILPGridSolver() = delete;

    virtual ~ILPGridSolver() = default;

    BinaryILPSolution solve(const SchedulingProblem &schedProblem) const override;

    BinaryILPSolution solve(const SchedulingProblem &schedProblem, const BinaryILPSolution &lowerBound) const override;

private:
    LPSolutionWithEarlyRange getEarlyRange(const SchedulingProblem &schedProblem) const;

    BinaryILPSolution solve_impl(const SchedulingProblem &schedProblem, BinaryILPSolution bestSolution) const;

    BinaryILPSolution solveWithEarlyConstrain(const SchedulingProblem &schedProblem, int earlyJobCount, EarlyJobConstrainType constrainType, int64_t currentBestObjextive) const;
};


#endif //ILPGRIDSOLVER_H
