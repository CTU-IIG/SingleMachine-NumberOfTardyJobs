#include "ImprovedUpperBound.h"
#include "../../solver/lp/mcf_class/MCFClassLPSolver.h"
#include "../../reduction/JobReduce.h"

LPSolution ImprovedUpperBound::solve(const SchedulingProblem &schedProblem) const {
    MCFClassLPSolver lpSolver(MCFClassSolver::RELAX_IV);
    LPSolution lpSolution = lpSolver.solve(schedProblem);

    double newObjective = this->recursiveSearch(schedProblem, lpSolution, this->startDepth);
    std::cerr << "[ImprovedUpperBound::solve]: UB: " << newObjective << std::endl;
    return LPSolution(newObjective, lpSolution.getVariables(), lpSolution.getReducedCosts(), lpSolution.getShadowPrices());
}

LPSolutionWithBranchPenalty ImprovedUpperBound::solveWithBranchPenalty(const SchedulingProblem &schedProblem) const {
    MCFClassLPSolver lpSolver(MCFClassSolver::RELAX_IV);
    LPSolutionWithBranchPenalty lpSolution = lpSolver.solveWithBranchPenalty(schedProblem);

    double newObjective = this->recursiveSearch(schedProblem, lpSolution, this->startDepth);
    std::cerr << "[ImprovedUpperBound::solveWithBranchPenalty]: UB: " << newObjective << std::endl;
    return LPSolutionWithBranchPenalty(newObjective, lpSolution.getVariables(), lpSolution.getReducedCosts(), lpSolution.getShadowPrices(), lpSolution.getBranchPenalties());
}

double ImprovedUpperBound::recursiveSearch(const SchedulingProblem &schedProblem, const LPSolution &lpSolution, int depth, int64_t objectiveSum) const {
    if(depth <= 0) {
        return lpSolution.getObjective() + objectiveSum;
    }

    if(lpSolution.getNonIntegerValueCount() == 0) {
        return lpSolution.getObjective() + objectiveSum;
    }

    int branchIndex = this->selectJobWithMostProcessingTimeJob(schedProblem, lpSolution);

    if(branchIndex == -1) {
        return lpSolution.getObjective() + objectiveSum;
    }

    JobReduce jobReduceDownID1(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->startDepth);
    JobReduce jobReduceUpID2(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->startDepth);

    ReducedSchedulingProblem schedulingProblemBranchDownW = jobReduceDownID1.reduceJob(schedProblem, branchIndex, false);
    ReducedSchedulingProblem schedulingProblemBranchUpW = jobReduceUpID2.reduceJob(schedProblem, branchIndex, true);
    SchedulingProblem schedulingProblemBranchDown = schedulingProblemBranchDownW.getReducedSchedulingProblemCopy();
    SchedulingProblem schedulingProblemBranchUp = schedulingProblemBranchUpW.getReducedSchedulingProblemCopy();

    MCFClassLPSolver lpSolverDown(MCFClassSolver::RELAX_IV);
    MCFClassLPSolver lpSolverUp(MCFClassSolver::RELAX_IV);

    LPSolution lpSolutionDown = lpSolverDown.solve(schedulingProblemBranchDown);
    LPSolution lpSolutionUp = lpSolverUp.solve(schedulingProblemBranchUp);

    return std::max(this->recursiveSearch(schedulingProblemBranchDown, lpSolutionDown, depth - 1, objectiveSum), recursiveSearch(schedulingProblemBranchUp, lpSolutionUp, depth - 1, objectiveSum + schedProblem[branchIndex].weight));
}

int ImprovedUpperBound::selectJobWithMostProcessingTimeJob(const SchedulingProblem &schedProblem, const LPSolution &lpSolution) const {
    std::vector<int> nonIntegerIndexes;

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        if((0 < lpSolution[i]) && (lpSolution[i] < 1)) {
            nonIntegerIndexes.push_back(i);
        }
    }

    if(nonIntegerIndexes.empty()) {
        std::cerr << "[ImprovedUpperBound::selectJobWithMostProcessingTimeJob]: Integer solution" << std::endl;
        return -1;
    }

    int mostProcTimeIndex = nonIntegerIndexes[0];
    int mostProcTime = schedProblem[nonIntegerIndexes[0]].processingTime;

    for (size_t i = 0; i < nonIntegerIndexes.size(); ++i) {
        if(schedProblem[i].processingTime > mostProcTime) {
            mostProcTimeIndex = nonIntegerIndexes[i];
            mostProcTime = schedProblem[i].processingTime;
        }
    }

    return mostProcTimeIndex;
}

ImprovedUpperBound::ImprovedUpperBound(ILPSolverImpl ilpSolverImpl,
                                       LowerBoundImpl lowerBoundImpl, UpperBoundImpl upperBoundImpl, int startDepth) : ilpSolverImpl(ilpSolverImpl),
                                                                                        lowerBoundImpl(lowerBoundImpl),
                                                                                        upperBoundImpl(upperBoundImpl),
                                                                                        startDepth(startDepth) {}
