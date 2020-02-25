#include "ImprovedLowerBound.h"

std::vector<int>
ImprovedLowerBound::computeCoreProblem(const SchedulingProblem &schedProblem, const LPSolution &lpSolution) const {
    std::vector<int> coreSet;

    for (size_t j = 0; j < lpSolution.size(); ++j) {
        bool isJobInCoreSet = false;

        if (0 < lpSolution[j] && lpSolution[j] < 1) {
            coreSet.push_back(j);
            isJobInCoreSet = true;
        } else if (lpSolution[j] == 0) {
            bool isDominatingJobJ = false;

            for (size_t i = 0; i < lpSolution.size(); ++i) {
                if (i != j && lpSolution[i] == 0 &&
                    schedProblem[i].processingTime <= schedProblem[j].processingTime &&
                    schedProblem[i].dueDate >= schedProblem[j].dueDate &&
                    schedProblem[i].deadline <= schedProblem[j].deadline &&
                    schedProblem[i].weight >= schedProblem[j].weight) {
                    if (schedProblem[i].processingTime != schedProblem[j].processingTime ||
                        schedProblem[i].dueDate != schedProblem[j].dueDate ||
                        schedProblem[i].deadline != schedProblem[j].deadline ||
                        schedProblem[i].weight != schedProblem[j].weight) {
                        isDominatingJobJ = true;
                        break;
                    }
                }
            }

            if (!isDominatingJobJ) {
                coreSet.push_back(j);
                isJobInCoreSet = true;
            }
        } else if (lpSolution[j] == 1) {
            bool isDominatingJobJ = false;

            for (size_t i = 0; i < lpSolution.size(); ++i) {
                if (i != j && lpSolution[i] == 1 &&
                    schedProblem[j].processingTime <= schedProblem[i].processingTime &&
                    schedProblem[j].dueDate >= schedProblem[i].dueDate &&
                    schedProblem[j].deadline <= schedProblem[i].deadline &&
                    schedProblem[j].weight >= schedProblem[i].weight) {
                    if (schedProblem[j].processingTime != schedProblem[i].processingTime ||
                        schedProblem[j].dueDate != schedProblem[i].dueDate ||
                        schedProblem[j].deadline != schedProblem[i].deadline ||
                        schedProblem[j].weight != schedProblem[i].weight) {
                        isDominatingJobJ = true;
                        break;
                    }
                }
            }

            if (!isDominatingJobJ) {
                coreSet.push_back(j);
                isJobInCoreSet = true;
            }
        }

        if (!isJobInCoreSet) {
            //this->jobsToReduce.emplace_back(j, (lpSolution[j] == 1));
        }
    }

    return coreSet;
}

BinaryILPSolution
ImprovedLowerBound::localSearch(const SchedulingProblem &schedProblem, const BinaryILPSolution &solution) const {
    return solution;
}

ImprovedLowerBound::ImprovedLowerBound(ILPSolverImpl ilpSolverImpl, LowerBoundImpl lowerBoundImpl, UpperBoundImpl upperBoundImpl, int upperBoundMaxDepth) : OriginalLowerBound(
        ilpSolverImpl, lowerBoundImpl, upperBoundImpl, upperBoundMaxDepth) {}
