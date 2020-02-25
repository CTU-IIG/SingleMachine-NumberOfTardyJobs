#include <iostream>
#include <algorithm>
#include <numeric>
#include "OriginalLowerBound.h"
#include "../../reduction/JobReduce.h"
#include "../../solver/ilp/OriginILPSolver.h"
#include "../../solver/ilp/ILPGridSolver.h"

std::vector<int>
OriginalLowerBound::computeCoreProblem(const SchedulingProblem &schedProblem, const LPSolution &lpSolution) const {
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
            bool isDominatedByJobJ = false;

            for (size_t i = 0; i < lpSolution.size(); ++i) {
                if (i != j && lpSolution[i] == 1 &&
                    schedProblem[i].processingTime <= schedProblem[j].processingTime &&
                    schedProblem[i].dueDate >= schedProblem[j].dueDate &&
                    schedProblem[i].deadline <= schedProblem[j].deadline &&
                    schedProblem[i].weight >= schedProblem[j].weight) {
                    if (schedProblem[i].processingTime != schedProblem[j].processingTime ||
                        schedProblem[i].dueDate != schedProblem[j].dueDate ||
                        schedProblem[i].deadline != schedProblem[j].deadline ||
                        schedProblem[i].weight != schedProblem[j].weight) {
                        isDominatedByJobJ = true;
                        break;
                    }
                }
            }

            if (!isDominatedByJobJ) {
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
OriginalLowerBound::localSearch(const SchedulingProblem &schedProblem, const BinaryILPSolution &solution) const {
    BinaryILPSolution bestSolution = solution;
    int64_t solutionObj = solution.getObjective();
    std::vector<bool> solutionVar = solution.getVariables();

    std::vector<bool> jobsForcedByDominance;
    for (size_t j = 0; j < schedProblem.size(); ++j) {
        if (solutionVar[j]) {
            bool jobJForceByDominance = false;

            for (size_t i = 0; i < schedProblem.size(); ++i) {
                if (i != j && solutionVar[i] == 1 &&
                    schedProblem[j].processingTime <= schedProblem[i].processingTime &&
                    schedProblem[j].dueDate >= schedProblem[i].dueDate &&
                    schedProblem[j].deadline <= schedProblem[i].deadline &&
                    schedProblem[j].weight >= schedProblem[i].weight) {
                    if (schedProblem[j].processingTime != schedProblem[i].processingTime ||
                        schedProblem[j].dueDate != schedProblem[i].dueDate ||
                        schedProblem[j].deadline != schedProblem[i].deadline ||
                        schedProblem[j].weight != schedProblem[i].weight) {
                        jobJForceByDominance = true;
                        break;
                    }
                }
            }

            jobsForcedByDominance.push_back(jobJForceByDominance);
        } else {
            bool jobJForceByDominance = false;

            for (size_t i = 0; i < schedProblem.size(); ++i) {
                if (i != j && solutionVar[i] == 0 &&
                    schedProblem[i].processingTime <= schedProblem[j].processingTime &&
                    schedProblem[i].dueDate >= schedProblem[j].dueDate &&
                    schedProblem[i].deadline <= schedProblem[j].deadline &&
                    schedProblem[i].weight >= schedProblem[j].weight) {
                    if (schedProblem[i].processingTime != schedProblem[j].processingTime ||
                        schedProblem[i].dueDate != schedProblem[j].dueDate ||
                        schedProblem[i].deadline != schedProblem[j].deadline ||
                        schedProblem[i].weight != schedProblem[j].weight) {
                        jobJForceByDominance = true;
                        break;
                    }
                }
            }

            jobsForcedByDominance.push_back(jobJForceByDominance);
        }
    }

    std::vector<int> preSortedSchedule = schedProblem.createScheduleFromSolution(solutionVar);

    for (size_t j = 0; j < schedProblem.size(); ++j) {
        if (jobsForcedByDominance[j]) {
            continue;
        }

        for (size_t i = j + 1; i < schedProblem.size(); ++i) {
            if (jobsForcedByDominance[i]) {
                continue;
            }

            if (solutionVar[j] == solutionVar[i]) {
                continue;
            }

            if (solutionVar[i] == false && solutionVar[j] == true && schedProblem[i].weight <= schedProblem[j].weight) {
                continue;
            }

            if (solutionVar[i] == true && solutionVar[j] == false && schedProblem[i].weight >= schedProblem[j].weight) {
                continue;
            }

            bool oldStateJobI = solutionVar[i];
            bool oldStateJobJ = solutionVar[j];

            solutionVar[i] = oldStateJobJ;
            solutionVar[j] = oldStateJobI;

            if (schedProblem.isFeasibleSolution(solutionVar, preSortedSchedule)) {
                if (solution[i] == false && solution[j] == true &&
                    (solutionObj + schedProblem[i].weight - schedProblem[j].weight) > bestSolution.getObjective()) {
                    int64_t newSolutionObjective = (solutionObj + schedProblem[i].weight - schedProblem[j].weight);
                    bestSolution = BinaryILPSolution(newSolutionObjective, solutionVar);

                    std::cout << "[OriginalLowerBound::localSearch]: Solution improved from " << solution.getObjective()
                              << ", to " << newSolutionObjective << std::endl;
                }

                if (solution[i] == true && solution[j] == false &&
                    (solutionObj - schedProblem[i].weight + schedProblem[j].weight) > bestSolution.getObjective()) {
                    int64_t newSolutionObjective = (solutionObj - schedProblem[i].weight + schedProblem[j].weight);
                    bestSolution = BinaryILPSolution(newSolutionObjective, solutionVar);

                    std::cout << "[OriginalLowerBound::localSearch]: Solution improved from " << solution.getObjective()
                              << ", to " << newSolutionObjective << std::endl;
                }
            }

            solutionVar[i] = oldStateJobI;
            solutionVar[j] = oldStateJobJ;
        }
    }

    return bestSolution;
}

BinaryILPSolution OriginalLowerBound::solve(const SchedulingProblem &schedProblem, const LPSolution &lpSolution) const {
    JobReduce jobReduce(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);
    std::vector<int> coreSet = this->computeCoreProblem(schedProblem, lpSolution);
    std::vector<std::pair<int, bool>> reduceJobs = this->getReduceJobsFromCoreSet(schedProblem, lpSolution, coreSet);
    ReducedSchedulingProblem reducedSchedProblem = jobReduce.reduceJobs(schedProblem, reduceJobs, true);
    reduceJobs = reducedSchedProblem.getReduceJobs();

    SchedulingProblem schedProblemReduced = reducedSchedProblem.getReducedSchedulingProblemCopy();

    std::cerr << "[OriginalLowerBound::solve]: LowerBound reducedJobSet size for ILP: " << reducedSchedProblem.size() << std::endl;

    BinaryILPSolution ilpSolution;

    if(this->ilpSolverImpl == ILPSolverImpl::ORIGINAL) {
        OriginILPSolver ilpSolver(1, false, false);
        ilpSolution = ilpSolver.solve(schedProblemReduced);
    } else if(this->ilpSolverImpl == ILPSolverImpl::IMPROVED) {
        ILPGridSolver ilpSolver(1, false, false);
        ilpSolution = ilpSolver.solve(schedProblemReduced);
    }

    BinaryILPSolution lowerBoundSolution = reducedSchedProblem.reconstructSchedulingProblem(ilpSolution);

    std::cerr << "[OriginalLowerBound::solve]: LB: " << lowerBoundSolution.getObjective() << std::endl;

    // LocalSearch
    BinaryILPSolution lowerBoundSolutionLS = this->localSearch(schedProblem, lowerBoundSolution);

    if (lowerBoundSolutionLS.getObjective() > lowerBoundSolution.getObjective()) {
        lowerBoundSolution = lowerBoundSolutionLS;
    }

    return lowerBoundSolution;
}

std::vector<std::pair<int, bool>>
OriginalLowerBound::getReduceJobsFromCoreSet(const SchedulingProblem &schedProblem, const LPSolution &lpSolution,
                                             const std::vector<int> &coreSet) const {
    if (lpSolution.size() != schedProblem.size()) {
        std::cerr << "[OriginalLowerBound::getReduceJobsFromCoreSet]: schedProblem.size() != lpSolution.size()"
                  << std::endl;
        return std::vector<std::pair<int, bool>>();
    }

    std::vector<std::pair<int, bool>> reduceJobs;
    std::vector<int> coreSetCopy = coreSet;

    std::sort(coreSetCopy.begin(), coreSetCopy.end());

    int lastInsertedIndex = -1;
    for (const int &coreJob : coreSetCopy) {
        for (int i = (lastInsertedIndex + 1); i < coreJob; ++i) {
            reduceJobs.emplace_back(i, lpSolution[i]);
        }

        lastInsertedIndex = coreJob;
    }

    for (size_t i = (lastInsertedIndex + 1); i < schedProblem.size(); ++i) {
        reduceJobs.emplace_back(i, lpSolution[i]);
    }

    return reduceJobs;
}

OriginalLowerBound::OriginalLowerBound(ILPSolverImpl ilpSolverImpl, LowerBoundImpl lowerBoundImpl, UpperBoundImpl upperBoundImpl, int upperBoundMaxDepth) : ilpSolverImpl(ilpSolverImpl), lowerBoundImpl(lowerBoundImpl), upperBoundImpl(upperBoundImpl), upperBoundMaxDepth(upperBoundMaxDepth) {}
