#include <iostream>
#include <RelaxIV.h>
#include "ImprovedUpperBoundOpt.h"

double ImprovedUpperBoundOpt::recursiveSearch(MCFClass *mcfSolver, const SchedulingProblem &schedProblem, const LPSolution &lpSolution, int depth) const {
    if(depth <= 0) {
        return lpSolution.getObjective();
    }

    if(lpSolution.getNonIntegerValueCount() == 0) {
        return lpSolution.getObjective();
    }

    int branchIndex = -1;

    if(schedProblem.hasDeadlines()) {
        branchIndex = this->selectJobByDeadline(schedProblem, lpSolution, false);
    } else {
        branchIndex = this->selectJobByDueDate(schedProblem, lpSolution, true);
    }

//    branchIndex = this->selectJobByDeadline(schedProblem, lpSolution, false);
//    branchIndex = this->selectJobByDeadline(schedProblem, lpSolution, true);
//    branchIndex = this->selectJobByProcessingTime(schedProblem, lpSolution, false);
//    branchIndex = this->selectJobByProcessingTime(schedProblem, lpSolution, true);
//    branchIndex = this->selectJobByDueDate(schedProblem, lpSolution, false);
//    branchIndex = this->selectJobByDueDate(schedProblem, lpSolution, true);
//    branchIndex = this->selectJobByMaxMinPseudoCost(this->computeBranchPenalties(schedProblem, lpSolution, mcfSolver));

    if(branchIndex == -1) {
        return lpSolution.getObjective();
    }

    double recursUpperObj = std::numeric_limits<double>::lowest();
    double recursLowerObj = std::numeric_limits<double>::lowest();

    int64_t duedateOldCapacity = mcfSolver->MCFUCap(2 * branchIndex);
    int64_t deadlineOldCapacity =  mcfSolver->MCFUCap(2 * branchIndex + 1);

    mcfSolver->ChgUCap(2 * branchIndex, duedateOldCapacity);
    mcfSolver->ChgUCap(2 * branchIndex + 1, 0);


    mcfSolver->SolveMCF();
    if (mcfSolver->MCFGetStatus() == MCFClass::kOK) {
        double solutionObj = (-1) * mcfSolver->MCFGetFO();
        std::vector<double> solutionVars;

        MCFClass::FRow arcFlow = new MCFClass::FNumber[mcfSolver->MCFm()];
        mcfSolver->MCFGetX(arcFlow);

        for (size_t i = 0; i < schedProblem.size(); ++i) {
            double jobDueDateFlow = arcFlow[i * 2];
            solutionVars.push_back((double) (jobDueDateFlow / (double) schedProblem[i].processingTime));
        }

        delete[] arcFlow;

        recursUpperObj = this->recursiveSearch(mcfSolver, schedProblem, LPSolution(solutionObj, solutionVars), depth - 1);
    } else {
        std::cerr << "[ImprovedUpperBoundOpt::recursiveSearch]: Solver not find solution for branching up on variable " << branchIndex << std::endl;
    }

    mcfSolver->ChgUCap(2 * branchIndex, 0);
    mcfSolver->ChgUCap(2 * branchIndex + 1, deadlineOldCapacity);

    mcfSolver->SolveMCF();
    if (mcfSolver->MCFGetStatus() == MCFClass::kOK) {
        double solutionObj = (-1) * mcfSolver->MCFGetFO();
        std::vector<double> solutionVars;

        MCFClass::FRow arcFlow = new MCFClass::FNumber[mcfSolver->MCFm()];
        mcfSolver->MCFGetX(arcFlow);

        for (size_t i = 0; i < schedProblem.size(); ++i) {
            double jobDueDateFlow = arcFlow[i * 2];
            solutionVars.push_back((double) (jobDueDateFlow / (double) schedProblem[i].processingTime));
        }

        delete[] arcFlow;

        recursLowerObj = this->recursiveSearch(mcfSolver, schedProblem, LPSolution(solutionObj, solutionVars), depth - 1);
    } else {
        std::cerr << "[ImprovedUpperBoundOpt::recursiveSearch]: Solver not find solution for branching down on variable " << branchIndex << std::endl;
    }

    mcfSolver->ChgUCap(2 * branchIndex, duedateOldCapacity);
    mcfSolver->ChgUCap(2 * branchIndex + 1, deadlineOldCapacity);

    return std::max(recursLowerObj, recursUpperObj);
}

LPSolution ImprovedUpperBoundOpt::solve(const SchedulingProblem &schedProblem) const {
    if (schedProblem.empty()) {
        std::cerr << "[ImprovedUpperBoundOpt::solve]: schedProblem is empty" << std::endl;

        return LPSolution();
    }

    std::unique_ptr<mfc_class_graph_t> graph = MCFClassLPSolver::buildGraph(schedProblem);
    std::unique_ptr<MCFClass> mcfSolver = std::make_unique<RelaxIV>();

    mcfSolver->LoadNet(graph->nodeCount, graph->arcCount, graph->nodeCount, graph->arcCount, graph->arcCapacity, graph->arcCost, graph->nodeDeficit, graph->arcStartNodeId, graph->arcEndNodeId);
    mcfSolver->SolveMCF();

    if (mcfSolver->MCFGetStatus() != MCFClass::kOK) {
        std::cerr << "[ImprovedUpperBoundOpt::solve]: Solver not find solution. Status: "<< mcfSolver->MCFGetStatus() << std::endl;
        return LPSolution();
    }

    mcfSolver->MCFGetX(graph->arcFlow);
    mcfSolver->MCFGetRC(graph->arcReducedCost);
    mcfSolver->MCFGetPi(graph->nodePotential);

    double solutionObj = (-1) * mcfSolver->MCFGetFO();
    std::vector<double> solutionVars;

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        double jobDueDateFlow = graph->arcFlow[i * 2];
        solutionVars.push_back((double) (jobDueDateFlow / (double) schedProblem[i].processingTime));
    }

    std::vector<double> solutionReducedCost = this->computeReducedCosts(schedProblem, *graph);
    LPSolution lpSolution(solutionObj, solutionVars, solutionReducedCost);

    double newObjective = this->recursiveSearch(mcfSolver.get(), schedProblem, lpSolution, this->startDepth);
    std::cerr << "[ImprovedUpperBoundOpt::solve]: UB: " << newObjective << std::endl;
    return LPSolution(newObjective, lpSolution.getVariables(), lpSolution.getReducedCosts(), lpSolution.getShadowPrices());
}

LPSolutionWithBranchPenalty ImprovedUpperBoundOpt::solveWithBranchPenalty(const SchedulingProblem &schedProblem) const {
    MCFClassLPSolver lpSolver(MCFClassSolver::RELAX_IV);
    LPSolutionWithBranchPenalty lpSolutionWithBP = lpSolver.solveWithBranchPenalty(schedProblem);

    LPSolution improvedSolution = this->solve(schedProblem);
    return LPSolutionWithBranchPenalty(improvedSolution.getObjective(), lpSolutionWithBP.getVariables(), lpSolutionWithBP.getReducedCosts(), lpSolutionWithBP.getShadowPrices(), lpSolutionWithBP.getBranchPenalties());
}

ImprovedUpperBoundOpt::ImprovedUpperBoundOpt(ILPSolverImpl ilpSolverImpl, LowerBoundImpl lowerBoundImpl, UpperBoundImpl upperBoundImpl, int startDepth) : ilpSolverImpl(
        ilpSolverImpl), lowerBoundImpl(lowerBoundImpl), upperBoundImpl(upperBoundImpl), startDepth(startDepth) {}

std::vector<double>
ImprovedUpperBoundOpt::computeReducedCosts(const SchedulingProblem &schedProblem, mfc_class_graph_t &graph) const {
    std::vector<double> reducedCosts;

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        double jobDueDateFlow = graph.arcFlow[i * 2];
        double dueDateArcReducedCost = graph.arcReducedCost[i * 2];

        double jobDeadlineFlow = graph.arcFlow[i * 2 + 1];
        double deadlineArcReducedCost = graph.arcReducedCost[i * 2 + 1];

        double originalProblemReducedCost =
                (-1) * (dueDateArcReducedCost - deadlineArcReducedCost) * std::abs((jobDueDateFlow - jobDeadlineFlow));

        reducedCosts.push_back(originalProblemReducedCost);
    }

    return reducedCosts;
}

std::vector<std::pair<int, branch_penalty_t>>
ImprovedUpperBoundOpt::computeBranchPenalties(const SchedulingProblem &schedProblem, const LPSolution &lpSolution, MCFClass *mcfSolver) const {
    std::vector<std::pair<int, branch_penalty_t>> branchingPenalties;

    for (size_t i = 0; i < lpSolution.size(); ++i) {
        if (0 < lpSolution[i] && lpSolution[i] < 1) {
            int64_t duedateOldCapacity = mcfSolver->MCFUCap(2 * i);
            int64_t deadlineOldCapacity =  mcfSolver->MCFUCap(2 * i + 1);

            mcfSolver->ChgUCap(2 * i, duedateOldCapacity);
            mcfSolver->ChgUCap(2 * i + 1, 0);

            double upperChange = std::numeric_limits<double>::lowest();
            double lowerChange = std::numeric_limits<double>::lowest();

            mcfSolver->SolveMCF();
            if (mcfSolver->MCFGetStatus() == MCFClass::kOK) {
                double upperChangeObj = (-1) * mcfSolver->MCFGetFO();
                upperChange = upperChangeObj - lpSolution.getObjective();
            } else {
                std::cerr
                        << "[ImprovedUpperBoundOpt::computeBranchPenalties]: Solver not find solution for branching up on variable "
                        << i << std::endl;
            }

            mcfSolver->ChgUCap(2 * i, 0);
            mcfSolver->ChgUCap(2 * i + 1, deadlineOldCapacity);

            mcfSolver->SolveMCF();
            if (mcfSolver->MCFGetStatus() == MCFClass::kOK) {
                double lowerChangeObj = (-1) * mcfSolver->MCFGetFO();
                lowerChange = lowerChangeObj - lpSolution.getObjective();
            } else {
                std::cerr
                        << "[ImprovedUpperBoundOpt::computeBranchPenalties]: Solver not find solution for branching down on variable "
                        << i << std::endl;
            }

            mcfSolver->ChgUCap(2 * i, duedateOldCapacity);
            mcfSolver->ChgUCap(2 * i + 1, deadlineOldCapacity);

            branchingPenalties.push_back(std::make_pair<int, branch_penalty_t>(i, {upperChange, lowerChange}));
        }
    }

    return branchingPenalties;
}

int ImprovedUpperBoundOpt::selectJobByProcessingTime(const SchedulingProblem &schedProblem, const LPSolution &lpSolution, bool max) const {
    std::vector<int> nonIntegerIndexes;

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        if((0 < lpSolution[i]) && (lpSolution[i] < 1)) {
            nonIntegerIndexes.push_back(i);
        }
    }

    if(nonIntegerIndexes.empty()) {
        std::cerr << "[ImprovedUpperBoundOpt::selectJobByProcessingTime]: Integer solution" << std::endl;
        return -1;
    }

    int bestAtribIndex = nonIntegerIndexes[0];
    int bestAtribValue = schedProblem[nonIntegerIndexes[0]].processingTime;

    for (size_t i = 0; i < nonIntegerIndexes.size(); ++i) {
        int jobAtribVal = schedProblem[i].processingTime;

        if(max && (jobAtribVal > bestAtribValue)) {
            bestAtribIndex = nonIntegerIndexes[i];
            bestAtribValue = jobAtribVal;
        } else if(!max && (jobAtribVal < bestAtribValue)) {
            bestAtribIndex = nonIntegerIndexes[i];
            bestAtribValue = jobAtribVal;
        }
    }

    return bestAtribIndex;
}

int ImprovedUpperBoundOpt::selectJobByDueDate(const SchedulingProblem &schedProblem, const LPSolution &lpSolution, bool max) const {
    std::vector<int> nonIntegerIndexes;

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        if((0 < lpSolution[i]) && (lpSolution[i] < 1)) {
            nonIntegerIndexes.push_back(i);
        }
    }

    if(nonIntegerIndexes.empty()) {
        std::cerr << "[ImprovedUpperBoundOpt::selectJobByDueDate]: Integer solution" << std::endl;
        return -1;
    }

    int bestAtribIndex = nonIntegerIndexes[0];
    int bestAtribValue = schedProblem[nonIntegerIndexes[0]].dueDate;

    for (size_t i = 0; i < nonIntegerIndexes.size(); ++i) {
        int jobAtribVal = schedProblem[i].dueDate;

        if(max && (jobAtribVal > bestAtribValue)) {
            bestAtribIndex = nonIntegerIndexes[i];
            bestAtribValue = jobAtribVal;
        } else if(!max && (jobAtribVal < bestAtribValue)) {
            bestAtribIndex = nonIntegerIndexes[i];
            bestAtribValue = jobAtribVal;
        }
    }

    return bestAtribIndex;
}

int ImprovedUpperBoundOpt::selectJobByDeadline(const SchedulingProblem &schedProblem, const LPSolution &lpSolution, bool max) const {
    std::vector<int> nonIntegerIndexes;

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        if((0 < lpSolution[i]) && (lpSolution[i] < 1)) {
            nonIntegerIndexes.push_back(i);
        }
    }

    if(nonIntegerIndexes.empty()) {
        std::cerr << "[ImprovedUpperBoundOpt::selectJobByDeadline]: Integer solution" << std::endl;
        return -1;
    }

    int bestAtribIndex = nonIntegerIndexes[0];
    int bestAtribValue = schedProblem[nonIntegerIndexes[0]].deadline;

    for (size_t i = 0; i < nonIntegerIndexes.size(); ++i) {
        int jobAtribVal = schedProblem[i].deadline;

        if(max && (jobAtribVal > bestAtribValue)) {
            bestAtribIndex = nonIntegerIndexes[i];
            bestAtribValue = jobAtribVal;
        } else if(!max && (jobAtribVal < bestAtribValue)) {
            bestAtribIndex = nonIntegerIndexes[i];
            bestAtribValue = jobAtribVal;
        }
    }

    return bestAtribIndex;
}

int ImprovedUpperBoundOpt::selectJobByMaxMinPseudoCost(const std::vector<std::pair<int, branch_penalty_t>> &branchingPenalties) const {
    double maxMinPseudoCost = std::numeric_limits<double>::lowest();
    int maxMinPseudoCostIndex = -1;

    for (const std::pair<int, branch_penalty_t> &branchPenaltyPair : branchingPenalties) {
        double minPseudoCost = std::min(std::abs(branchPenaltyPair.second.lowerChange), std::abs(branchPenaltyPair.second.upperChange));

        if (minPseudoCost > maxMinPseudoCost) {
            maxMinPseudoCost = minPseudoCost;
            maxMinPseudoCostIndex = branchPenaltyPair.first;
        }
    }

    return maxMinPseudoCostIndex;
}