#include <iostream>
#include <cmath>
#include "MCFClassLPSolver.h"
#include "../../../libs/mcf_class/MCFSimplex/MCFSimplex.h"
#include "../../../libs/mcf_class/RelaxIV/RelaxIV.h"
#include "../lemon/LemonLPSolver.h"
#include "../GurobiLPSolver.h"

MCFClassLPSolver::MCFClassLPSolver(const MCFClassSolver solverType) : solverType(solverType) {}

std::unique_ptr<mfc_class_graph_t> MCFClassLPSolver::buildGraph(const SchedulingProblem &schedProblem) {
    std::unique_ptr<mfc_class_graph_t> graph = std::make_unique<mfc_class_graph_t>();
    std::vector<int64_t> timePointsOrdered = schedProblem.getOrderedTimePoints();

    graph->nodeCount = schedProblem.size() + timePointsOrdered.size();
    graph->arcCount = 2 * schedProblem.size() + timePointsOrdered.size() - 1;

    graph->arcCapacity = new MCFClass::FNumber[graph->arcCount];
    graph->arcCost = new MCFClass::CNumber[graph->arcCount];
    graph->nodeDeficit = new MCFClass::FNumber[graph->nodeCount];

    graph->arcStartNodeId = new MCFClass::Index[graph->arcCount];
    graph->arcEndNodeId = new MCFClass::Index[graph->arcCount];

    graph->arcFlow = new MCFClass::FNumber[graph->arcCount];
    graph->arcReducedCost = new MCFClass::CNumber[graph->arcCount];
    graph->nodePotential = new MCFClass::CNumber[graph->nodeCount];

    for (size_t i = 0; i < graph->nodeCount; ++i) {
        graph->nodeDeficit[i] = 0;
        graph->nodePotential[i] = 0;
    }

    for (size_t i = 0; i < graph->arcCount; ++i) {
        graph->arcCapacity[i] = 0;
        graph->arcCost[i] = 0;
        graph->arcStartNodeId[i] = 0;
        graph->arcEndNodeId[i] = 0;
        graph->arcFlow[i] = 0;
        graph->arcReducedCost[i] = 0;
    }

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        graph->mapJobIndexToNode[i] = i + 1;
    }

    for (size_t i = 0; i < timePointsOrdered.size(); ++i) {
        graph->mapTimePointToNode[timePointsOrdered[i]] = schedProblem.size() + i + 1;
    }

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        // Duedate
        graph->arcStartNodeId[2 * i] = graph->mapJobIndexToNode[i];
        graph->arcEndNodeId[2 * i] = graph->mapTimePointToNode[schedProblem[i].dueDate];

        // Deadline
        graph->arcStartNodeId[2 * i + 1] = graph->mapJobIndexToNode[i];
        graph->arcEndNodeId[2 * i + 1] = graph->mapTimePointToNode[schedProblem[i].deadline];
    }

    for (size_t i = 0; i < (timePointsOrdered.size() - 1); ++i) {
        graph->arcStartNodeId[2 * schedProblem.size() + i] = graph->mapTimePointToNode[timePointsOrdered[i]];
        graph->arcEndNodeId[2 * schedProblem.size() + i] = graph->mapTimePointToNode[timePointsOrdered[i + 1]];
    }

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        graph->nodeDeficit[graph->mapJobIndexToNode[i] - 1] = -schedProblem[i].processingTime;
    }

    graph->nodeDeficit[graph->mapTimePointToNode[timePointsOrdered[timePointsOrdered.size() - 1]] -
                       1] = schedProblem.getTotalProcessingTime();

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        graph->arcCapacity[2 * i] = schedProblem[i].processingTime;
        graph->arcCapacity[2 * i + 1] = schedProblem[i].processingTime;
    }

    for (size_t i = 0; i < (timePointsOrdered.size() - 1); ++i) {
        graph->arcCapacity[i + 2 * schedProblem.size()] = timePointsOrdered[i];
    }

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        graph->arcCost[i * 2] = -((double) schedProblem[i].weight / (double) schedProblem[i].processingTime);
    }

    return graph;
}

LPSolutionWithBranchPenalty
MCFClassLPSolver::solve_impl(const SchedulingProblem &schedProblem, bool computeBranchPenalties) const {
    if (schedProblem.empty()) {
        std::cerr << "[MCFClassLPSolver::solve_impl]: schedProblem is empty" << std::endl;

        return LPSolutionWithBranchPenalty();
    }

    std::unique_ptr<mfc_class_graph_t> graph = this->buildGraph(schedProblem);
    std::unique_ptr<MCFClass> mcfSolver = nullptr;

    if (this->solverType == MCFClassSolver::RELAX_IV) {
        mcfSolver = std::make_unique<RelaxIV>();
    } else if (this->solverType == MCFClassSolver::NETWORK_SIMPLEX) {
        mcfSolver = std::make_unique<MCFSimplex>();
    }

    mcfSolver->LoadNet(graph->nodeCount, graph->arcCount, graph->nodeCount, graph->arcCount, graph->arcCapacity,
                       graph->arcCost, graph->nodeDeficit, graph->arcStartNodeId, graph->arcEndNodeId);
    mcfSolver->SolveMCF();

    if (mcfSolver->MCFGetStatus() != MCFClass::kOK) {
        std::cerr << "[MCFClassLPSolver::solve_impl]: Solver not find solution. Status: "<< mcfSolver->MCFGetStatus() << std::endl;
        return LPSolutionWithBranchPenalty();
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
    std::vector<std::pair<int, branch_penalty_t>> branchPenalties;

    if (computeBranchPenalties) {
        branchPenalties = this->computeBranchPenalties(schedProblem,
                                                       LPSolution(solutionObj, solutionVars, solutionReducedCost),
                                                       *graph, mcfSolver.get());
    }

    std::cerr << "[MCFClassLPSolver::solve_impl]: UB: " << solutionObj << std::endl;

    return LPSolutionWithBranchPenalty(solutionObj, solutionVars, solutionReducedCost, branchPenalties);
}

std::vector<double>
MCFClassLPSolver::computeReducedCosts(const SchedulingProblem &schedProblem, mfc_class_graph_t &graph) const {
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
MCFClassLPSolver::computeBranchPenalties(const SchedulingProblem &schedProblem, const LPSolution &lpSolution,
                                         mfc_class_graph_t &graph, MCFClass *mcfSolver) const {
    std::vector<std::pair<int, branch_penalty_t>> branchingPenalties;

    for (size_t i = 0; i < lpSolution.size(); ++i) {
        if (0 < lpSolution[i] && lpSolution[i] < 1) {
            int64_t duedateOldCapacity = graph.arcCapacity[2 * i];
            int64_t deadlineOldCapacity = graph.arcCapacity[2 * i + 1];

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
                        << "[MCFClassLPSolver::computeBranchPenalties]: Solver not find solution for branching up on variable "
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
                        << "[MCFClassLPSolver::computeBranchPenalties]: Solver not find solution for branching down on variable "
                        << i << std::endl;
            }

            mcfSolver->ChgUCap(2 * i, duedateOldCapacity);
            mcfSolver->ChgUCap(2 * i + 1, deadlineOldCapacity);

            branchingPenalties.push_back(std::make_pair<int, branch_penalty_t>(i, {upperChange, lowerChange}));
        }
    }

    return branchingPenalties;
}

LPSolution MCFClassLPSolver::solve(const SchedulingProblem &schedProblem) const {
    LPSolutionWithBranchPenalty solutionWithBP = this->solve_impl(schedProblem, false);
    return LPSolution(solutionWithBP.getObjective(), solutionWithBP.getVariables(), solutionWithBP.getReducedCosts(),
                      solutionWithBP.getShadowPrices());
}

LPSolutionWithBranchPenalty MCFClassLPSolver::solveWithBranchPenalty(const SchedulingProblem &schedProblem) const {
    return this->solve_impl(schedProblem, true);
}
