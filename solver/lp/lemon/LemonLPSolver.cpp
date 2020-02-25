#include <lemon/network_simplex.h>
#include "LemonLPSolver.h"

std::unique_ptr<lemon_graph_t> LemonLPSolver::buildGraph(const SchedulingProblem &schedProblem) const {
    std::unique_ptr<lemon_graph_t> graph = std::make_unique<lemon_graph_t>();
    std::vector<int64_t> timePointsOrdered = schedProblem.getOrderedTimePoints();

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        graph->mapJobIndexToNode[i] = graph->networkGraph.addNode();
    }

    for (long timePoint : timePointsOrdered) {
        graph->mapTimePointToNode[timePoint] = graph->networkGraph.addNode();
    }

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        graph->arcs.push_back(graph->networkGraph.addArc(graph->mapJobIndexToNode[i],
                                                         graph->mapTimePointToNode[schedProblem[i].dueDate]));
        graph->arcs.push_back(graph->networkGraph.addArc(graph->mapJobIndexToNode[i],
                                                         graph->mapTimePointToNode[schedProblem[i].deadline]));
    }

    for (size_t i = 0; i < (timePointsOrdered.size() - 1); ++i) {
        graph->arcs.push_back(graph->networkGraph.addArc(graph->mapTimePointToNode[timePointsOrdered[i]],
                                                         graph->mapTimePointToNode[timePointsOrdered[i + 1]]));
    }

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        graph->mapSupplyToNode[graph->mapJobIndexToNode[i]] = schedProblem[i].processingTime;
    }

    graph->mapSupplyToNode[graph->mapTimePointToNode[timePointsOrdered[timePointsOrdered.size() - 1]]] = ((-1) *
                                                                                                          schedProblem.getTotalProcessingTime());

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        graph->mapCapacityToArc[graph->arcs[2 * i]] = schedProblem[i].processingTime;
        graph->mapCapacityToArc[graph->arcs[2 * i + 1]] = schedProblem[i].processingTime;
    }

    for (size_t i = 0; i < (timePointsOrdered.size() - 1); ++i) {
        graph->mapCapacityToArc[graph->arcs[i + 2 * schedProblem.size()]] = timePointsOrdered[i];
    }

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        graph->mapCostToArc[graph->arcs[i * 2]] = -((double) schedProblem[i].weight / (double) schedProblem[i].processingTime);
        graph->mapCostToArc[graph->arcs[i * 2 + 1]] = 0;
    }

    for (size_t i = 0; i < (timePointsOrdered.size() - 1); ++i) {
        graph->mapCostToArc[graph->arcs[i + 2 * schedProblem.size()]] = 0;
    }

    for (size_t i = 0; i < (timePointsOrdered.size() - 1); ++i) {
        graph->mapSupplyToNode[graph->mapTimePointToNode[timePointsOrdered[i]]] = 0;
    }

    return graph;
}

LPSolutionWithBranchPenalty
LemonLPSolver::solve_impl(const SchedulingProblem &schedProblem, bool computeBranchPenalties) const {
    if (schedProblem.empty()) {
        std::cerr << "[LemonLPSolver::solve_impl]: schedProblem is empty" << std::endl;

        return LPSolutionWithBranchPenalty();
    }

    std::unique_ptr<lemon_graph_t> graph = this->buildGraph(schedProblem);

    lemon::NetworkSimplex<lemon::SmartDigraph, int64_t, double> networkSimplex(graph->networkGraph);

    networkSimplex.upperMap(graph->mapCapacityToArc);
    networkSimplex.costMap(graph->mapCostToArc);
    networkSimplex.supplyMap(graph->mapSupplyToNode);

    lemon::NetworkSimplex<lemon::SmartDigraph, int64_t, double>::ProblemType solverResult = networkSimplex.run();
    if (solverResult != lemon::NetworkSimplex<lemon::SmartDigraph, int64_t, double>::OPTIMAL) {
        std::cerr << "[LemonLPSolver::solve_impl]: Solver not find solution" << std::endl;
        return LPSolutionWithBranchPenalty();
    }

    networkSimplex.flowMap(graph->mapFlowToArc);
    networkSimplex.potentialMap(graph->mapPotentialToNode);

    double solutionObj = (-1) * networkSimplex.totalCost();
    std::vector<double> solutionVars;

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        double jobDueDateFlow = graph->mapFlowToArc[graph->arcs[i * 2]];
        solutionVars.push_back((double) (jobDueDateFlow / (double) schedProblem[i].processingTime));
    }

    std::vector<double> solutionReducedCost = this->computeReducedCosts(schedProblem, *graph);
    std::vector<std::pair<int, branch_penalty_t>> branchPenalties;

    if(computeBranchPenalties) {
        branchPenalties = this->computeBranchPenalties(schedProblem,LPSolution(solutionObj, solutionVars, solutionReducedCost));
    }

    std::cerr << "[LemonLPSolver::solve_impl]: UB: " << solutionObj << std::endl;

    return LPSolutionWithBranchPenalty(solutionObj, solutionVars, solutionReducedCost, branchPenalties);
}

std::vector<double>
LemonLPSolver::computeReducedCosts(const SchedulingProblem &schedProblem, lemon_graph_t &graph) const {
    std::vector<double> reducedCosts;

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        double jobDueDateFlow = graph.mapFlowToArc[graph.arcs[i * 2]];
        double jobDueDateCost = graph.mapCostToArc[graph.arcs[i * 2]];
        double jobDueDatePotential = graph.mapPotentialToNode[graph.mapTimePointToNode[schedProblem[i].dueDate]];

        double jobDeadlineFlow = graph.mapFlowToArc[graph.arcs[i * 2 + 1]];
        double jobDeadlineCost = graph.mapCostToArc[graph.arcs[i * 2 + 1]];
        double jobDeadlinePotential = graph.mapPotentialToNode[graph.mapTimePointToNode[schedProblem[i].deadline]];

        double jobNodePotential = graph.mapPotentialToNode[graph.mapJobIndexToNode[i]];

        double networkFlowReducedCostDueDateArc = (jobDueDateCost + jobNodePotential - jobDueDatePotential);
        double networkFlowReducedCostDeadlineArc = (jobDeadlineCost + jobNodePotential - jobDeadlinePotential);

        double originalProblemReducedCost =
                (-1) * (networkFlowReducedCostDueDateArc - networkFlowReducedCostDeadlineArc) *
                std::abs((jobDueDateFlow - jobDeadlineFlow));

        reducedCosts.push_back(originalProblemReducedCost);
    }

    return reducedCosts;
}

std::vector<std::pair<int, branch_penalty_t>>
LemonLPSolver::computeBranchPenalties(const SchedulingProblem &schedProblem, const LPSolution &lpSolution) const {
    std::vector<std::pair<int, branch_penalty_t>> branchingPenalties;

    for (size_t i = 0; i < lpSolution.size(); ++i) {
        if (0 < lpSolution[i] && lpSolution[i] < 1) {
            std::unique_ptr<lemon_graph_t> graph = this->buildGraph(schedProblem);

            int64_t duedateOldCapacity = graph->mapCapacityToArc[graph->arcs[2 * i]];
            int64_t deadlineOldCapacity = graph->mapCapacityToArc[graph->arcs[2 * i + 1]];

            graph->mapCapacityToArc[graph->arcs[2 * i]] = duedateOldCapacity;
            graph->mapCapacityToArc[graph->arcs[2 * i + 1]] = 0;

            lemon::NetworkSimplex<lemon::SmartDigraph, int64_t, double> networkSimplex(graph->networkGraph);
            double upperChange = std::numeric_limits<double>::lowest();
            double lowerChange = std::numeric_limits<double>::lowest();

            networkSimplex.upperMap(graph->mapCapacityToArc);
            networkSimplex.costMap(graph->mapCostToArc);
            networkSimplex.supplyMap(graph->mapSupplyToNode);

            lemon::NetworkSimplex<lemon::SmartDigraph, int64_t, double>::ProblemType solverResult = networkSimplex.run();
            if (solverResult == lemon::NetworkSimplex<lemon::SmartDigraph, int64_t, double>::OPTIMAL) {
                double upperChangeObj = (-1) * networkSimplex.totalCost();
                upperChange = upperChangeObj - lpSolution.getObjective();
            } else {
                std::cerr
                        << "[LemonLPSolver::computeBranchPenalties]: Solver not find solution for branching up on variable "
                        << i << std::endl;
            }

            graph->mapCapacityToArc[graph->arcs[2 * i]] = 0;
            graph->mapCapacityToArc[graph->arcs[2 * i + 1]] = deadlineOldCapacity;

            networkSimplex.upperMap(graph->mapCapacityToArc);
            networkSimplex.costMap(graph->mapCostToArc);
            networkSimplex.supplyMap(graph->mapSupplyToNode);

            solverResult = networkSimplex.run();
            if (solverResult == lemon::NetworkSimplex<lemon::SmartDigraph, int64_t, double>::OPTIMAL) {
                double lowerChangeObj = (-1) * networkSimplex.totalCost();
                lowerChange = lowerChangeObj - lpSolution.getObjective();
            } else {
                std::cerr
                        << "[LemonLPSolver::computeBranchPenalties]: Solver not find solution for branching down on variable "
                        << i << std::endl;
            }

            graph->mapCapacityToArc[graph->arcs[2 * i]] = duedateOldCapacity;
            graph->mapCapacityToArc[graph->arcs[2 * i + 1]] = deadlineOldCapacity;

            branchingPenalties.push_back(std::make_pair<int, branch_penalty_t>(i, {upperChange, lowerChange}));
        }
    }

    return branchingPenalties;
}

LPSolution LemonLPSolver::solve(const SchedulingProblem &schedProblem) const {
    LPSolutionWithBranchPenalty solutionWithBP = this->solve_impl(schedProblem, false);
    return LPSolution(solutionWithBP.getObjective(), solutionWithBP.getVariables(), solutionWithBP.getReducedCosts(),
                      solutionWithBP.getShadowPrices());
}

LPSolutionWithBranchPenalty LemonLPSolver::solveWithBranchPenalty(const SchedulingProblem &schedProblem) const {
    return this->solve_impl(schedProblem, true);
}
