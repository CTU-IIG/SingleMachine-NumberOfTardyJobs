#ifndef LEMONLPSOLVER_H
#define LEMONLPSOLVER_H


#include <lemon/smart_graph.h>
#include <unordered_map>
#include "../LPSolver.h"
#include "../LPSolutionWithBranchPenalty.h"

typedef struct lemon_graph_struct {
    lemon::SmartDigraph networkGraph;
    std::vector<lemon::SmartDigraph::Arc> arcs;

    std::unordered_map<int32_t, lemon::SmartDigraph::Node> mapJobIndexToNode;
    std::unordered_map<int32_t, lemon::SmartDigraph::Node> mapTimePointToNode;

    lemon::SmartDigraph::ArcMap<int64_t> mapCapacityToArc;
    lemon::SmartDigraph::ArcMap<double> mapCostToArc;
    lemon::SmartDigraph::NodeMap<int64_t> mapSupplyToNode;

    lemon::SmartDigraph::ArcMap<int64_t> mapFlowToArc;
    lemon::SmartDigraph::NodeMap<double> mapPotentialToNode;

    lemon_graph_struct() : mapCapacityToArc(networkGraph), mapCostToArc(networkGraph), mapSupplyToNode(networkGraph),
                           mapFlowToArc(networkGraph), mapPotentialToNode(networkGraph) {}
} lemon_graph_t;

class LemonLPSolver : public LPSolver {
public:
    LemonLPSolver() = default;

    virtual ~LemonLPSolver() = default;

    LPSolution solve(const SchedulingProblem &schedProblem) const override;

    LPSolutionWithBranchPenalty solveWithBranchPenalty(const SchedulingProblem &schedProblem) const override;

private:
    LPSolutionWithBranchPenalty solve_impl(const SchedulingProblem &schedProblem, bool computeBranchPenalties = false) const;

    std::unique_ptr<lemon_graph_t> buildGraph(const SchedulingProblem &schedProblem) const;

    std::vector<double> computeReducedCosts(const SchedulingProblem &schedProblem, lemon_graph_t &graph) const;

    std::vector<std::pair<int, branch_penalty_t>>
    computeBranchPenalties(const SchedulingProblem &schedProblem, const LPSolution &lpSolution) const;
};


#endif //LEMONLPSOLVER_H
