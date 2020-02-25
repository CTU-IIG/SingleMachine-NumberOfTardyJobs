#ifndef MCFCLASSLPSOLVER_H
#define MCFCLASSLPSOLVER_H

#include <memory>
#include <unordered_map>
#include "../LPSolver.h"
#include "../LPSolutionWithBranchPenalty.h"
#include "../../../libs/mcf_class/MCFClass/MCFClass.h"

using namespace MCFClass_di_unipi_it;

enum class MCFClassSolver {
    RELAX_IV, NETWORK_SIMPLEX
};

typedef struct mfc_class_graph_struct {
    MCFClass::Index nodeCount = 0;
    MCFClass::Index arcCount = 0;

    MCFClass::FRow arcCapacity = nullptr;
    MCFClass::CRow arcCost = nullptr;
    MCFClass::FRow nodeDeficit = nullptr;

    MCFClass::Index_Set arcStartNodeId = nullptr;
    MCFClass::Index_Set arcEndNodeId = nullptr;

    std::unordered_map<uint32_t, uint32_t> mapJobIndexToNode;
    std::unordered_map<uint32_t, uint32_t> mapTimePointToNode;

    MCFClass::FRow arcFlow = nullptr;
    MCFClass::CRow arcReducedCost = nullptr;
    MCFClass::CRow nodePotential = nullptr;

    mfc_class_graph_struct() = default;

    ~mfc_class_graph_struct() {
        delete[] this->arcCapacity;
        this->arcCapacity = nullptr;

        delete[] this->arcCost;
        this->arcCost = nullptr;

        delete[] this->nodeDeficit;
        this->nodeDeficit = nullptr;

        delete[] this->arcStartNodeId;
        this->arcStartNodeId = nullptr;

        delete[] this->arcEndNodeId;
        this->arcEndNodeId = nullptr;

        delete[] this->arcFlow;
        this->arcFlow = nullptr;

        delete[] this->arcReducedCost;
        this->arcReducedCost = nullptr;

        delete[] this->nodePotential;
        this->nodePotential = nullptr;
    }
} mfc_class_graph_t;

class MCFClassLPSolver : public LPSolver {
private:
    const MCFClassSolver solverType;

public:
    explicit MCFClassLPSolver(const MCFClassSolver solverType);

    virtual ~MCFClassLPSolver() = default;

    LPSolution solve(const SchedulingProblem &schedProblem) const override;

    LPSolutionWithBranchPenalty solveWithBranchPenalty(const SchedulingProblem &schedProblem) const override;

    static std::unique_ptr<mfc_class_graph_t> buildGraph(const SchedulingProblem &schedProblem);

private:
    LPSolutionWithBranchPenalty
    solve_impl(const SchedulingProblem &schedProblem, bool computeBranchPenalties = false) const;

    std::vector<double> computeReducedCosts(const SchedulingProblem &schedProblem, mfc_class_graph_t &graph) const;

    std::vector<std::pair<int, branch_penalty_t>>
    computeBranchPenalties(const SchedulingProblem &schedProblem, const LPSolution &lpSolution,
                           mfc_class_graph_t &graph,
                           MCFClass *mcfSolver) const;
};


#endif //MCFCLASSLPSOLVER_H
