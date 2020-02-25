#ifndef IMPROVEDUPPERBOUNDOPT_H
#define IMPROVEDUPPERBOUNDOPT_H


#include "../../solver/lp/LPSolver.h"
#include "../../solver/lp/mcf_class/MCFClassLPSolver.h"

class ImprovedUpperBoundOpt : public LPSolver {
private:
    ILPSolverImpl ilpSolverImpl = ILPSolverImpl::ORIGINAL;
    LowerBoundImpl lowerBoundImpl = LowerBoundImpl::ORIGINAL;
    UpperBoundImpl upperBoundImpl = UpperBoundImpl::ORIGINAL;
    int startDepth = 1;
public:
    ImprovedUpperBoundOpt() = delete;

    explicit ImprovedUpperBoundOpt(ILPSolverImpl ilpSolverImpl, LowerBoundImpl lowerBoundImpl, UpperBoundImpl upperBoundImpl, int startDepth);

    virtual ~ImprovedUpperBoundOpt() = default;

    LPSolution solve(const SchedulingProblem &schedProblem) const override;

    LPSolutionWithBranchPenalty solveWithBranchPenalty(const SchedulingProblem &schedProblem) const override;

private:
    double recursiveSearch(MCFClass *mcfSolver, const SchedulingProblem &schedProblem, const LPSolution &lpSolution, int depth = 1) const;

    std::vector<double> computeReducedCosts(const SchedulingProblem &schedProblem, mfc_class_graph_t &graph) const;

    std::vector<std::pair<int, branch_penalty_t>> computeBranchPenalties(const SchedulingProblem &schedProblem, const LPSolution &lpSolution, MCFClass *mcfSolver) const;

    int selectJobByProcessingTime(const SchedulingProblem &schedProblem, const LPSolution &lpSolution, bool max = true) const;

    int selectJobByDueDate(const SchedulingProblem &schedProblem, const LPSolution &lpSolution, bool max = true) const;

    int selectJobByDeadline(const SchedulingProblem &schedProblem, const LPSolution &lpSolution, bool max = true) const;

    int selectJobByMaxMinPseudoCost(const std::vector<std::pair<int, branch_penalty_t>> &branchingPenalties) const;
};


#endif //IMPROVEDUPPERBOUNDOPT_H
