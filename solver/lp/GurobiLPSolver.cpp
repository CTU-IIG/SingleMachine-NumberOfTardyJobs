#include <iostream>
#include <gurobi_c++.h>
#include "GurobiLPSolver.h"

GurobiLPSolver::GurobiLPSolver(uint32_t threadCount, bool printOutput, bool usePresolve) : threadCount(threadCount),
                                                                                           printOutput(printOutput),
                                                                                           usePresolve(usePresolve) {}

LPSolution GurobiLPSolver::solve(const SchedulingProblem &schedulingProblem) const {
    return solve_impl(schedulingProblem);
}

LPSolution GurobiLPSolver::solve_impl(const SchedulingProblem &schedProblem) const {
    if (schedProblem.empty()) {
        std::cerr << "[GurobiLPSolver::solve_impl]: schedProblem is empty" << std::endl;
        return LPSolution();
    }

    GRBEnv environment;
    GRBModel model(environment);

    model.set(GRB_DoubleParam_MIPGap, 0);
    model.set(GRB_DoubleParam_PSDTol, 0);
    model.set(GRB_IntParam_Threads, this->threadCount);
    model.set(GRB_IntParam_Presolve, ((this->usePresolve) ? GRB_PRESOLVE_AUTO : GRB_PRESOLVE_OFF));
    model.set(GRB_IntParam_OutputFlag, (int) this->printOutput);

    model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

    std::vector<int64_t> timePointsOrdered = schedProblem.getOrderedTimePoints();
    GRBVar xVar[schedProblem.size()];
    GRBConstr constrains[timePointsOrdered.size()];

    try {
        for (size_t i = 0; i < schedProblem.size(); ++i) {
            xVar[i] = model.addVar(0.0, 1.0, schedProblem[i].weight, GRB_CONTINUOUS);
        }

        for (size_t i = 0; i < timePointsOrdered.size(); ++i) {
            int64_t timePoint = timePointsOrdered[i];

            // Fix for negative timePoint
            if (timePoint < 0) {
                timePoint = 0;
            }

            GRBLinExpr constrainExpr;

            for (size_t j = 0; j < schedProblem.size(); ++j) {
                if (schedProblem[j].deadline <= timePoint) {
                    constrainExpr += schedProblem[j].processingTime;
                }

                if ((schedProblem[j].deadline > timePoint) && (schedProblem[j].dueDate <= timePoint)) {
                    constrainExpr += xVar[j] * schedProblem[j].processingTime;
                }
            }

            constrains[i] = model.addConstr(constrainExpr <= timePoint);
        }

        model.optimize();
    } catch (GRBException &ex) {
        std::cerr << "GRBException: " << ex.getMessage() << std::endl;
        return LPSolution();
    }

    if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
        std::cerr << "[GurobiLPSolver::solve_impl]: not optimal solution found" << std::endl;
        return LPSolution();
    }

    double objective = model.get(GRB_DoubleAttr_ObjVal);
    std::vector<double> varValues;
    std::vector<double> reducedPrices;
    std::vector<double> shadowCosts;

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        varValues.push_back(xVar[i].get(GRB_DoubleAttr_X));
    }

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        reducedPrices.push_back(xVar[i].get(GRB_DoubleAttr_RC));
    }

    for (size_t i = 0; i < timePointsOrdered.size(); ++i) {
        shadowCosts.push_back(constrains[i].get(GRB_DoubleAttr_Pi));
    }

    std::cerr << "[GurobiLPSolver::solve_impl]: UB: " << objective << std::endl;

    return LPSolution(objective, varValues, reducedPrices, shadowCosts);
}

LPSolutionWithBranchPenalty GurobiLPSolver::solveWithBranchPenalty(const SchedulingProblem &schedProblem) const {
    std::cerr << "[GurobiLPSolver::solveWithBranchPenalty]: Unimplemented method" << std::endl;
    exit(1);
    return LPSolutionWithBranchPenalty();
}
