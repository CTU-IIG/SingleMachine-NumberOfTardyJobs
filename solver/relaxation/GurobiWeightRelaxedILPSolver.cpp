#include <iostream>
#include <gurobi_c++.h>
#include <cmath>
#include "GurobiWeightRelaxedILPSolver.h"

GurobiWeightRelaxedILPSolver::GurobiWeightRelaxedILPSolver(uint32_t threadCount, bool printOutput, bool usePresolve)
        : threadCount(threadCount),
          printOutput(printOutput),
          usePresolve(usePresolve) {}

BinaryILPSolution GurobiWeightRelaxedILPSolver::solve(const SchedulingProblem &schedProblem) const {
    return solve_impl(schedProblem);
}

BinaryILPSolution GurobiWeightRelaxedILPSolver::solve_impl(const SchedulingProblem &schedProblem) const {
    if (schedProblem.empty()) {
        std::cerr << "[GurobiWeightRelaxedILPSolver::solve_impl]: schedProblem is empty" << std::endl;
        return BinaryILPSolution();
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

    try {
        for (size_t i = 0; i < schedProblem.size(); ++i) {
            xVar[i] = model.addVar(0.0, 1.0, 1.0, GRB_BINARY);
        }

        for (int64_t timePoint : timePointsOrdered) {
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

            model.addConstr(constrainExpr <= timePoint);
        }

        model.optimize();
    } catch (GRBException &ex) {
        std::cerr << "GurobiWeightRelaxedILPSolver::solve_impl]: GRBException: " << ex.getMessage() << std::endl;
        return BinaryILPSolution();
    }

    if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
        std::cerr << "GurobiWeightRelaxedILPSolver::solve_impl]: not optimal solution found" << std::endl;
        return BinaryILPSolution();
    }

    int64_t objective = (int64_t) model.get(GRB_DoubleAttr_ObjVal);
    std::vector<bool> varValues;

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        varValues.push_back((bool) std::lround(xVar[i].get(GRB_DoubleAttr_X)));
    }

    return BinaryILPSolution(objective, varValues);
}

BinaryILPSolution
GurobiWeightRelaxedILPSolver::solve(const SchedulingProblem &schedProblem, const BinaryILPSolution &lowerBound) const {
    std::cerr << "[GurobiWeightRelaxedILPSolver::solve]: Unimplemented method" << std::endl;
    exit(1);
    return BinaryILPSolution();
}
