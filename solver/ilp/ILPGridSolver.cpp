#include <gurobi_c++.h>
#include <cmath>
#include "ILPGridSolver.h"
#include "../lp/GurobiLPSolver.h"
#include "../relaxation/GurobiWeightRelaxedILPSolver.h"

ILPGridSolver::ILPGridSolver(uint32_t threadCount, bool printOutput, bool usePresolve) : threadCount(threadCount),
                                                                                             printOutput(printOutput),
                                                                                             usePresolve(usePresolve) {}

BinaryILPSolution ILPGridSolver::solve(const SchedulingProblem &schedProblem) const {
    return solve_impl(schedProblem, BinaryILPSolution());
}

BinaryILPSolution ILPGridSolver::solve_impl(const SchedulingProblem &schedProblem, BinaryILPSolution bestSolution) const {
    LPSolutionWithEarlyRange lpSolution = this->getEarlyRange(schedProblem);

    if(lpSolution.getEarlyRangeLower() == lpSolution.getEarlyRangeUpper()) {
        return lpSolution.convertToBinarySolution();
    }

    if(bestSolution.empty()) {
        bestSolution = schedProblem.constructILPSolutionFromLPSolution(lpSolution);
    }

    GurobiWeightRelaxedILPSolver relaxedWeighSolver(this->threadCount, this->printOutput, this->usePresolve);
    BinaryILPSolution relaxedWeightSolution = relaxedWeighSolver.solve(schedProblem);
    int relaxedWeightObjective = relaxedWeightSolution.getObjective();

    int earlyRangeLower = std::max(lpSolution.getEarlyRangeLower(), 1);
    int earlyRangeUpper = std::min(lpSolution.getEarlyRangeUpper(), std::min(relaxedWeightObjective, (int)schedProblem.size()));

    std::cerr << "[ILPGridSolver::solve_impl]: GurobiWeightRelaxedILPSolver: " << relaxedWeightObjective << ", earlyRangeUpper: " << lpSolution.getEarlyRangeUpper() << std::endl;

    if(earlyRangeLower > 1) {
        BinaryILPSolution newSolution = this->solveWithEarlyConstrain(schedProblem, earlyRangeLower - 1, EarlyJobConstrainType::LOWER_EQUAL, bestSolution.getObjective());

        if(!newSolution.empty() && newSolution.getObjective() > bestSolution.getObjective()) {
            bestSolution = newSolution;
        }
    }

    for (int k = earlyRangeLower; k <= earlyRangeUpper; ++k) {
        BinaryILPSolution newSolution = this->solveWithEarlyConstrain(schedProblem, k, EarlyJobConstrainType::EQUAL, bestSolution.getObjective());

        if(!newSolution.empty() && newSolution.getObjective() > bestSolution.getObjective()) {
            bestSolution = newSolution;
        }
    }

    if(relaxedWeightObjective > earlyRangeUpper) {
        BinaryILPSolution newSolution = this->solveWithEarlyConstrain(schedProblem, earlyRangeUpper + 1, EarlyJobConstrainType::GREATER_EQUAL, bestSolution.getObjective());

        if(!newSolution.empty() && newSolution.getObjective() > bestSolution.getObjective()) {
            bestSolution = newSolution;
        }
    }

    return bestSolution;
}

LPSolutionWithEarlyRange ILPGridSolver::getEarlyRange(const SchedulingProblem &schedProblem) const {
    GurobiLPSolver gurobiLpSolver(this->threadCount, this->printOutput, this->usePresolve);
    LPSolution solution = gurobiLpSolver.solve(schedProblem);

    if(solution.empty()) {
        return LPSolutionWithEarlyRange(0, 0);
    } else {
        int earlyRangeLower = 0;
        int earlyRangeUpper = 0;

        std::vector<double> solutionVals = solution.getVariables();
        for (size_t i = 0; i < schedProblem.size(); ++i) {
            double value = solutionVals[i];

            if(value > 0) {
                earlyRangeUpper++;

                if(value == 1) {
                    earlyRangeLower++;
                }
            }
        }

        return LPSolutionWithEarlyRange(solution.getObjective(), solution.getVariables(), earlyRangeLower, earlyRangeUpper);
    }
}

BinaryILPSolution ILPGridSolver::solveWithEarlyConstrain(const SchedulingProblem &schedProblem, int earlyJobCount,
                                                         EarlyJobConstrainType constrainType, int64_t currentBestObjextive) const {
    std::vector<int64_t> timePointsOrdered = schedProblem.getOrderedTimePoints();
    GRBVar xVar[schedProblem.size()];

    GRBEnv environment;
    GRBModel model(environment);

    model.set(GRB_DoubleParam_MIPGap, 0);
    model.set(GRB_DoubleParam_PSDTol, 0);
    model.set(GRB_IntParam_Threads, this->threadCount);
    model.set(GRB_IntParam_Presolve, ((this->usePresolve) ? GRB_PRESOLVE_AUTO : GRB_PRESOLVE_OFF));
    model.set(GRB_IntParam_OutputFlag, (int) this->printOutput);

    model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

    try {
        for (size_t i = 0; i < schedProblem.size(); ++i) {
            xVar[i] = model.addVar(0.0, 1.0, schedProblem[i].weight, GRB_BINARY);
        }

        for (int64_t timePoint : timePointsOrdered) {
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

        GRBLinExpr constrainEarlyJobs;
        for (size_t i = 0; i < schedProblem.size(); ++i) {
            constrainEarlyJobs += xVar[i];
        }

        if(constrainType == EarlyJobConstrainType::EQUAL) {
            model.addConstr(constrainEarlyJobs == earlyJobCount);
        } else if(constrainType == EarlyJobConstrainType::LOWER_EQUAL) {
            model.addConstr(constrainEarlyJobs <= earlyJobCount);
        } else if(constrainType == EarlyJobConstrainType::GREATER_EQUAL) {
            model.addConstr(constrainEarlyJobs >= earlyJobCount);
        } else if(constrainType == EarlyJobConstrainType::NONE) {

        }

        model.set(GRB_DoubleParam_Cutoff, (double)(currentBestObjextive + 1));
        model.optimize();
    } catch (GRBException &ex) {
        std::cerr << "[ILPGridSolver::solveWithEarlyConstrain]: GRBException: " << ex.getMessage() << std::endl;
        return BinaryILPSolution();
    }

    if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
        return BinaryILPSolution();
    } else if (model.get(GRB_IntAttr_Status) == GRB_CUTOFF) {
        return BinaryILPSolution();
    }

    int64_t objective = (int64_t) model.get(GRB_DoubleAttr_ObjVal);
    std::vector<bool> varValues;

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        varValues.push_back((bool) std::lround(xVar[i].get(GRB_DoubleAttr_X)));
    }

    return BinaryILPSolution(objective, varValues);
}

BinaryILPSolution ILPGridSolver::solve(const SchedulingProblem &schedProblem, const BinaryILPSolution &lowerBound) const {
    return solve_impl(schedProblem, lowerBound);
}
