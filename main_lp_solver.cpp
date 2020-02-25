#include <gurobi_c++.h>
#include <iostream>
#include <map>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <csignal>
#include <chrono>
#include "scheduling/SchedulingProblem.h"
#include "solver/lp/GurobiLPSolver.h"
#include "solver/lp/lemon/LemonLPSolver.h"
#include "solver/lp/mcf_class/MCFClassLPSolver.h"

const bool PRINT_REDUCED_COST = true;
const bool PRINT_SHADOW_PRICE = false;

std::string inputFilePath = "";
std::string outputFilePath = "";

volatile int64_t programStartTimeUS;
volatile int64_t totalRuntime = -1;
volatile double rootObjective = -1;

int threadLimit = 0;
int lpSolverType = 0;

void signalHandler(int signum) {
    int64_t endTimeUS = std::chrono::time_point_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now()).time_since_epoch().count();
    int64_t totalRunTimeDiff = endTimeUS - programStartTimeUS;

    totalRuntime = totalRunTimeDiff;

    std::cout << "#KILLED" << std::endl;
    std::cout << "#total_runtime_us";
    std::cout << " root_objective";
    std::cout << std::endl;

    std::cout << totalRuntime << " ";
    std::cout << rootObjective;
    std::cout << std::endl;

    exit(signum);
}

int main(int argc, char **argv) {
    std::signal(SIGTERM, signalHandler);
    programStartTimeUS = std::chrono::time_point_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now()).time_since_epoch().count();

    if (argc < 3) {
        std::cerr << "Invalid argument count" << std::endl;
        exit(1);
    }

    try {
        inputFilePath = argv[1];
        outputFilePath = argv[2];

        if (argc == 4) {
            threadLimit = atoi(argv[3]);
        }

        if (argc == 5) {
            lpSolverType = atoi(argv[4]);
        }

        std::ifstream inputFile;
        std::ofstream outputFile;

        inputFile.open(inputFilePath);

        SchedulingProblem schedulingProblem = SchedulingProblem::createFromStream(inputFile);

        inputFile.close();

        std::unique_ptr<LPSolver> lpSolver = nullptr;
        if (lpSolverType == 0) {
            lpSolver = std::make_unique<GurobiLPSolver>(threadLimit, false, false);
        } else if (lpSolverType == 1) {
            lpSolver = std::make_unique<LemonLPSolver>();
        } else if (lpSolverType == 2) {
            lpSolver = std::make_unique<MCFClassLPSolver>(MCFClassSolver::NETWORK_SIMPLEX);
        } else if (lpSolverType == 3) {
            lpSolver = std::make_unique<MCFClassLPSolver>(MCFClassSolver::RELAX_IV);
        } else {
            std::cerr << "Invalid LP solver type." << std::endl;
            return 1;
        }

        LPSolution solution = lpSolver->solve(schedulingProblem);

        outputFile.open(outputFilePath);
        if (solution.empty()) {
            outputFile << "-1" << std::endl;
            rootObjective = -1;
        } else {
            outputFile << solution.getObjective() << std::endl;
            rootObjective = solution.getObjective();
            std::vector<double> solutionVal = solution.getVariables();

            for (size_t i = 0; i < schedulingProblem.size(); ++i) {
                double value = solutionVal[i];

                outputFile << value;

                if ((i + 1) != schedulingProblem.size()) {
                    outputFile << " ";
                }
            }

            outputFile << std::endl;

            if (PRINT_REDUCED_COST) {
                std::vector<double> reducedCosts = solution.getReducedCosts();
                for (size_t i = 0; i < reducedCosts.size(); ++i) {
                    double reducedCost = reducedCosts[i];

                    outputFile << reducedCost;

                    if ((i + 1) != schedulingProblem.size()) {
                        outputFile << " ";
                    }
                }

                outputFile << std::endl;
            }

            if (PRINT_SHADOW_PRICE) {
                std::vector<double> shadowPrices = solution.getShadowPrices();
                for (size_t i = 0; i < shadowPrices.size(); ++i) {
                    double shadowPrice = shadowPrices[i];

                    outputFile << shadowPrice;

                    if ((i + 1) != schedulingProblem.size()) {
                        outputFile << " ";
                    }
                }

                outputFile << std::endl;
            }
        }

        outputFile.close();
    } catch (GRBException &ex) {
        std::cerr << "GRBException: " << ex.getMessage() << std::endl;
    }

    int64_t endTimeUS = std::chrono::time_point_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now()).time_since_epoch().count();
    int64_t totalRunTimeDiff = endTimeUS - programStartTimeUS;

    totalRuntime = totalRunTimeDiff;

    std::cout << "#total_runtime_us";
    std::cout << " root_objective";
    std::cout << std::endl;

    std::cout << totalRuntime << " ";
    std::cout << rootObjective;
    std::cout << std::endl;

    return 0;
}