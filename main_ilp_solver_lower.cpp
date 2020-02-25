#include <gurobi_c++.h>
#include <iostream>
#include <map>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <chrono>
#include <csignal>
#include "scheduling/SchedulingProblem.h"
#include "solver/ilp/OriginILPSolver.h"
#include "solver/relaxation/GurobiWeightRelaxedILPSolver.h"
#include "solver/relaxation/GurobiMaxProcessingTimeILPSolver.h"
#include "solver/ilp/ILPGridSolver.h"
#include "solver/lp/LPSolver.h"
#include "solver/lp/GurobiLPSolver.h"

std::string inputFilePath = "";
std::string outputFilePath = "";

volatile int64_t programStartTimeUS;
volatile int64_t totalRuntime = -1;
volatile int64_t rootObjective = -1;
volatile int64_t minEarly = -1;
volatile int64_t lpEarly = -1;

int threadLimit = 0;

void signalHandler(int signum) {
    int64_t endTimeUS = std::chrono::time_point_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now()).time_since_epoch().count();
    int64_t totalRunTimeDiff = endTimeUS - programStartTimeUS;

    totalRuntime = totalRunTimeDiff;

    std::cout << "#KILLED" << std::endl;
    std::cout << "#total_runtime_us";
    std::cout << " root_objective";
    std::cout << " min_early";
    std::cout << " min_lp_early";
    std::cout << std::endl;

    std::cout << totalRuntime << " ";
    std::cout << rootObjective << " ";
    std::cout << minEarly << " ";
    std::cout << lpEarly;
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

        std::ifstream inputFile;
        std::ofstream outputFile;

        inputFile.open(inputFilePath);

        SchedulingProblem schedulingProblem = SchedulingProblem::createFromStream(inputFile);

        inputFile.close();

        GurobiMaxProcessingTimeILPSolver originIlpSolver(threadLimit, false, false);
        BinaryILPSolution solution = originIlpSolver.solve(schedulingProblem);

        outputFile.open(outputFilePath);
        if (solution.empty()) {
            outputFile << "-1" << std::endl;
            rootObjective = -1;
        } else {
            outputFile << solution.getObjective() << std::endl;
            rootObjective = solution.getObjective();
            std::vector<bool> solutionVal = solution.getVariables();
            minEarly = 0;

            for (size_t i = 0; i < schedulingProblem.size(); ++i) {
                int value = solutionVal[i];

                outputFile << value;

                if ((i + 1) != schedulingProblem.size()) {
                    outputFile << " ";
                }

                if(value == 1) {
                    ++minEarly;
                }
            }

            outputFile << std::endl;
        }

        outputFile.close();


        GurobiLPSolver lpSolver(threadLimit, false, false);
        LPSolution lpSolution = lpSolver.solve(schedulingProblem);

        if (!lpSolution.empty()) {
            lpEarly = 0;

            for (size_t i = 0; i < schedulingProblem.size(); ++i) {
                if(lpSolution[i] == 1) {
                    ++lpEarly;
                }
            }
        }
    } catch (GRBException &ex) {
        std::cerr << "GRBException: " << ex.getMessage() << std::endl;
    }

    int64_t endTimeUS = std::chrono::time_point_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now()).time_since_epoch().count();
    int64_t totalRunTimeDiff = endTimeUS - programStartTimeUS;

    totalRuntime = totalRunTimeDiff;

    std::cout << "#total_runtime_us";
    std::cout << " root_objective";
    std::cout << " min_early";
    std::cout << " min_lp_early";
    std::cout << std::endl;

    std::cout << totalRuntime << " ";
    std::cout << rootObjective << " ";
    std::cout << minEarly << " ";
    std::cout << lpEarly;
    std::cout << std::endl;

    return 0;
}