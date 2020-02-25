#include <istream>
#include <vector>
#include <fstream>
#include <gurobi_c++.h>
#include <cmath>
#include <chrono>
#include "common/types.h"
#include "scheduling/SchedulingProblem.h"
#include "reduction/JobReduce.h"
#include "algorithm/Baptiste.h"
#include <iomanip>
#include <csignal>
#include <memory>
#include <iostream>
#include <cstring>

const bool PRINT_OPTIMAL_SCHEDULE = false;
const int ILP_THRESHOLD = 14000000;

std::string inputFilePath = "";
std::string outputFilePath = "";

volatile int64_t programStartTimeUS;
volatile int64_t totalRuntime = -1;
volatile int64_t rootObjective = -1;

int threadLimit = 0;

int upperBoundMaxDepth = 8;
ILPSolverImpl ilpSolverImpl = ILPSolverImpl::ORIGINAL;
LowerBoundImpl lowerBoundImpl = LowerBoundImpl::ORIGINAL;
UpperBoundImpl upperBoundImpl = UpperBoundImpl::ORIGINAL;

std::unique_ptr<Baptiste> baptiste = nullptr;

void signalHandler(int signum) {
    if (baptiste == nullptr) {
        return;
    }

    int64_t endTimeUS = std::chrono::time_point_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now()).time_since_epoch().count();
    int64_t totalRunTimeDiff = endTimeUS - programStartTimeUS;

    totalRuntime = totalRunTimeDiff;

    std::cout << "#KILLED" << std::endl;
    std::cout << "#total_runtime_us";
    std::cout << " root_objective";
    std::cout << " root_lower_bound";
    std::cout << " root_upper_bound";
    std::cout << " total_node_visited";
    std::cout << " root_red_prob_size";
    std::cout << std::endl;

    std::cout << totalRuntime << " ";
    std::cout << rootObjective << " ";
    std::cout << baptiste->getRootLowerBound() << " ";
    std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << baptiste->getRootUpperBound()
              << " ";
    std::cout << baptiste->getTotalNodeVisited() << " ";
    std::cout << baptiste->getRootReducedProbSize();
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

        if (argc >= 4) {
            threadLimit = atoi(argv[3]);
        }

        if (argc >= 5) {
            if(strcmp(argv[4], "--ilp-orig") == 0) {
                ilpSolverImpl = ILPSolverImpl::ORIGINAL;
            } else if(strcmp(argv[4], "--ilp-imprv") == 0) {
                ilpSolverImpl = ILPSolverImpl::IMPROVED;
            } else {
                std::cerr << "Invalid argument" << std::endl;
                exit(1);
            }
        }

        if (argc >= 6) {
            if(strcmp(argv[5], "--lb-orig") == 0) {
                lowerBoundImpl = LowerBoundImpl::ORIGINAL;
            } else if(strcmp(argv[5], "--lb-imprv") == 0) {
                lowerBoundImpl = LowerBoundImpl::IMPROVED;
            } else {
                std::cerr << "Invalid argument" << std::endl;
                exit(1);
            }
        }

        if (argc >= 7) {
            if(strcmp(argv[6], "--ub-orig") == 0) {
                upperBoundImpl = UpperBoundImpl::ORIGINAL;
            } else if(strcmp(argv[6], "--ub-imprv") == 0) {
                upperBoundImpl = UpperBoundImpl::IMPROVED;
            } else {
                std::cerr << "Invalid argument" << std::endl;
                exit(1);
            }
        }

        std::cout << "Configuration:" << std::endl;
        std::cout << "ILPSolver: " << static_cast<int>(ilpSolverImpl) << std::endl;
        std::cout << "LowerBound: " << static_cast<int>(lowerBoundImpl) << std::endl;
        std::cout << "UpperBound: " << static_cast<int>(upperBoundImpl) << std::endl;
        std::cout << "Threads: " << threadLimit << std::endl;
        std::cout << "ILP threshold: " << ILP_THRESHOLD << std::endl;

        std::ifstream inputFile;
        std::ofstream outputFile;

        inputFile.open(inputFilePath);

        SchedulingProblem schedulingProblem = SchedulingProblem::createFromStream(inputFile, false);

        inputFile.close();

        baptiste = std::make_unique<Baptiste>(schedulingProblem, threadLimit, ILP_THRESHOLD, ilpSolverImpl, lowerBoundImpl, upperBoundImpl, upperBoundMaxDepth);
        BinaryILPSolution solution = baptiste->solve();

        outputFile.open(outputFilePath);
        if (solution.empty()) {
            outputFile << "-1" << std::endl;
            rootObjective = -1;
        } else {
            outputFile << solution.getObjective() << std::endl;
            rootObjective = solution.getObjective();
            std::vector<bool> solutionVal = solution.getVariables();

            for (size_t i = 0; i < schedulingProblem.size(); ++i) {
                int value = solutionVal[i];

                outputFile << value;

                if ((i + 1) != schedulingProblem.size()) {
                    outputFile << " ";
                }
            }

            outputFile << std::endl;

            if (PRINT_OPTIMAL_SCHEDULE) {
                std::vector<int> schedule = schedulingProblem.createScheduleFromSolution(solutionVal);

                for (size_t i = 0; i < schedule.size(); ++i) {
                    int scheduledJob = schedule[i];

                    outputFile << scheduledJob;

                    if ((i + 1) != schedulingProblem.size()) {
                        outputFile << " ";
                    }
                }

                outputFile << std::endl;
            }

            outputFile.close();

            outputFile << std::endl;
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
    std::cout << " root_lower_bound";
    std::cout << " root_upper_bound";
    std::cout << " total_node_visited";
    std::cout << " root_red_prob_size";
    std::cout << std::endl;

    std::cout << totalRuntime << " ";
    std::cout << rootObjective << " ";
    std::cout << baptiste->getRootLowerBound() << " ";
    std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << baptiste->getRootUpperBound()
              << " ";
    std::cout << baptiste->getTotalNodeVisited() << " ";
    std::cout << baptiste->getRootReducedProbSize();
    std::cout << std::endl;

    return 0;
}
