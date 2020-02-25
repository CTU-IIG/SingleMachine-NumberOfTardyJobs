#include <tuple>
#include <vector>
#include <algorithm>
#include <random>
#include <iostream>
#include "generator.h"

std::random_device randomDevice;
std::mt19937 mt(randomDevice());

std::uniform_int_distribution<int64_t> processingTimeDist;
std::uniform_int_distribution<int64_t> weightDist;
std::uniform_int_distribution<int64_t> duedateDist;
std::uniform_int_distribution<int64_t> jobWeightDist;

bool jobComparator(const job_t &first, const job_t &second) {
    return std::tie(first.processingTime, first.dueDate, first.deadline, first.weight) <
           std::tie(second.processingTime, second.dueDate, second.deadline, second.weight);
}

bool compareByDeadline(const job_t &first, const job_t &second) {
    return first.deadline < second.deadline;
}

bool feasibilityTest(std::vector<job_t> &jobSet) {
    bool feasible = true;
    int64_t totalProcessingTime = 0;
    std::vector<job_t> copyOfJobSet = jobSet;

    std::sort(copyOfJobSet.begin(), copyOfJobSet.end(), compareByDeadline);

    for (auto &job : copyOfJobSet) {
        if ((totalProcessingTime + job.processingTime) <= job.deadline) {
            totalProcessingTime += job.processingTime;
        } else {
            feasible = false;
            break;
        }
    }

    return feasible;
}

generator_config_t createDefaultConfig() {
    generator_config_t config = {};
    config.duedateLower = 0.3;
    config.duedateUpper = 0.7;
    config.processingTimeLower = 1;
    config.processingTimeUpper = 100;
    config.weightLower = 1;
    config.weightUpper = 100;
    config.deadlineGap = 1.1;
    config.jobCount = 10;
    config.correlated = Correlated::NONE;
    config.weightCorrelationConstant = 20;
    config.feasibleInstance = false;
    config.noDeadlines = false;
    config.twoDuedates = false;
    config.clampDeadline = false;

    return config;
}

generator_config_t parseArguments(std::vector<std::string> &arguments) {
    generator_config_t config = createDefaultConfig();

    for (size_t i = 0; i < arguments.size(); ++i) {
        if (arguments[i] == "--duedate-l" && (i + 1) < arguments.size()) {
            config.duedateLower = std::stof(arguments[i + 1]);
        } else if (arguments[i] == "--duedate-u" && (i + 1) < arguments.size()) {
            config.duedateUpper = std::stof(arguments[i + 1]);
        } else if (arguments[i] == "--correlated" && (i + 1) < arguments.size()) {
            config.correlated = static_cast<Correlated>(std::stoi(arguments[i + 1]));
        } else if (arguments[i] == "--processing-time-l" && (i + 1) < arguments.size()) {
            config.processingTimeLower = std::stoi(arguments[i + 1]);
        } else if (arguments[i] == "--processing-time-u" && (i + 1) < arguments.size()) {
            config.processingTimeUpper = std::stoi(arguments[i + 1]);
        } else if (arguments[i] == "--weight-l" && (i + 1) < arguments.size()) {
            config.weightLower = std::stoi(arguments[i + 1]);
        } else if (arguments[i] == "--weight-u" && (i + 1) < arguments.size()) {
            config.weightUpper = std::stoi(arguments[i + 1]);
        } else if (arguments[i] == "--deadline-gap" && (i + 1) < arguments.size()) {
            config.deadlineGap = std::stof(arguments[i + 1]);
        } else if (arguments[i] == "--job-count" && (i + 1) < arguments.size()) {
            config.jobCount = std::stoi(arguments[i + 1]);
        } else if (arguments[i] == "--weight-correlation" && (i + 1) < arguments.size()) {
            config.weightCorrelationConstant = std::stoi(arguments[i + 1]);
        } else if (arguments[i] == "--feasible") {
            config.feasibleInstance = true;
        } else if (arguments[i] == "--no-deadline") {
            config.noDeadlines = true;
        } else if (arguments[i] == "--two-duedates") {
            config.twoDuedates = true;
        } else if (arguments[i] == "--clamp-deadline") {
            config.clampDeadline = true;
        }
    }

    return config;
}

int main(int argc, char **argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);

    generator_config_t config = parseArguments(arguments);

    processingTimeDist = std::uniform_int_distribution<int64_t>(config.processingTimeLower, config.processingTimeUpper);
    weightDist = std::uniform_int_distribution<int64_t>(config.weightLower, config.weightUpper);

    std::vector<job_t> jobSet((unsigned long) config.jobCount);
    bool instanceIsFeasible = false;

    do {
        jobSet.clear();
        jobSet.resize((unsigned long) config.jobCount);

        for (int i = 0; i < jobSet.size(); ++i) {
            jobSet[i].processingTime = (int32_t) processingTimeDist(mt);
        }

        if (config.correlated == Correlated::STRONGLY) {
            for (int i = 0; i < jobSet.size(); ++i) {
                jobSet[i].weight = (int32_t) (jobSet[i].processingTime + config.weightCorrelationConstant);
            }
        } else if (config.correlated == Correlated::WEAKLY) {
            for (int i = 0; i < jobSet.size(); ++i) {
                jobWeightDist = std::uniform_int_distribution<int64_t>(jobSet[i].processingTime,
                                                                       jobSet[i].processingTime +
                                                                       config.weightCorrelationConstant);
                jobSet[i].weight = (int32_t) jobWeightDist(mt);
            }
        } else if (config.correlated == Correlated::NONE) {
            for (int i = 0; i < jobSet.size(); ++i) {
                jobSet[i].weight = (int32_t) weightDist(mt);
            }
        }

        int64_t totalProcessingTime = 0;
        for (int i = 0; i < jobSet.size(); ++i) {
            totalProcessingTime += jobSet[i].processingTime;
        }

        duedateDist = std::uniform_int_distribution<int64_t>(
                std::lround((double) totalProcessingTime * (double) config.duedateLower),
                std::lround((double) totalProcessingTime * (double) config.duedateUpper));

        if (config.twoDuedates) {
            int64_t firstDuedate = (int64_t) duedateDist(mt);
            int64_t secondDuedate = (int64_t) duedateDist(mt);

            for (int i = 0; i < (jobSet.size() / 2); ++i) {
                jobSet[i].dueDate = firstDuedate;
            }

            for (int i = (jobSet.size() / 2); i < jobSet.size(); ++i) {
                jobSet[i].dueDate = secondDuedate;
            }
        } else {
            for (int i = 0; i < jobSet.size(); ++i) {
                jobSet[i].dueDate = (int64_t) duedateDist(mt);
            }
        }

        for (int i = 0; i < jobSet.size(); ++i) {
            if (config.noDeadlines) {
                jobSet[i].deadline = totalProcessingTime;
            } else {
                std::uniform_int_distribution<int64_t> jobDeadlineDist(jobSet[i].dueDate, std::lround(
                        (double) totalProcessingTime * (double) config.deadlineGap));
                jobSet[i].deadline = (int64_t) jobDeadlineDist(mt);

                if(config.clampDeadline) {
                    jobSet[i].deadline = std::min(jobSet[i].deadline, totalProcessingTime);
                }
            }
        }

        instanceIsFeasible = feasibilityTest(jobSet);

        if (!config.feasibleInstance || instanceIsFeasible) {
            std::sort(jobSet.begin(), jobSet.end(), jobComparator);

            std::cout << "#job_count" << std::endl;
            std::cout << "#processing_time due_date deadline weight" << std::endl;
            std::cout << std::endl;

            std::cout << jobSet.size() << std::endl;
            for (int i = 0; i < jobSet.size(); ++i) {
                std::cout << jobSet[i].processingTime << " " << jobSet[i].dueDate << " " << jobSet[i].deadline << " "
                          << jobSet[i].weight << std::endl;
            }
        }
    } while (config.feasibleInstance && !instanceIsFeasible);
}
