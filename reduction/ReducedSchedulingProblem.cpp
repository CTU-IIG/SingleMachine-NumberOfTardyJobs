#include <iostream>
#include <algorithm>
#include "ReducedSchedulingProblem.h"

SchedulingProblem ReducedSchedulingProblem::getReducedSchedulingProblemCopy() const {
    return this->reducedSchedProblem;
}

const std::vector<job_t> &ReducedSchedulingProblem::getJobs() const {
    return this->reducedSchedProblem.getJobs();
}

bool ReducedSchedulingProblem::empty() const {
    return this->reducedSchedProblem.empty();
}

size_t ReducedSchedulingProblem::size() const {
    return this->reducedSchedProblem.size();
}

const job_t &ReducedSchedulingProblem::operator[](const int index) const {
    return this->reducedSchedProblem[index];
}

BinaryILPSolution ReducedSchedulingProblem::reconstructSchedulingProblem(const BinaryILPSolution &reducedSolution) const {
    if (/*reducedSolution.empty() || */reducedSolution.getVariables().size() != this->reducedSchedProblem.size()) {
        std::cerr << "[ReducedSchedulingProblem::reconstructSchedulingProblem]: Invalid reducedSolution" << std::endl;
        return BinaryILPSolution();
    }

    size_t newVarCount = this->reducedSchedProblem.size() + this->reduceJobs.size();
    if(newVarCount != this->originalSchedProblem.size()) {
        std::cerr << "[ReducedSchedulingProblem::reconstructSchedulingProblem]: Invalid state of instance" << std::endl;
    }

    int64_t newObjective = reducedSolution.getObjective() +
                           this->computeReducedJobsObjective();
    std::vector<bool> newVars;

    bool insertedVars[newVarCount] = {};
    newVars.resize(newVarCount);
    for (const std::pair<int, bool> &reduceJob : this->reduceJobs) {
        newVars[reduceJob.first] = reduceJob.second;
        insertedVars[reduceJob.first] = true;
    }

    std::vector<bool> reducedSolVars = reducedSolution.getVariables();

    size_t reducedSolIndex = 0;
    for (size_t i = 0; i < newVarCount; ++i) {
        if (!insertedVars[i]) {
            newVars[i] = reducedSolVars[reducedSolIndex];
            insertedVars[i] = true;
            ++reducedSolIndex;
        }
    }

    if (reducedSolIndex != reducedSolVars.size()) {
        std::cerr
                << "[ReducedSchedulingProblem::reconstructSchedulingProblem]: reducedSolIndex != reducedSolVars.size()"
                << std::endl;
    }

    return BinaryILPSolution(newObjective, newVars);
}

int64_t ReducedSchedulingProblem::computeReducedJobsObjective() const {
    int64_t reducedJobsObj = 0;

    for (const std::pair<int, bool> &reduceJob : this->reduceJobs) {
        if (reduceJob.second) {
            reducedJobsObj += this->originalSchedProblem[reduceJob.first].weight;
        }
    }

    return reducedJobsObj;
}

int64_t ReducedSchedulingProblem::getReducedJobsObjective() const {
    if (this->reducedJobsObjective == -1) {
        this->reducedJobsObjective = this->computeReducedJobsObjective();
    }

    return this->reducedJobsObjective;
}

ReducedSchedulingProblem::ReducedSchedulingProblem(const SchedulingProblem &originalSchedProblem,
                                                   const SchedulingProblem &reducedSchedProblem,
                                                   const std::vector<std::pair<int, bool>> &reduceJobs)
        : originalSchedProblem(originalSchedProblem), reducedSchedProblem(reducedSchedProblem),
          reduceJobs(reduceJobs) {}

bool ReducedSchedulingProblem::isFeasible() const {
    return this->reducedSchedProblem.isFeasible();
}

std::vector<int> ReducedSchedulingProblem::getNonReducedJobs() const {
    std::vector<std::pair<int, bool>> reducedJobsCopy = this->reduceJobs;
    std::sort(reducedJobsCopy.begin(), reducedJobsCopy.end());

    std::vector<int> nonReducedJobs;

    int lastInsertedIndex = -1;
    for (const std::pair<int, bool> &reduceJob : reducedJobsCopy) {
        for (int i = (lastInsertedIndex + 1); i < reduceJob.first; ++i) {
            nonReducedJobs.push_back(i);
        }

        lastInsertedIndex = reduceJob.first;
    }

    for (size_t i = (lastInsertedIndex + 1); i < this->originalSchedProblem.size(); ++i) {
        nonReducedJobs.push_back(i);
    }

    return nonReducedJobs;
}

std::vector<std::pair<int, bool>> ReducedSchedulingProblem::getReduceJobs() const {
    return this->reduceJobs;
}

int ReducedSchedulingProblem::mapReducedIndexToOriginal(int reducedIndex) const {
    if(reducedIndex >= this->reducedSchedProblem.size()) {
        std::cerr << "[ReducedSchedulingProblem::mapReducedIndexToOriginal]: Invalid index" << std::endl;
        return -1;
    }

    if(this->indexMapping.empty() && !this->reducedSchedProblem.empty()) {
        size_t newVarCount = this->reducedSchedProblem.size() + this->reduceJobs.size();
        if(newVarCount != this->originalSchedProblem.size()) {
            std::cerr << "[ReducedSchedulingProblem::mapReducedIndexToOriginal]: Invalid state of instance" << std::endl;
        }

        bool insertedVars[newVarCount] = {};
        for (const std::pair<int, bool> &reduceJob : this->reduceJobs) {
            insertedVars[reduceJob.first] = true;
        }

        size_t reducedSolIndex = 0;
        for (size_t i = 0; i < newVarCount; ++i) {
            if (!insertedVars[i]) {
                this->indexMapping[reducedSolIndex] = i;
                insertedVars[i] = true;
                ++reducedSolIndex;
            }
        }
    }

    return this->indexMapping[reducedIndex];
}

int64_t ReducedSchedulingProblem::getNonZeroCoefCountInModel() const {
    return this->reducedSchedProblem.getNonZeroCoefCountInModel();
}

BinaryILPSolution ReducedSchedulingProblem::transformSolutionToReduced(const BinaryILPSolution &originalSolution) const {
    if (/*reducedSolution.empty() || */originalSolution.getVariables().size() != this->originalSchedProblem.size()) {
        std::cerr << "[ReducedSchedulingProblem::transformSolutionToReduced]: Invalid originalSolution" << std::endl;
        return BinaryILPSolution();
    }

    size_t newVarCount = this->reducedSchedProblem.size();

    int64_t newObjective = originalSolution.getObjective() - this->computeReducedJobsObjective();
    std::vector<bool> newVars;

    bool fixedVars[originalSolution.size()] = {};
    newVars.resize(newVarCount);
    for (const std::pair<int, bool> &reduceJob : this->reduceJobs) {
        fixedVars[reduceJob.first] = true;
    }

    std::vector<bool> origSolVals = originalSolution.getVariables();

    size_t reducedSolIndex = 0;
    for (size_t i = 0; i < origSolVals.size(); ++i) {
        if (!fixedVars[i]) {
            newVars[reducedSolIndex] = origSolVals[i];
            fixedVars[i] = true;
            ++reducedSolIndex;
        }
    }

    if (reducedSolIndex != this->reducedSchedProblem.size()) {
        std::cerr << "[ReducedSchedulingProblem::transformSolutionToReduced]: reducedSolIndex != reducedSolVars.size()" << std::endl;
    }

    return BinaryILPSolution(newObjective, newVars);
}
