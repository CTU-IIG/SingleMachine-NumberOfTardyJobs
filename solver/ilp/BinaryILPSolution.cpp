#include "BinaryILPSolution.h"

#include <utility>
#include <limits>

int64_t BinaryILPSolution::getObjective() const {
    return objective;
}

const std::vector<bool> &BinaryILPSolution::getVariables() const {
    return variables;
}

BinaryILPSolution::BinaryILPSolution() : objective(std::numeric_limits<int64_t>::min()) {}

BinaryILPSolution::BinaryILPSolution(int64_t objective) : objective(objective) {}

BinaryILPSolution::BinaryILPSolution(int64_t objective, std::vector<bool> variables) : objective(objective), variables(
        std::move(variables)) {}

bool BinaryILPSolution::empty() const {
    return this->variables.empty();
}

bool BinaryILPSolution::operator[](const int index) const {
    return this->variables[index];
}

size_t BinaryILPSolution::size() const {
    return this->variables.size();
}

void BinaryILPSolution::setObjective(int64_t objective) {
    this->objective = objective;
}

void BinaryILPSolution::increaseObjective(int64_t increaseBy) {
    this->objective += increaseBy;
}
