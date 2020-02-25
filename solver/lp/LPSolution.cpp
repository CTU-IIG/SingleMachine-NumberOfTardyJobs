#include <limits>
#include <utility>
#include <cmath>
#include "LPSolution.h"

LPSolution::LPSolution() : objective(std::numeric_limits<double>::lowest()) {}


LPSolution::LPSolution(double objective, std::vector<double> variables, std::vector<double> reducedCosts)
        : objective(objective), variables(std::move(variables)), reducedCosts(std::move(reducedCosts)) {}

LPSolution::LPSolution(double objective, std::vector<double> variables, std::vector<double> reducedCosts,
                       std::vector<double> shadowPrices) : objective(objective), variables(std::move(variables)),
                                                           reducedCosts(std::move(reducedCosts)),
                                                           shadowPrices(std::move(shadowPrices)) {}

double LPSolution::getObjective() const {
    return objective;
}

const std::vector<double> &LPSolution::getVariables() const {
    return variables;
}

LPSolution::LPSolution(double objective, std::vector<double> variables) : objective(objective), variables(
        std::move(variables)) {}

const std::vector<double> &LPSolution::getReducedCosts() const {
    return reducedCosts;
}

const std::vector<double> &LPSolution::getShadowPrices() const {
    return shadowPrices;
}

bool LPSolution::empty() const {
    return this->variables.empty();
}

BinaryILPSolution LPSolution::convertToBinarySolution() const {
    int64_t intObjective = this->objective;
    std::vector<bool> binVariables;

    for (double variable : this->variables) {
        binVariables.push_back(std::lround(variable));
    }

    return BinaryILPSolution(intObjective, binVariables);
}

size_t LPSolution::size() const {
    return this->variables.size();
}

const double &LPSolution::operator[](const int index) const {
    return this->variables[index];
}

size_t LPSolution::getNonIntegerValueCount() const {
    size_t nonIntegerVariableCount = 0;

    for (double variable : this->variables) {
        if((0 < variable) && (variable < 1)) {
            nonIntegerVariableCount++;
        }
    }

    return nonIntegerVariableCount;
}
