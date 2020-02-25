#ifndef LPSOLUTION_H
#define LPSOLUTION_H


#include <vector>
#include <cstddef>
#include "../ilp/BinaryILPSolution.h"

class LPSolution {
private:
    double objective;
    std::vector<double> variables;
    std::vector<double> reducedCosts;
    std::vector<double> shadowPrices;

public:
    virtual ~LPSolution() = default;

    LPSolution();

    explicit LPSolution(double objective, std::vector<double> variables);

    LPSolution(double objective, std::vector<double> variables, std::vector<double> reducedCosts);

    LPSolution(double objective, std::vector<double> variables, std::vector<double> reducedCosts,
               std::vector<double> shadowPrices);

    double getObjective() const;

    const std::vector<double> &getVariables() const;

    const std::vector<double> &getReducedCosts() const;

    const std::vector<double> &getShadowPrices() const;

    bool empty() const;

    size_t size() const;

    const double &operator[](const int index) const;

    BinaryILPSolution convertToBinarySolution() const;

    size_t getNonIntegerValueCount() const;
};


#endif //LPSOLUTION_H
