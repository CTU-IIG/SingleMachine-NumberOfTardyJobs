#ifndef BINARYILPSOLUTION_H
#define BINARYILPSOLUTION_H


#include <vector>
#include <cstdint>
#include <cstddef>

class BinaryILPSolution {
private:
    int64_t objective;
    std::vector<bool> variables;

public:
    virtual ~BinaryILPSolution() = default;

    BinaryILPSolution();

    BinaryILPSolution(int64_t objective);

    explicit BinaryILPSolution(int64_t objective, std::vector<bool> variables);

    int64_t getObjective() const;

    const std::vector<bool> &getVariables() const;

    bool empty() const;

    size_t size() const;

    bool operator[](const int index) const;

    void setObjective(int64_t objective);

    void increaseObjective(int64_t increaseBy);
};


#endif //BINARYILPSOLUTION_H
