#ifndef TYPES_H
#define TYPES_H

enum class LPSolverType {GUROBI, LEMON, MCFClass_SIMPLEX, MCFClass_RELAX_IV};

enum class ILPSolverImpl {ORIGINAL, IMPROVED};

enum class LowerBoundImpl {ORIGINAL, IMPROVED};

enum class UpperBoundImpl {ORIGINAL, IMPROVED};

typedef struct {
    int32_t processingTime;
    int32_t dueDate;
    int32_t deadline;
    int32_t weight;
} job_t;
#endif
