#include <iostream>
#include "Baptiste.h"
#include "../solver/lp/lemon/LemonLPSolver.h"
#include "../heuristic/lowerbound/OriginalLowerBound.h"
#include "../solver/ilp/OriginILPSolver.h"
#include "../reduction/JobReduce.h"
#include "../solver/lp/mcf_class/MCFClassLPSolver.h"
#include "../heuristic/lowerbound/ImprovedLowerBound.h"
#include "../solver/ilp/ILPGridSolver.h"
#include "../heuristic/upperbound/ImprovedUpperBoundOpt.h"

Baptiste::Baptiste(const SchedulingProblem &schedulingProblem, uint32_t threadCount, int64_t ilpModelNonZeroLimit, ILPSolverImpl ilpSolverImpl, LowerBoundImpl lowerBoundImpl, UpperBoundImpl upperBoundImpl, int upperBoundMaxDepth) : schedulingProblem(schedulingProblem), threadCount(threadCount), ilpModelNonZeroLimit(ilpModelNonZeroLimit), ilpSolverImpl(ilpSolverImpl), lowerBoundImpl(lowerBoundImpl), upperBoundImpl(upperBoundImpl), upperBoundMaxDepth(upperBoundMaxDepth) {}

BinaryILPSolution Baptiste::solve() {
    return this->branchAndBounds(this->schedulingProblem, 0);
}

BinaryILPSolution Baptiste::branchAndBounds(const SchedulingProblem &schedProblem, const int64_t objectiveAdd) {
    std::cerr << "[branchAndBounds]: Iteration, globalBestObjective: " << this->globalBestObjective<< std::endl;

    if(schedProblem.empty()) {
        std::cerr << "[branchAndBounds]: schedProblem is empty" << std::endl;
        return BinaryILPSolution(0);
    }

    if(!schedProblem.isFeasible()) {
        std::cerr << "[branchAndBounds]: schedProblem not have feasible solution" << std::endl;
        return BinaryILPSolution(-100000000000);
    }

    if(schedProblem.size() <= 100) {
        ++this->totalNodeVisited;
        std::cerr << "[branchAndBounds]: schedProblem.size() = 100" << std::endl;

        BinaryILPSolution ilpSolution;
        if(this->ilpSolverImpl == ILPSolverImpl::ORIGINAL) {
            OriginILPSolver ilpSolver(this->threadCount, false, false);
            ilpSolution = ilpSolver.solve(schedProblem);
        } else if(this->ilpSolverImpl == ILPSolverImpl::IMPROVED) {
            ILPGridSolver ilpSolver(this->threadCount, false, false);
            ilpSolution = ilpSolver.solve(schedProblem);
        }

        if((ilpSolution.getObjective() + objectiveAdd) > this->globalBestObjective) {
            this->globalBestObjective = (ilpSolution.getObjective() + objectiveAdd);
        }

        return ilpSolution;
    }

    std::unique_ptr<LPSolver> classLpSolver = nullptr;
    if(this->upperBoundImpl == UpperBoundImpl::ORIGINAL) {
        classLpSolver = std::make_unique<MCFClassLPSolver>(MCFClassSolver::RELAX_IV);
    } else if(this->upperBoundImpl == UpperBoundImpl::IMPROVED) {
        classLpSolver = std::make_unique<ImprovedUpperBoundOpt>(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);
    } else {
        std::cerr << "[branchAndBounds]: Unknow lower bound implementation" << std::endl;
        exit(1);
    }

    LPSolution classLpSolution = classLpSolver->solve(schedProblem);

    if(classLpSolution.empty()) {
        std::cerr << "[branchAndBounds]: classLpSolution.empty()" << std::endl;
        return BinaryILPSolution(-100000000000);
    }

    ++this->totalNodeVisited;

    if(this->totalNodeVisited == 1) {
        this->rootUpperBound = classLpSolution.getObjective();
    }

    if(classLpSolution.getNonIntegerValueCount() == 0) {
        std::cerr << "[branchAndBounds]: Solved by LP" << std::endl;
        return classLpSolution.convertToBinarySolution();
    }

    if((classLpSolution.getObjective() + objectiveAdd) < (double)(this->globalBestObjective + 1)) {
        std::cerr << "[branchAndBounds]: UB is worst than this->globalBestObjective + 1" << std::endl;
        return schedProblem.constructILPSolutionFromLPSolution(classLpSolution);
    }

    std::unique_ptr<OriginalLowerBound> lowerBoundHeuristic = nullptr;
    if(this->lowerBoundImpl == LowerBoundImpl::ORIGINAL) {
        lowerBoundHeuristic = std::make_unique<OriginalLowerBound>(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);
    } else if(this->lowerBoundImpl == LowerBoundImpl::IMPROVED) {
        lowerBoundHeuristic = std::make_unique<ImprovedLowerBound>(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);
    } else {
        std::cerr << "[JobReduce::reduceJobsByVarFixing]: Unknow lower bound implementation" << std::endl;
        exit(1);
    }

    BinaryILPSolution lowerBoundHSolution = lowerBoundHeuristic->solve(schedProblem, classLpSolution);

    if(this->totalNodeVisited == 1) {
        this->rootLowerBound = lowerBoundHSolution.getObjective();
    }

    if((lowerBoundHSolution.getObjective() + objectiveAdd) > this->globalBestObjective) {
        this->globalBestObjective = (lowerBoundHSolution.getObjective() + objectiveAdd);
    }

    if(classLpSolution.getObjective() < (double)(lowerBoundHSolution.getObjective() + 1)) {
        std::cerr << "[branchAndBounds]: UB is worst than LB + 1" << std::endl;
        return lowerBoundHSolution;
    }

    JobReduce jobReduce(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);
    ReducedSchedulingProblem reducedSchedProblem = jobReduce.reduceJobsByVarFixing(schedProblem, lowerBoundHSolution, classLpSolution);

    std::cerr << "[branchAndBounds]: Reduced set: " << reducedSchedProblem.size() << std::endl;
    std::cerr << "[branchAndBounds]: Fixed jobSet: " << reducedSchedProblem.getReduceJobs().size() << std::endl;

    if(this->totalNodeVisited == 1) {
        this->rootReducedProbSize = reducedSchedProblem.size();
    }

    if(reducedSchedProblem.empty()) {
        std::cerr << "[branchAndBounds]: reducedSchedProblem is empty" << std::endl;
        BinaryILPSolution simpleSolution = BinaryILPSolution(0);
        std::cout << reducedSchedProblem.reconstructSchedulingProblem(simpleSolution).getObjective() << std::endl;
        return reducedSchedProblem.reconstructSchedulingProblem(simpleSolution);
    }

    if(!reducedSchedProblem.isFeasible()) {
        std::cerr << "[branchAndBounds]: reducedSchedProblem not have feasible solution" << std::endl;
        return BinaryILPSolution(-100000000000);
    }

    if(reducedSchedProblem.getNonZeroCoefCountInModel() <= ilpModelNonZeroLimit) {
        std::cerr << "[branchAndBounds]: reducedJobSet <= ilpModelNonZeroLimit" << std::endl;
        SchedulingProblem schedProblemReducedJobs = reducedSchedProblem.getReducedSchedulingProblemCopy();

        BinaryILPSolution ilpSolutionOfReduced;
        if(this->ilpSolverImpl == ILPSolverImpl::ORIGINAL) {
            OriginILPSolver ilpSolver(this->threadCount, false, false);
            ilpSolutionOfReduced = ilpSolver.solve(schedProblemReducedJobs);
        } else if(this->ilpSolverImpl == ILPSolverImpl::IMPROVED) {
            ILPGridSolver ilpSolver(this->threadCount, false, false);
            BinaryILPSolution lowerBoundHSolutionReduced = reducedSchedProblem.transformSolutionToReduced(lowerBoundHSolution);
            ilpSolutionOfReduced = ilpSolver.solve(schedProblemReducedJobs, lowerBoundHSolutionReduced);
        }

        if(ilpSolutionOfReduced.empty()) {
            std::cerr << "[branchAndBounds]: reducedJobSet <= ilpModelNonZeroLimit but it is infeasible" << std::endl;
            return BinaryILPSolution(-100000000000);
        }

        BinaryILPSolution reconstructedSolution = reducedSchedProblem.reconstructSchedulingProblem(ilpSolutionOfReduced);

        if((reconstructedSolution.getObjective() + objectiveAdd) > this->globalBestObjective) {
            this->globalBestObjective = (reconstructedSolution.getObjective() + objectiveAdd);
        }

        if(reconstructedSolution.getObjective() > lowerBoundHSolution.getObjective()) {
            return reconstructedSolution;
        } else {
            return lowerBoundHSolution;
        }
    }

    std::unique_ptr<LPSolver> classLpSolverR2 = nullptr;
    if(this->upperBoundImpl == UpperBoundImpl::ORIGINAL) {
        classLpSolverR2 = std::make_unique<MCFClassLPSolver>(MCFClassSolver::RELAX_IV);
    } else if(this->upperBoundImpl == UpperBoundImpl::IMPROVED) {
        classLpSolverR2 = std::make_unique<ImprovedUpperBoundOpt>(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);
    } else {
        std::cerr << "[branchAndBounds]: Unknow lower bound implementation" << std::endl;
        exit(1);
    }

    LPSolutionWithBranchPenalty classLpSolutionBPR2 = classLpSolverR2->solveWithBranchPenalty(reducedSchedProblem.getReducedSchedulingProblemCopy());

    if(classLpSolutionBPR2.getNonIntegerValueCount() == 0) {
        std::cerr << "[branchAndBounds]: Solved by LP after reduction" << std::endl;
        BinaryILPSolution convertedSolutionToILP = classLpSolutionBPR2.convertToBinarySolution();
        return reducedSchedProblem.reconstructSchedulingProblem(convertedSolutionToILP);
    }

    if(classLpSolutionBPR2.getBranchPenalties().size() != classLpSolutionBPR2.getNonIntegerValueCount() ) {
        std::cerr << "[branchAndBounds]: Inconsistent state of application, because classLpSolutionBPR2.getBranchPenalties() and classLpSolutionBPR2.getNonIntegerValueCount() are differ.";
        exit(111);
    }

    int maxMinPseudoCostIndex = classLpSolutionBPR2.getMaxMinPseudoCost().first;

    if(maxMinPseudoCostIndex == -1) {
        std::cerr << "[branchAndBounds]: Inconsistent state of application, because maxMinPseudoCostIndex == -1";
        exit(112);
    }


    JobReduce jobReduceDownID1(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);
    JobReduce jobReduceUpID2(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);

    SchedulingProblem schedProblemBeforeBran = reducedSchedProblem.getReducedSchedulingProblemCopy();
    ReducedSchedulingProblem schedulingProblemBranchDown = jobReduceDownID1.reduceJob(schedProblemBeforeBran, maxMinPseudoCostIndex, false);
    ReducedSchedulingProblem schedulingProblemBranchUp = jobReduceUpID2.reduceJob(schedProblemBeforeBran, maxMinPseudoCostIndex, true);

    BinaryILPSolution branchDownSolution = this->branchAndBounds(schedulingProblemBranchDown.getReducedSchedulingProblemCopy(), objectiveAdd + reducedSchedProblem.getReducedJobsObjective() + 0);
    int64_t branchDownObjective = branchDownSolution.getObjective() + reducedSchedProblem.getReducedJobsObjective() + 0;
    if(branchDownObjective > lowerBoundHSolution.getObjective()) {
        BinaryILPSolution reconstrBranchDownSol = reducedSchedProblem.reconstructSchedulingProblem(schedulingProblemBranchDown.reconstructSchedulingProblem(branchDownSolution));

        if((reconstrBranchDownSol.getObjective() + objectiveAdd) > this->globalBestObjective) {
            this->globalBestObjective = (reconstrBranchDownSol.getObjective() + objectiveAdd);
        }

        lowerBoundHSolution = reconstrBranchDownSol;
    }


    if(schedProblemBeforeBran[maxMinPseudoCostIndex].processingTime <= schedProblemBeforeBran[maxMinPseudoCostIndex].dueDate) {
        BinaryILPSolution branchUpSolution = this->branchAndBounds(schedulingProblemBranchUp.getReducedSchedulingProblemCopy(), objectiveAdd + reducedSchedProblem.getReducedJobsObjective() + schedProblemBeforeBran[maxMinPseudoCostIndex].weight);
        int64_t branchUpObjective = branchUpSolution.getObjective() + reducedSchedProblem.getReducedJobsObjective() + schedProblemBeforeBran[maxMinPseudoCostIndex].weight;
        if(branchUpObjective > lowerBoundHSolution.getObjective()) {
            BinaryILPSolution reconstrBranchUpSol = reducedSchedProblem.reconstructSchedulingProblem(schedulingProblemBranchUp.reconstructSchedulingProblem(branchUpSolution));

            if((reconstrBranchUpSol.getObjective() + objectiveAdd) > this->globalBestObjective) {
                this->globalBestObjective = (reconstrBranchUpSol.getObjective() + objectiveAdd);
            }

            lowerBoundHSolution = reconstrBranchUpSol;
        }
    }

    return lowerBoundHSolution;
}

volatile int64_t Baptiste::getGlobalBestObjective() const {
    return globalBestObjective;
}

volatile int64_t Baptiste::getTotalNodeVisited() const {
    return totalNodeVisited;
}

volatile int64_t Baptiste::getRootLowerBound() const {
    return rootLowerBound;
}

volatile double Baptiste::getRootUpperBound() const {
    return rootUpperBound;
}

volatile int Baptiste::getRootReducedProbSize() const {
    return rootReducedProbSize;
}
