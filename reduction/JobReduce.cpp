#include <iostream>
#include <unordered_set>
#include "JobReduce.h"
#include "../common/util.h"
#include "../solver/lp/lemon/LemonLPSolver.h"
#include "../solver/ilp/OriginILPSolver.h"
#include "../heuristic/lowerbound/OriginalLowerBound.h"
#include "../solver/lp/mcf_class/MCFClassLPSolver.h"
#include "../heuristic/lowerbound/ImprovedLowerBound.h"
#include "../heuristic/upperbound/ImprovedUpperBound.h"
#include "../heuristic/upperbound/ImprovedUpperBoundOpt.h"


SchedulingProblem JobReduce::reduceJob_impl(const SchedulingProblem &schedProblem, int jobIndex, bool jobEarly) const {
    std::vector<job_t> reducedInstance;

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        if(i != jobIndex) {
            reducedInstance.push_back(schedProblem[i]);
        }
    }

    int Dj = ((jobEarly) ? schedProblem[jobIndex].dueDate : schedProblem[jobIndex].deadline);

    for (job_t &job : reducedInstance) {
        if(job.dueDate <= Dj) {
            int possibleLatestStartTimeJ = Dj - schedProblem[jobIndex].processingTime;

            if(possibleLatestStartTimeJ < job.dueDate) {
                job.dueDate = possibleLatestStartTimeJ;
            }
        } else {
            job.dueDate -= schedProblem[jobIndex].processingTime;
        }

        if(job.deadline <= Dj) {
            int possibleLatestStartTimeJ = Dj - schedProblem[jobIndex].processingTime;

            if(possibleLatestStartTimeJ < job.deadline) {
                job.deadline = possibleLatestStartTimeJ;
            }
        } else {
            job.deadline -= schedProblem[jobIndex].processingTime;
        }
    }

    return SchedulingProblem(reducedInstance, schedProblem.hasDeadlines());
}

SchedulingProblem JobReduce::reduceJobs_impl(const SchedulingProblem &schedProblem, const std::vector<std::pair<int, bool>> &reduceJobs) const {
    std::vector<std::pair<int, bool>> reduceJobsCopy = reduceJobs;
    SchedulingProblem schedProblemCopy = schedProblem;

    for (size_t i = 0; i < reduceJobsCopy.size(); ++i) {
        schedProblemCopy = this->reduceJob_impl(schedProblemCopy, reduceJobsCopy[i].first, reduceJobsCopy[i].second);

        for (size_t j = i; j < reduceJobsCopy.size(); ++j) {
            --reduceJobsCopy[j].first;
        }
    }

    return schedProblemCopy;
}

SchedulingProblem JobReduce::reduceJob_impl(const SchedulingProblem &schedProblem, const std::vector<int> &jobIndexes, const std::vector<bool> &jobEarly) const {
    if(jobIndexes.size() != jobEarly.size()) {
        std::cerr << "[JobReduce::reduceJobs]: Invalid input" << std::endl;
    }

    std::vector<std::pair<int, bool>> reduceJobs;
    for (size_t i = 0; i < jobIndexes.size(); ++i) {
        reduceJobs.emplace_back(jobIndexes[i], jobEarly[i]);
    }

    return this->reduceJobs_impl(schedProblem, reduceJobs);
}

ReducedSchedulingProblem JobReduce::reduceJob(const SchedulingProblem &schedProblem, int jobIndex, bool jobEarly) const {
    SchedulingProblem reducedSchedProblem = this->reduceJob_impl(schedProblem, jobIndex, jobEarly);
    std::vector<std::pair<int, bool>> reducedJobs;
    reducedJobs.emplace_back(jobIndex, jobEarly);

    return ReducedSchedulingProblem(schedProblem, reducedSchedProblem, reducedJobs);
}

ReducedSchedulingProblem JobReduce::reduceJobs(const SchedulingProblem &schedProblem, const std::vector<std::pair<int, bool>> &reduceJobs, bool allowExtendedReduction) const {
    SchedulingProblem reducedSchedProblem = this->reduceJobs_impl(schedProblem, reduceJobs);
    std::vector<std::pair<int, bool>> reducedJobs = reduceJobs;

    if(allowExtendedReduction) {
        int nextReduceJobsCount = 0;

        do {
            std::vector<std::pair<int, bool>> nextReduceJobs = this->reductionPreprocessing(reducedSchedProblem);

            if(!nextReduceJobs.empty()) {
                ReducedSchedulingProblem redSchedProbTemp(schedProblem, reducedSchedProblem, reducedJobs);

                for (const std::pair<int, bool> &reduceJob : nextReduceJobs) {
                    reducedJobs.emplace_back(redSchedProbTemp.mapReducedIndexToOriginal(reduceJob.first), reduceJob.second);
                }

                std::sort(reducedJobs.begin(), reducedJobs.end());

                reducedSchedProblem = this->reduceJobs_impl(schedProblem, reducedJobs);
            }

            nextReduceJobsCount = nextReduceJobs.size();
        } while(nextReduceJobsCount > 0);

        return ReducedSchedulingProblem(schedProblem, reducedSchedProblem, reducedJobs);
    }

    return ReducedSchedulingProblem(schedProblem, reducedSchedProblem, reducedJobs);
}

std::vector<std::pair<int, bool>> JobReduce::getForceJobByDominance(const SchedulingProblem &schedProblem, int jobIndex, bool jobEarly) const {
    std::vector<std::pair<int, bool>> reducedJobs;

    if(jobEarly) {
        for (size_t i = 0; i < schedProblem.size(); ++i) {
            bool domCond = (schedProblem[i].processingTime <= schedProblem[jobIndex].processingTime)
                           && (schedProblem[i].dueDate >= schedProblem[jobIndex].dueDate)
                           && (schedProblem[i].deadline <= schedProblem[jobIndex].deadline)
                           && (schedProblem[i].weight >= schedProblem[jobIndex].weight);

            if(domCond) {
                int strictEq = (schedProblem[i].processingTime != schedProblem[jobIndex].processingTime)
                               || (schedProblem[i].dueDate != schedProblem[jobIndex].dueDate )
                               || (schedProblem[i].deadline != schedProblem[jobIndex].deadline)
                               || (schedProblem[i].weight != schedProblem[jobIndex].weight);

                if(strictEq) {
                    reducedJobs.emplace_back(i, jobEarly);
                }
            }
        }
    } else {
        for (size_t j = 0; j < schedProblem.size(); ++j) {
            bool domCond = (schedProblem[jobIndex].processingTime <= schedProblem[j].processingTime)
                        && (schedProblem[jobIndex].dueDate >= schedProblem[j].dueDate)
                        && (schedProblem[jobIndex].deadline <= schedProblem[j].deadline)
                        && (schedProblem[jobIndex].weight >= schedProblem[j].weight);

            if(domCond) {
                int strictEq = (schedProblem[jobIndex].processingTime != schedProblem[j].processingTime)
                            || (schedProblem[jobIndex].dueDate != schedProblem[j].dueDate)
                            || (schedProblem[jobIndex].deadline != schedProblem[j].deadline)
                            || (schedProblem[jobIndex].weight != schedProblem[j].weight);

                if(strictEq) {
                    reducedJobs.emplace_back(j, jobEarly);
                }
            }
        }
    }

    return reducedJobs;
}

std::vector<std::pair<int, bool>> JobReduce::getForceJobByDominance(const SchedulingProblem &schedProblem, const std::vector<std::pair<int, bool>> &reduceJobs) const {
    std::vector<std::pair<int, bool>> newReduceJobsVec;
    std::unordered_set<std::pair<int, bool>> newReduceJobsSet;

    for (const std::pair<int, bool> &reduceJob : reduceJobs) {
        std::vector<std::pair<int, bool>> currentReducedJobs = this->getForceJobByDominance(schedProblem, reduceJob.first, reduceJob.second);
        if(!currentReducedJobs.empty()) {
            newReduceJobsSet.insert(currentReducedJobs.begin(), currentReducedJobs.end());
        }
    }

    newReduceJobsVec.insert(newReduceJobsVec.end(), newReduceJobsSet.begin(), newReduceJobsSet.end());
    std::sort(newReduceJobsVec.begin(), newReduceJobsVec.end());
    return newReduceJobsVec;
}

ReducedSchedulingProblem JobReduce::reduceJobsByVarFixing(const SchedulingProblem &schedProblem, const BinaryILPSolution &ilpSolution, const LPSolution &lpSolution) const {
    if(schedProblem.size() != ilpSolution.getVariables().size() || schedProblem.size() != lpSolution.getVariables().size()) {
        std::cerr << "[JobReduce::reduceJobsByVarFixing]: Size error" << std::endl;
    }

    std::vector<std::pair<int, bool>> reducedJobs;
    std::vector<double> reducedCosts = lpSolution.getReducedCosts();
    double upperBound = lpSolution.getObjective();
    double lowerBound = ilpSolution.getObjective();

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        if((reducedCosts[i] != 0) && (lpSolution[i] == 0 || lpSolution[i] == 1)) {
            if((reducedCosts[i] < 0) && ((upperBound + reducedCosts[i]) <= lowerBound)) {
                //R1
                reducedJobs.emplace_back(i, false);
            }

            if((reducedCosts[i] > 0) && ((upperBound - reducedCosts[i]) <= lowerBound)) {
                //R2
                reducedJobs.emplace_back(i, true);
            }
        }
    }

    JobReduce jobReduceID1(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);
    std::unique_ptr<LPSolver> lemonLpSolverID1 = nullptr;
    if(this->upperBoundImpl == UpperBoundImpl::ORIGINAL) {
        lemonLpSolverID1 = std::make_unique<MCFClassLPSolver>(MCFClassSolver::RELAX_IV);
    } else if(this->upperBoundImpl == UpperBoundImpl::IMPROVED) {
        lemonLpSolverID1 = std::make_unique<ImprovedUpperBoundOpt>(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);
    } else {
        std::cerr << "[JobReduce::reduceJobsByVarFixing]: Unknow lower bound implementation" << std::endl;
        exit(1);
    }

    std::unique_ptr<OriginalLowerBound> lowerBoundID1 = nullptr;
    if(this->lowerBoundImpl == LowerBoundImpl::ORIGINAL) {
        lowerBoundID1 = std::make_unique<OriginalLowerBound>(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);
    } else if(this->lowerBoundImpl == LowerBoundImpl::IMPROVED) {
        lowerBoundID1 = std::make_unique<ImprovedLowerBound>(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);
    } else {
        std::cerr << "[JobReduce::reduceJobsByVarFixing]: Unknow lower bound implementation" << std::endl;
        exit(1);
    }

    ReducedSchedulingProblem reducedSchedProblemID1Wrapper = jobReduceID1.reduceJobs(schedProblem, reducedJobs, true);
    reducedJobs = reducedSchedProblemID1Wrapper.getReduceJobs();

    SchedulingProblem reducedSchedProblemID1 = reducedSchedProblemID1Wrapper.getReducedSchedulingProblemCopy();
    LPSolutionWithBranchPenalty reducedUpperBoundID1 = lemonLpSolverID1->solveWithBranchPenalty(reducedSchedProblemID1);
    BinaryILPSolution reducedLowerBoundID1 = lowerBoundID1->solve(reducedSchedProblemID1, reducedUpperBoundID1);


    double reducedUpperBound = reducedUpperBoundID1.getObjective();
    double reducedLowerBound = reducedLowerBoundID1.getObjective();

    bool performNextIteration = false;
    std::vector<std::pair<int, branch_penalty_t>> reducedBranchPenalties = reducedUpperBoundID1.getBranchPenalties();

    for (size_t branchPenaltyIndex = 0; branchPenaltyIndex < reducedBranchPenalties.size(); ++branchPenaltyIndex) {
        int jobIndex = reducedBranchPenalties[branchPenaltyIndex].first;
        branch_penalty_t branchPenalty = reducedBranchPenalties[branchPenaltyIndex].second;

        if ((0 < reducedUpperBoundID1[jobIndex]) && (reducedUpperBoundID1[jobIndex] < 1)) {
            JobReduce jobReduceRID2(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);

            if((reducedUpperBound + branchPenalty.upperChange <= reducedLowerBound) && (reducedUpperBound + branchPenalty.lowerChange <= reducedLowerBound)) {
                std::cerr << "[JobReduce::reduceJobsByVarFixing]: Inconsistent state, double fixed basis variable, caused by numerical errors. Skipping variable: " << jobIndex << std::endl;
            } else if(reducedUpperBound + branchPenalty.upperChange <= reducedLowerBound) {
                int remapedIndex = reducedSchedProblemID1Wrapper.mapReducedIndexToOriginal(jobIndex);
                reducedJobs.emplace_back(remapedIndex, false);
                ReducedSchedulingProblem reducedSchedProblem = jobReduceRID2.reduceJobs(schedProblem, reducedJobs);

                if(!reducedSchedProblem.isFeasible()) {
                    std::cerr << "[JobReduce::reduceJobsByVarFixing]]: Inconsistent state: Infeasible solution by rounding down variable: " << jobIndex << std::endl;
                    reducedJobs.pop_back();
                } else {
                    performNextIteration = true;
                }
            } else if(reducedUpperBound + branchPenalty.lowerChange <= reducedLowerBound) {
                int remapedIndex = reducedSchedProblemID1Wrapper.mapReducedIndexToOriginal(jobIndex);
                reducedJobs.emplace_back(remapedIndex, true);
                ReducedSchedulingProblem reducedSchedProblem = jobReduceRID2.reduceJobs(schedProblem, reducedJobs);

                if(!reducedSchedProblem.isFeasible()) {
                    std::cerr << "[JobReduce::reduceJobsByVarFixing]]: Inconsistent state: Infeasible solution by rounding up variable: " << jobIndex << std::endl;
                    reducedJobs.pop_back();
                } else {
                    performNextIteration = true;
                }
            }
        }
    }

    std::sort(reducedJobs.begin(), reducedJobs.end());

    JobReduce jobReduceFinID1(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);
    ReducedSchedulingProblem reducedSchedProblemFinID1 = jobReduceFinID1.reduceJobs(schedProblem, reducedJobs, true);
    std::vector<int> nonReducedJobs = reducedSchedProblemFinID1.getNonReducedJobs();

    std::cerr << "performNextIteration: " << performNextIteration << std::endl;

    if (performNextIteration && !reducedSchedProblemFinID1.empty()) {
        SchedulingProblem schedProblemNext = reducedSchedProblemFinID1.getReducedSchedulingProblemCopy();

        std::unique_ptr<OriginalLowerBound> lowerBoundID2 = nullptr;
        if(this->lowerBoundImpl == LowerBoundImpl::ORIGINAL) {
            lowerBoundID2 = std::make_unique<OriginalLowerBound>(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);
        } else if(this->lowerBoundImpl == LowerBoundImpl::IMPROVED) {
            lowerBoundID2 = std::make_unique<ImprovedLowerBound>(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);
        } else {
            std::cerr << "[JobReduce::reduceJobsByVarFixing]: Unknow lower bound implementation" << std::endl;
            exit(1);
        }

        std::unique_ptr<LPSolver> lpSolverNext = nullptr;
        if(this->upperBoundImpl == UpperBoundImpl::ORIGINAL) {
            lpSolverNext = std::make_unique<MCFClassLPSolver>(MCFClassSolver::RELAX_IV);
        } else if(this->upperBoundImpl == UpperBoundImpl::IMPROVED) {
            lpSolverNext = std::make_unique<ImprovedUpperBoundOpt>(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);
        } else {
            std::cerr << "[JobReduce::reduceJobsByVarFixing]: Unknow lower bound implementation" << std::endl;
            exit(1);
        }

        JobReduce jobReduceNext(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);

        LPSolutionWithBranchPenalty lpSolutionNext = lpSolverNext->solveWithBranchPenalty(schedProblemNext);
        BinaryILPSolution ilpSolutionNext = lowerBoundID2->solve(schedProblemNext, lpSolutionNext);

        ReducedSchedulingProblem reducedSchedProblemNext = jobReduceNext.reduceJobsByVarFixing(schedProblemNext, ilpSolutionNext, lpSolutionNext);
        std::vector<std::pair<int, bool>> reducedJobsNext = reducedSchedProblemNext.getReduceJobs();

        for (size_t i = 0; i < reducedJobsNext.size(); ++i) {
            reducedJobs.emplace_back(nonReducedJobs[reducedJobsNext[i].first], reducedJobsNext[i].second);
        }

        std::sort(reducedJobs.begin(), reducedJobs.end());

        JobReduce jobReduceOut(this->ilpSolverImpl, this->lowerBoundImpl, this->upperBoundImpl, this->upperBoundMaxDepth);
        return jobReduceOut.reduceJobs(schedProblem, reducedJobs, true);
    }

    return reducedSchedProblemFinID1;
}

std::vector<std::pair<int, bool>> JobReduce::reductionPreprocessing(const SchedulingProblem &schedProblem) const {
    std::vector<std::pair<int, bool>> reduceJobs;

    for (size_t i = 0; i < schedProblem.size(); ++i) {
        if(schedProblem[i].dueDate < 0) {
            reduceJobs.emplace_back(i, false);
        } else if(schedProblem[i].dueDate >= 0 && schedProblem[i].processingTime == 0) {
            reduceJobs.emplace_back(i, true);
        } else if(schedProblem[i].dueDate >= 0 && schedProblem[i].dueDate == schedProblem[i].deadline && schedProblem[i].processingTime <= schedProblem[i].dueDate) {
            reduceJobs.emplace_back(i, true);
        }
    }

    return reduceJobs;
}

JobReduce::JobReduce(ILPSolverImpl ilpSolverImpl, LowerBoundImpl lowerBoundImpl, UpperBoundImpl upperBoundImpl, int upperBoundMaxDepth) : ilpSolverImpl(ilpSolverImpl), lowerBoundImpl(lowerBoundImpl), upperBoundImpl(upperBoundImpl), upperBoundMaxDepth(upperBoundMaxDepth) {}
