#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <numeric>
#include "SchedulingProblem.h"
#include "../common/util.h"

SchedulingProblem::SchedulingProblem(const std::vector<job_t> &jobs, bool containDeadlines) : jobs(jobs), containDeadlines(containDeadlines) {}

SchedulingProblem::SchedulingProblem(const SchedulingProblem &schedulingProblem) {
    this->jobs = schedulingProblem.jobs;
    this->orderedTimePoints = schedulingProblem.orderedTimePoints;
    this->totalProcessingTime = schedulingProblem.totalProcessingTime;
    this->containDeadlines = schedulingProblem.containDeadlines;
}

int64_t SchedulingProblem::computeTotalProcessingTime() const {
    int64_t processingTimeSum = 0;

    for (const job_t &job : this->jobs) {
        processingTimeSum += job.processingTime;
    }

    return processingTimeSum;
}

std::vector<int64_t> SchedulingProblem::computeOrderedTimePoints() const {
    std::unordered_set<int64_t> timePointsSet;
    std::vector<int64_t> orderedTimePointsVec;

    for (const job_t &job : this->jobs) {
        timePointsSet.insert(job.dueDate);
        timePointsSet.insert(job.deadline);
    }

    orderedTimePointsVec.insert(orderedTimePointsVec.end(), timePointsSet.begin(), timePointsSet.end());
    std::sort(orderedTimePointsVec.begin(), orderedTimePointsVec.end());

    return orderedTimePointsVec;
}

const std::vector<job_t> &SchedulingProblem::getJobs() const {
    return jobs;
}

const std::vector<int64_t> &SchedulingProblem::getOrderedTimePoints() const {
    if (this->orderedTimePoints.empty() && !this->jobs.empty()) {
        this->orderedTimePoints = this->computeOrderedTimePoints();
    }

    return this->orderedTimePoints;
}

int64_t SchedulingProblem::getTotalProcessingTime() const {
    if (this->totalProcessingTime == -1) {
        this->totalProcessingTime = this->computeTotalProcessingTime();
    }

    return this->totalProcessingTime;
}

bool SchedulingProblem::empty() const {
    return this->jobs.empty();
}

size_t SchedulingProblem::size() const {
    return this->jobs.size();
}

const job_t &SchedulingProblem::operator[](const int index) const {
    return this->jobs[index];
}

void SchedulingProblem::sortByDeadline() {
    std::sort(this->jobs.begin(), this->jobs.end(), SchedulingProblem::compareByDeadline);
}

bool SchedulingProblem::compareByDeadline(const job_t &first, const job_t &second) {
    return first.deadline < second.deadline;
}

SchedulingProblem SchedulingProblem::createFromStream(std::istream &input, bool clampDeadline) {
    std::ios::sync_with_stdio(false);

    std::string inputLine;

    std::getline(input, inputLine);
    std::getline(input, inputLine);
    std::getline(input, inputLine);

    std::vector<job_t> jobs;
    uint32_t jobCount = 0;

    input >> jobCount;

    bool sameDeadline = true;
    int64_t totalProcessingTime = 0;
    for (size_t i = 0; i < jobCount; ++i) {
        job_t job = {};

        input >> job.processingTime >> job.dueDate >> job.deadline >> job.weight;
        totalProcessingTime += job.processingTime;

        jobs.push_back(job);

        if(sameDeadline && i > 0 && job.deadline != jobs[0].deadline) {
            sameDeadline = false;
        }
    }

    if (clampDeadline) {
        for (job_t &job : jobs) {
            if (job.deadline > totalProcessingTime) {
                job.deadline = totalProcessingTime;
            }
        }
    }

    return SchedulingProblem(jobs, !sameDeadline);
}

void SchedulingProblem::setDeadline(int index, int32_t deadline) {
    if (0 <= index && index < this->jobs.size()) {
        this->jobs[index].deadline = deadline;
    }
}

BinaryILPSolution SchedulingProblem::constructILPSolutionFromLPSolution(const LPSolution &lpSolution) const {
    if (this->size() != lpSolution.getVariables().size()) {
        std::cerr << "[SchedulingProblem::constructILPSolutionFromLPSolution]: this->size() != lpSolution.size()"
                  << std::endl;
        return BinaryILPSolution();
    }

    SchedulingProblem schedProblemCopy = *(this);
    std::vector<double> lpSolVar = lpSolution.getVariables();
    for (size_t i = 0; i < schedProblemCopy.size(); ++i) {
        if (lpSolVar[i] >= 1) {
            schedProblemCopy.setDeadline(i, schedProblemCopy[i].dueDate);
        }
    }

    return schedProblemCopy.findFeasibleSolutionByEDF();
}

BinaryILPSolution SchedulingProblem::findFeasibleSolutionByEDF() const {
    SchedulingProblem schedProblemCopy = *(this);
    schedProblemCopy.sortByDeadline();

    int64_t processingTimeSum = 0;
    int64_t objective = 0;
    std::vector<bool> solVar;

    for (size_t i = 0; i < schedProblemCopy.size(); ++i) {
        if ((processingTimeSum + schedProblemCopy[i].processingTime) <= schedProblemCopy[i].deadline) {
            processingTimeSum += schedProblemCopy[i].processingTime;

            if ((processingTimeSum + schedProblemCopy[i].processingTime) <= schedProblemCopy[i].dueDate) {
                objective += schedProblemCopy[i].weight;
                solVar.push_back(true);
            } else {
                solVar.push_back(false);
            }
        } else {
            return BinaryILPSolution();
        }
    }

    return BinaryILPSolution(objective, solVar);
}

bool SchedulingProblem::isFeasible() const {
    return !this->findFeasibleSolutionByEDF().empty();
}

bool SchedulingProblem::isFeasibleSolution(const std::vector<bool> &solution) const {
    return !this->createScheduleFromSolution(solution).empty();
}

bool SchedulingProblem::isFeasibleSolution(const BinaryILPSolution &solution) const {
    return !this->createScheduleFromSolution(solution.getVariables()).empty();
}

bool SchedulingProblem::isFeasibleSolution(const std::vector<bool> &solution, const std::vector<int> &presorted) const {
    return !this->createScheduleFromSolution(solution, presorted).empty();
}

std::vector<int> SchedulingProblem::createScheduleFromSolution(const std::vector<bool> &solution) const {
    std::vector<job_t> jobsCopy = this->jobs;
    std::vector<int> schedule(this->jobs.size());

    if (this->jobs.size() != solution.size()) {
        std::cerr << "[SchedulingProblem::createScheduleFromSolution]: this->jobs.size() != solution.size()"
                  << std::endl;
        return std::vector<int>();
    }

    for (size_t i = 0; i < jobsCopy.size(); ++i) {
        if (solution[i]) {
            jobsCopy[i].deadline = jobsCopy[i].dueDate;
        }
    }

    std::iota(schedule.begin(), schedule.end(), 0);
    sort(schedule.begin(), schedule.end(), [&jobsCopy](const size_t &indexFirst, const size_t &indexSecond) {
        return jobsCopy[indexFirst].deadline < jobsCopy[indexSecond].deadline;
    });

    int64_t processingTimeSum = 0;

    for (const int &scheduleIndex : schedule) {
        job_t scheduledJob = jobsCopy[scheduleIndex];

        if ((processingTimeSum + scheduledJob.processingTime) <= scheduledJob.deadline) {
            processingTimeSum += scheduledJob.processingTime;
        } else {
            schedule.clear();
            break;
        }
    }

    return schedule;
}

std::vector<int> SchedulingProblem::createScheduleFromSolution(const std::vector<bool> &solution, const std::vector<int> &presorted) const {
    std::vector<job_t> jobsCopy = this->jobs;
    std::vector<int> schedule = presorted;

    if (this->jobs.size() != solution.size()) {
        std::cerr << "[SchedulingProblem::createScheduleFromSolution]: this->jobs.size() != solution.size()" << std::endl;
        return std::vector<int>();
    }

    for (size_t i = 0; i < jobsCopy.size(); ++i) {
        if (solution[i]) {
            jobsCopy[i].deadline = jobsCopy[i].dueDate;
        }
    }

    sortings::InsertionSort::sort(schedule.begin(), schedule.end(), [&jobsCopy](const size_t &indexFirst, const size_t &indexSecond) {
        return jobsCopy[indexFirst].deadline < jobsCopy[indexSecond].deadline;
    });

    int64_t processingTimeSum = 0;

    for (const int &scheduleIndex : schedule) {
        job_t scheduledJob = jobsCopy[scheduleIndex];

        if ((processingTimeSum + scheduledJob.processingTime) <= scheduledJob.deadline) {
            processingTimeSum += scheduledJob.processingTime;
        } else {
            schedule.clear();
            break;
        }
    }

    return schedule;
}

int64_t SchedulingProblem::getNonZeroCoefCountInModel() const {
    if (this->nonZeroCoefCountInModel == -1) {
        int64_t nonZeroCount = 0;
        std::vector<int64_t> timePointsOrdered = this->getOrderedTimePoints();

        for (int64_t &timePoint : timePointsOrdered) {
            for (const job_t &job : this->jobs) {
                if((job.deadline > timePoint) && (job.dueDate <= timePoint) ) {
                    nonZeroCount++;
                }
            }
        }

        this->nonZeroCoefCountInModel = nonZeroCount;
    }

    return this->nonZeroCoefCountInModel;
}

bool SchedulingProblem::hasDeadlines() const {
    return this->containDeadlines;
}
