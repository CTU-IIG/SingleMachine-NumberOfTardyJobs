#ifndef GENERATOR_H
#define GENERATOR_H

enum Correlated {
    NONE = 0, WEAKLY = 1, STRONGLY = 2
};

typedef struct {
    int64_t processingTime;
    int64_t dueDate;
    int64_t deadline;
    int64_t weight;
} job_t;

typedef struct {
    float duedateLower = 0.3;
    float duedateUpper = 0.7;
    int64_t processingTimeLower = 1;
    int64_t processingTimeUpper = 100;
    int64_t weightLower = 1;
    int64_t weightUpper = 100;
    float deadlineGap = 1.1;
    int64_t jobCount = 10;
    Correlated correlated = Correlated::NONE;
    int64_t weightCorrelationConstant = 20;
    bool feasibleInstance = false;
    bool noDeadlines = false;
    bool twoDuedates = false;
    bool clampDeadline = false;
} generator_config_t;

bool jobComparator(const job_t &first, const job_t &second);

bool compareByDeadline(const job_t &first, const job_t &second);

bool feasibilityTest(std::vector<job_t> &jobSet);

generator_config_t parseArguments(std::vector<std::string> &arguments);

generator_config_t createDefaultConfig();

#endif
