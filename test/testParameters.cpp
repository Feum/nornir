/**
 *  Different tests on parameters.
 **/
#include "parametersLoader.hpp"
#include <algorithm>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include "../src/manager.hpp"
#include "gtest/gtest.h"

using namespace nornir;

TEST(ParametersTest, SimpleTest) {
    Parameters p = getParameters("repara");
    // Set a min just to let things work.
    p.requirements.powerConsumption = NORNIR_REQUIREMENT_MIN;

    // MaxUtilization
    p.requirements.maxUtilization = 200;
    EXPECT_EQ(p.validate(), VALIDATION_WRONG_REQUIREMENT);

    // Correct maxUtilization
    p.requirements.maxUtilization = 90;
    EXPECT_EQ(p.validate(), VALIDATION_OK);
    p.requirements.maxUtilization = NORNIR_REQUIREMENT_UNDEF;

    // MinUtilization
    p.requirements.minUtilization = -1;
    EXPECT_EQ(p.validate(), VALIDATION_WRONG_REQUIREMENT);
    p.requirements.minUtilization = NORNIR_REQUIREMENT_UNDEF;

    // Correct minUtilization
    p.requirements.minUtilization = 80;
    EXPECT_EQ(p.validate(), VALIDATION_OK);
    p.requirements.minUtilization = NORNIR_REQUIREMENT_UNDEF;

    // Overlapping min/maxUtilization
    p.requirements.minUtilization = 90;
    p.requirements.maxUtilization = 80;
    EXPECT_EQ(p.validate(), VALIDATION_WRONG_REQUIREMENT);
    p.requirements.minUtilization = NORNIR_REQUIREMENT_UNDEF;
    p.requirements.maxUtilization = NORNIR_REQUIREMENT_UNDEF;

    // Bandwidth negative
    p.requirements.bandwidth = -1;
    EXPECT_EQ(p.validate(), VALIDATION_WRONG_REQUIREMENT);

    // Bandwidth ok
    p.requirements.bandwidth = 100;
    EXPECT_EQ(p.validate(), VALIDATION_OK);
    p.requirements.bandwidth = NORNIR_REQUIREMENT_UNDEF;

    // Execution time negative
    p.requirements.executionTime = -1;
    EXPECT_EQ(p.validate(), VALIDATION_WRONG_REQUIREMENT);

    // Correct execution time but no expectedTasks
    p.requirements.executionTime = 10;
    EXPECT_EQ(p.validate(), VALIDATION_WRONG_REQUIREMENT);

    // Correct execution time and expectedTasks
    p.requirements.expectedTasksNumber = 100;
    EXPECT_EQ(p.validate(), VALIDATION_OK);
    p.requirements.executionTime = NORNIR_REQUIREMENT_UNDEF;
    p.requirements.expectedTasksNumber = NORNIR_REQUIREMENT_UNDEF;

    // Wrong latency
    p.requirements.latency = -1;
    EXPECT_EQ(p.validate(), VALIDATION_WRONG_REQUIREMENT);
    p.requirements.latency = NORNIR_REQUIREMENT_UNDEF;

    // Ok latency
    p.requirements.latency = 100;
    EXPECT_EQ(p.validate(), VALIDATION_OK);
    p.requirements.latency = NORNIR_REQUIREMENT_UNDEF;

    // Wrong power
    p.requirements.powerConsumption = -1;
    EXPECT_EQ(p.validate(), VALIDATION_WRONG_REQUIREMENT);
    p.requirements.powerConsumption = NORNIR_REQUIREMENT_UNDEF;

    // Ok power
    p.requirements.powerConsumption = 100;
    p.requirements.bandwidth = NORNIR_REQUIREMENT_MAX;
    EXPECT_EQ(p.validate(), VALIDATION_OK);
    p.requirements.powerConsumption = NORNIR_REQUIREMENT_UNDEF;

    // All requirement and one min
    p.requirements.minUtilization = 80;
    p.requirements.maxUtilization = 90;
    p.requirements.bandwidth = 10;
    p.requirements.latency = 10;
    p.requirements.powerConsumption = NORNIR_REQUIREMENT_MIN;
    EXPECT_EQ(p.validate(), VALIDATION_OK);

    // Two min/max
    p.requirements.latency = NORNIR_REQUIREMENT_MIN;
    EXPECT_EQ(p.validate(), VALIDATION_WRONG_REQUIREMENT);
    p.requirements.minUtilization = NORNIR_REQUIREMENT_UNDEF;
    p.requirements.maxUtilization = NORNIR_REQUIREMENT_UNDEF;
    p.requirements.bandwidth = NORNIR_REQUIREMENT_UNDEF;
    p.requirements.latency = NORNIR_REQUIREMENT_UNDEF;
    p.requirements.powerConsumption = NORNIR_REQUIREMENT_UNDEF;
}
