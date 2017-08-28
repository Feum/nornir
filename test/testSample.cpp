/**
 *  Different tests on samples operators.
 **/
#include <algorithm>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include "../src/nornir.hpp"
#include "gtest/gtest.h"


// getMaximumBandwidth() function.
TEST(SamplesTest, MaximumBandwidth) {
    nornir::MonitoredSample s;
    // First check when saturated
    s.loadPercentage = MAX_RHO;
    s.bandwidth = 100;
    EXPECT_EQ(s.getMaximumBandwidth(), s.bandwidth);

    // Then check when not saturated
    EXPECT_GT(MAX_RHO, 1);
    s.loadPercentage = MAX_RHO - 1;
    s.bandwidth = 100;
    EXPECT_EQ(s.getMaximumBandwidth(), s.bandwidth / (s.loadPercentage / 100.0));
}

TEST(SamplesTest, Load) {
    nornir::MonitoredSample sample;
    std::string fieldStr = "[Watts: 99.9 Knarr Sample: [Load: 90 Bandwidth: 100 "
                           "Latency: 200 NumTasks: 300 "
                           "CustomField0: 0 CustomField1: 1 CustomField2: 2 "
                           "CustomField3: 3 CustomField4: 4 CustomField5: 5 "
                           "CustomField6: 6 CustomField7: 7 CustomField8: 8 "
                           "CustomField9: 9 ] ]";
    std::stringstream ss(fieldStr);
    ss >> sample;
    assert(sample.watts == 99.9);
    assert(sample.loadPercentage == 90);
    assert(sample.bandwidth == 100);
    assert(sample.latency == 200);
    assert(sample.numTasks == 300);
    for(size_t i = 0; i < KNARR_MAX_CUSTOM_FIELDS; i++){
        assert(sample.customFields[i] == i);
    }
}

TEST(SamplesTest, Operators) {
    srand(time(NULL));

    for(unsigned int i = 0; i < 10; i++){
        nornir::MonitoredSample sample;
        sample.bandwidth = rand() % 100;
        sample.latency = rand() % 100 + 1;
        sample.loadPercentage = rand() % 100 + 2;
        sample.numTasks = rand() % 100 + 3;
        for(size_t i = 0; i < KNARR_MAX_CUSTOM_FIELDS; i++){
            sample.customFields[i] = rand() % 100 + 4 + i + 1;
        }
        sample.watts = rand() % 100 + 10;

        // Assignment
        nornir::MonitoredSample sample2 = sample;

        // Multiplication by constant
        sample2 *= 10;
        EXPECT_EQ(sample2.bandwidth, sample.bandwidth*10);
        EXPECT_EQ(sample2.latency, sample.latency*10);
        EXPECT_EQ(sample2.loadPercentage, sample.loadPercentage*10);
        EXPECT_EQ(sample2.numTasks, sample.numTasks*10);
        for(size_t i = 0; i < KNARR_MAX_CUSTOM_FIELDS; i++){
            EXPECT_EQ(sample2.customFields[i], sample.customFields[i]*10);
        }
        EXPECT_EQ(sample2.watts, sample.watts*10);

        // Division by constant
        sample2 /= 10;
        EXPECT_EQ(sample2.bandwidth, sample.bandwidth);
        EXPECT_EQ(sample2.latency, sample.latency);
        EXPECT_EQ(sample2.loadPercentage, sample.loadPercentage);
        EXPECT_EQ(sample2.numTasks, sample.numTasks);
        for(size_t i = 0; i < KNARR_MAX_CUSTOM_FIELDS; i++){
            EXPECT_EQ(sample2.customFields[i], sample.customFields[i]);
        }
        EXPECT_EQ(sample2.watts, sample.watts);

        // Copy constructor and sum
        nornir::MonitoredSample r(sample + sample2);
        EXPECT_EQ(r.bandwidth, sample.bandwidth + sample2.bandwidth);
        EXPECT_EQ(r.latency, sample.latency + sample2.latency);
        EXPECT_EQ(r.loadPercentage, sample.loadPercentage + sample2.loadPercentage);
        EXPECT_EQ(r.numTasks, sample.numTasks + sample2.numTasks);
        for(size_t i = 0; i < KNARR_MAX_CUSTOM_FIELDS; i++){
            EXPECT_EQ(r.customFields[i], sample.customFields[i] + sample2.customFields[i]);
        }
        EXPECT_EQ(r.watts, sample.watts + sample2.watts);

        // Subtraction
        r = sample - sample2;
        EXPECT_EQ(r.bandwidth, 0);
        EXPECT_EQ(r.latency, 0);
        EXPECT_EQ(r.loadPercentage, 0);
        EXPECT_EQ(r.numTasks, 0);
        for(size_t i = 0; i < KNARR_MAX_CUSTOM_FIELDS; i++){
            EXPECT_EQ(r.customFields[i], 0);
        }
        EXPECT_EQ(r.watts, 0);

        // Multiplication
        r = sample * sample2;
        EXPECT_EQ(r.bandwidth, sample.bandwidth * sample2.bandwidth);
        EXPECT_EQ(r.latency, sample.latency * sample2.latency);
        EXPECT_EQ(r.loadPercentage, sample.loadPercentage * sample2.loadPercentage);
        EXPECT_EQ(r.numTasks, sample.numTasks * sample2.numTasks);
        for(size_t i = 0; i < KNARR_MAX_CUSTOM_FIELDS; i++){
            EXPECT_EQ(r.customFields[i], sample.customFields[i] * sample2.customFields[i]);
        }
        EXPECT_EQ(r.watts, sample.watts * sample2.watts);

        // Division
        // Adjust sample2 to avoid zeros before dividing.
        regularize(sample2);
        sample2.bandwidth += 1;
        sample2.latency += 1;
        sample2.loadPercentage += 1;
        sample2.numTasks += 1;
        for(size_t i = 0; i < KNARR_MAX_CUSTOM_FIELDS; i++){
            sample2.customFields[i] += 1;
        }
        sample2.watts += 1;
        r = sample / sample2;
        EXPECT_EQ(r.bandwidth, sample.bandwidth / sample2.bandwidth);
        EXPECT_EQ(r.latency, sample.latency / sample2.latency);
        EXPECT_EQ(r.loadPercentage, sample.loadPercentage / sample2.loadPercentage);
        EXPECT_EQ(r.numTasks, sample.numTasks / sample2.numTasks);
        for(size_t i = 0; i < KNARR_MAX_CUSTOM_FIELDS; i++){
            EXPECT_EQ(r.customFields[i], sample.customFields[i] / sample2.customFields[i]);
        }
        EXPECT_EQ(r.watts, sample.watts / sample2.watts);

        // Sqrt
        r = squareRoot(sample);
        EXPECT_EQ(r.bandwidth, sqrt(sample.bandwidth));
        EXPECT_EQ(r.latency, sqrt(sample.latency));
        EXPECT_EQ(r.loadPercentage, sqrt(sample.loadPercentage));
        EXPECT_EQ(r.numTasks, sqrt(sample.numTasks));
        for(size_t i = 0; i < KNARR_MAX_CUSTOM_FIELDS; i++){
            EXPECT_EQ(r.customFields[i], sqrt(sample.customFields[i]));
        }
        EXPECT_EQ(r.watts, sqrt(sample.watts));

        // Zero
        zero(r);
        EXPECT_EQ(r.bandwidth, 0);
        EXPECT_EQ(r.latency, 0);
        EXPECT_EQ(r.loadPercentage, 0);
        EXPECT_EQ(r.numTasks, 0);
        for(size_t i = 0; i < KNARR_MAX_CUSTOM_FIELDS; i++){
            EXPECT_EQ(r.customFields[i], 0);
        }
        EXPECT_EQ(r.watts, 0);

        // Sqrt on zero
        r = squareRoot(r);
        zero(r);
        EXPECT_EQ(r.bandwidth, 0);
        EXPECT_EQ(r.latency, 0);
        EXPECT_EQ(r.loadPercentage, 0);
        EXPECT_EQ(r.numTasks, 0);
        for(size_t i = 0; i < KNARR_MAX_CUSTOM_FIELDS; i++){
            EXPECT_EQ(r.customFields[i], 0);
        }
        EXPECT_EQ(r.watts, 0);

        // Regularize
        r = r - sample;
        // All the values in r should be < 0 so by regularizing them they should
        // become 0
        regularize(r);
        EXPECT_EQ(r.bandwidth, 0);
        EXPECT_EQ(r.latency, 0);
        EXPECT_EQ(r.loadPercentage, 0);
        EXPECT_EQ(r.numTasks, 0);
        for(size_t i = 0; i < KNARR_MAX_CUSTOM_FIELDS; i++){
            EXPECT_EQ(r.customFields[i], 0);
        }
        EXPECT_EQ(r.watts, 0);

        // Minimum
        // Be sure that sample2 is different from sample.
        sample2 = sample + sample;
        r = minimum(sample, sample2);
        EXPECT_EQ(r.bandwidth, std::min(sample.bandwidth, sample2.bandwidth));
        EXPECT_EQ(r.latency, std::min(sample.latency, sample2.latency));
        EXPECT_EQ(r.loadPercentage, std::min(sample.loadPercentage, sample2.loadPercentage));
        EXPECT_EQ(r.numTasks, std::min(sample.numTasks, sample2.numTasks));
        for(size_t i = 0; i < KNARR_MAX_CUSTOM_FIELDS; i++){
            EXPECT_EQ(r.customFields[i], std::min(sample.customFields[i], sample2.customFields[i]));
        }
        EXPECT_EQ(r.watts, std::min(sample.watts, sample2.watts));

        // Maximum
        r = maximum(sample, sample2);
        EXPECT_EQ(r.bandwidth, std::max(sample.bandwidth, sample2.bandwidth));
        EXPECT_EQ(r.latency, std::max(sample.latency, sample2.latency));
        EXPECT_EQ(r.loadPercentage, std::max(sample.loadPercentage, sample2.loadPercentage));
        EXPECT_EQ(r.numTasks, std::max(sample.numTasks, sample2.numTasks));
        for(size_t i = 0; i < KNARR_MAX_CUSTOM_FIELDS; i++){
            EXPECT_EQ(r.customFields[i], std::max(sample.customFields[i], sample2.customFields[i]));
        }
        EXPECT_EQ(r.watts, std::max(sample.watts, sample2.watts));
    }
}
