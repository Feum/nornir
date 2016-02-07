#!/bin/bash
#
# Fields
#
# 1: AvgSecondaryLoss
# 2: StddevSecondaryLoss
# 3: AvgCalibrationSteps
# 4: StddevCalibrationSteps
# 5: AvgCalibrationTime
# 6: StddevCalibrationTime
# 7: AvgCalibrationTimePerc
# 8: StddevCalibrationTimePerc
# 9: AvgCalibrationTasks
# 10: StddevCalibrationTasks
# 11: AvgCalibrationTasksPerc
# 12: StddevCalibrationTasksPerc
# 13: AvgReconfigurationTime
# 14: StddevReconfigurationTime

for bench in simple_mandelbrot canneal pbzip2 blackscholes #videoprocessing
do
    echo "======================================" $bench "======================================"
    for run in POWER_BUDGET_REGRESSION_LINEAR_HALTON POWER_BUDGET_REGRESSION_LINEAR_HALTON_FAST PERF_COMPLETION_TIME_LIMARTINEZ PERF_COMPLETION_TIME_REGRESSION_LINEAR_HALTON PERF_COMPLETION_TIME_REGRESSION_LINEAR_HALTON_FAST
    do
        echo "=============================" $run "================================="
        tail -n 2 $bench/$run/results.csv | cut -f 1,7,13 | column -t
    done
    echo "==================================================================================="
    echo 
    echo
done
