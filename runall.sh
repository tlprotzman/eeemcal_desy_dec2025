#!/bin/bash

# Trap Ctrl-C (SIGINT) and cleanup
trap 'echo -e "\n\033[31mInterrupted! Stopping all processes...\033[0m"; pkill -P $$; exit 130' INT TERM

echo -e "\033[32mStarting full calibration and analysis chain...\033[0m"
mkdir -p output
rm output/gain_factors.root # remove the old calibrations
rm output/adc_to_gev_calibration.root
rm output/tot_widths.root
rm output/tot_calibration_values.root

set -e

echo -e "\033[32mRunning gain_match.C\033[0m"
root -q -b -x -l gain_match.cxx
echo -e "\033[32mRunning adc_calibration.C\033[0m"
root -q -b -x -l adc_calibration.cxx
echo -e "\033[32mRunning tot_calibration.C\033[0m"
root -q -b -x -l tot_calibration.cxx
echo -e "\033[32mRunning tot_calibration.C\033[0m"
root -q -b -x -l tot_calibration.cxx # Needs to run twice to generate the cuts
echo -e "\033[32mRunning energy_resolution.C\033[0m"
root -q -b -x -l energy_resolution.cxx
