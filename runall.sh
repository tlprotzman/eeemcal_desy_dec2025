mkdir -p output
rm output/gain_factors.root # remove the old calibrations
rm output/adc_to_gev_calibration.root
rm output/tot_widths.root
rm output/tot_calibration_values.root
root -q -b -x -l gain_match.C
root -q -b -x -l adc_calibration.C
root -q -b -x -l tot_calibration.C
root -q -b -x -l tot_calibration.C # Needs to run twice to generate the cuts
root -q -b -x -l energy_resolution.C