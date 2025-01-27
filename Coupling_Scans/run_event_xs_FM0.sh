#!/bin/bash

# Navigate to the output directory
cd aa_ww_xs_NP || exit

# Define the FM0 values to scan
FM0_VALUES=(-1.0e-08 -8.0e-09 -6.0e-09 -4.0e-09 -2.0e-09 0.0 2.0e-09 4.0e-09 6.0e-09 8.0e-09 1.0e-08)

# Define the output file for cross-sections
CROSS_SECTION_FILE="cross_sections.txt"
echo "FM0 Value     Cross-section (pb)" > $CROSS_SECTION_FILE

# Loop over the FM0 values
for FM0 in "${FM0_VALUES[@]}"; do
  # Format FM0 for readable directory naming
  formatted_FM0=$(printf "%.1e" $FM0)

  echo "Running simulation for FM0 = $FM0"

  # Update the FM0 value in Block anoinputs
  sed -i '/Block anoinputs/,/^$/s/^\( *3 \).*$/\1 '"$FM0"' # FM0/' Cards/param_card.dat

  # Launch MadGraph and generate events
  ./bin/generate_events <<EOF
0
EOF

# Extract cross-section from the log file
run_log=$(ls Events/run_*/run_*_tag_1_banner.txt | tail -n 1)
if [ -f "$run_log" ]; then
  # Extract cross-section value and handle variations in spacing
  cross_section=$(grep "Integrated weight (pb)" "$run_log" | awk -F':' '{print $2}' | xargs)
  if [ -z "$cross_section" ]; then
    cross_section="N/A"
  fi
else
  cross_section="N/A"
  echo "Error: Log file $run_log not found for FM0 = $FM0"
fi


  # Save the FM0 value and cross-section to the output file
  echo "$FM0       $cross_section +- $cross_error" >> $CROSS_SECTION_FILE

  # Rename the results directory to include the FM0 value
  run_dir="Events/run_FM0_${formatted_FM0}"
  if [ -d "Events/run_01" ]; then
    mv Events/run_01 $run_dir
    echo "Results saved to $run_dir"
  else
    echo "No results directory found for FM0 = $FM0"
  fi
done

# Print completion message
echo "All simulations completed. Cross-sections saved to $CROSS_SECTION_FILE"
