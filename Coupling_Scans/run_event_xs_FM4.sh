#!/bin/bash

# Navigate to the output directory
cd aa_ww_xs_NP || exit

# Define the FM4 values to scan
FM4_VALUES=(-1.0e-08 -8.0e-09 -6.0e-09 -4.0e-09 -2.0e-09 0.0 2.0e-09 4.0e-09 6.0e-09 8.0e-09 1.0e-08)

# Define the output file for cross-sections
CROSS_SECTION_FILE="cross_sections_FM4.txt"
echo "FM4 Value     Cross-section (pb)" > $CROSS_SECTION_FILE

# Loop over the FM4 values
for FM4 in "${FM4_VALUES[@]}"; do
  # Format FM4 for readable directory naming
  formatted_FM4=$(printf "%.1e" $FM4)

  echo "Running simulation for FM4 = $FM4"

  # Update the FM4 value in Block anoinputs
  sed -i '/Block anoinputs/,/^$/s/^\( *7 \).*$/\1 '"$FM4"' # FM4/' Cards/param_card.dat

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
  echo "Error: Log file $run_log not found for FM4 = $FM4"
fi


  # Save the FM4 value and cross-section to the output file
  echo "$FM4       $cross_section +- $cross_error" >> $CROSS_SECTION_FILE

  # Rename the results directory to include the FM4 value
  run_dir="Events/run_FM4_${formatted_FM4}"
  if [ -d "Events/run_01" ]; then
    mv Events/run_01 $run_dir
    echo "Results saved to $run_dir"
  else
    echo "No results directory found for FM4 = $FM4"
  fi
done

# Print completion message
echo "All simulations completed. Cross-sections saved to $CROSS_SECTION_FILE"
