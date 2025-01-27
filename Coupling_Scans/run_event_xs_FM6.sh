#!/bin/bash

# Navigate to the output directory
cd aa_ww_xs_NP || exit

# Define the FM6 values to scan
FM6_VALUES=(-1.0e-08 -8.0e-09 -6.0e-09 -4.0e-09 -2.0e-09 0.0 2.0e-09 4.0e-09 6.0e-09 8.0e-09 1.0e-08)

# Define the output file for cross-sections
CROSS_SECTION_FILE="cross_sections_FM6.txt"
echo "FM6 Value     Cross-section (pb)" > $CROSS_SECTION_FILE

# Loop over the FM6 values
for FM6 in "${FM6_VALUES[@]}"; do
  # Format FM6 for readable directory naming
  formatted_FM6=$(printf "%.1e" $FM6)

  echo "Running simulation for FM6 = $FM6"

  # Update the FM6 value in Block anoinputs
  sed -i '/Block anoinputs/,/^$/s/^\( *9 \).*$/\1 '"$FM6"' # FM6/' Cards/param_card.dat

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
  echo "Error: Log file $run_log not found for FM6 = $FM6"
fi


  # Save the FM6 value and cross-section to the output file
  echo "$FM6       $cross_section +- $cross_error" >> $CROSS_SECTION_FILE

  # Rename the results directory to include the FM6 value
  run_dir="Events/run_FM6_${formatted_FM6}"
  if [ -d "Events/run_01" ]; then
    mv Events/run_01 $run_dir
    echo "Results saved to $run_dir"
  else
    echo "No results directory found for FM6 = $FM6"
  fi
done

# Print completion message
echo "All simulations completed. Cross-sections saved to $CROSS_SECTION_FILE"
