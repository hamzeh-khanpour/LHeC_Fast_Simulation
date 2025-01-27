#!/bin/bash

# Navigate to the output directory
cd aa_ww_xs_NP || exit

# Define the FT4 values to scan
FT4_VALUES=(-1.0e-08 -8.0e-09 -6.0e-09 -4.0e-09 -2.0e-09 0.0 2.0e-09 4.0e-09 6.0e-09 8.0e-09 1.0e-08)

# Define the output file for cross-sections
CROSS_SECTION_FILE="cross_sections_FT4.txt"
echo "FT4 Value     Cross-section (pb)" > $CROSS_SECTION_FILE

# Loop over the FT4 values
for FT4 in "${FT4_VALUES[@]}"; do
  # Format FT4 for readable directory naming
  formatted_FT4=$(printf "%.1e" $FT4)

  echo "Running simulation for FT4 = $FT4"

  # Update the FT4 value in Block anoinputs
  sed -i '/Block anoinputs/,/^$/s/^\( *15 \).*$/\1 '"$FT4"' # FT4/' Cards/param_card.dat

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
  echo "Error: Log file $run_log not found for FT4 = $FT4"
fi


  # Save the FT4 value and cross-section to the output file
  echo "$FT4       $cross_section +- $cross_error" >> $CROSS_SECTION_FILE

  # Rename the results directory to include the FT4 value
  run_dir="Events/run_FT4_${formatted_FT4}"
  if [ -d "Events/run_01" ]; then
    mv Events/run_01 $run_dir
    echo "Results saved to $run_dir"
  else
    echo "No results directory found for FT4 = $FT4"
  fi
done

# Print completion message
echo "All simulations completed. Cross-sections saved to $CROSS_SECTION_FILE"
