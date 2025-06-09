#!/bin/bash

# Navigate to the output directory
cd aa_ww_xs_NP || { echo "❌ Directory 'aa_ww_xs_NP' not found. Exiting."; exit 1; }

# Define the FM2 values to scan
FM2_VALUES=(
  -1000.0e-12 -800.0e-12 -600.0e-12 -400.0e-12 -200.0e-12
  -100.0e-12 -80.0e-12 -60.0e-12 -40.0e-12 -20.0e-12 -10.0e-12
  0.0
  10.0e-12 20.0e-12 40.0e-12 60.0e-12 80.0e-12 100.0e-12
  200.0e-12 400.0e-12 600.0e-12 800.0e-12 1000.0e-12
)

# Define the output file for cross-sections
CROSS_SECTION_FILE="cross_sections_FM2.txt"
echo -e "FM2 Value (e-12)\tCross-section (pb) ± Uncertainty" > "$CROSS_SECTION_FILE"

# Loop over the FM2 values
for FM2 in "${FM2_VALUES[@]}"; do
  # Format FM2 for readable naming
  formatted_FM2=$(printf "%.1e" "$FM2")

  echo "🔄 Running simulation for FM2 = $FM2"

  # Update the FM2 value in Block anoinputs, line with ID '5'
  sed -i 's/^\( *5 \).*/\1 '"$FM2"' # FM2/' Cards/param_card.dat

  # Run MadGraph simulation
  ./bin/generate_events <<EOF
0
EOF

  # Find the most recent banner log
  run_log=$(ls Events/run_*/run_*_tag_1_banner.txt 2>/dev/null | sort | tail -n 1)

  if [[ -f "$run_log" ]]; then
    # Extract the cross-section and uncertainty (if available)
    cs_line=$(grep "Integrated weight (pb)" "$run_log")
    cross_section=$(echo "$cs_line" | awk -F'[:±+]' '{print $2}' | xargs)
    cross_error=$(echo "$cs_line" | awk -F'[:±+]' '{print $3}' | xargs)

    # Handle missing values
    [[ -z "$cross_section" ]] && cross_section="N/A"
    [[ -z "$cross_error" ]] && cross_error="N/A"
  else
    echo "⚠️  Log file not found for FM2 = $FM2"
    cross_section="N/A"
    cross_error="N/A"
  fi

  # Save the result
  echo -e "$FM2\t$cross_section ± $cross_error" >> "$CROSS_SECTION_FILE"

  # Rename the run directory
  latest_run_dir=$(ls -d Events/run_* 2>/dev/null | sort | tail -n 1)
  if [[ -d "$latest_run_dir" ]]; then
    new_run_dir="Events/run_FM2_${formatted_FM2}"
    mv "$latest_run_dir" "$new_run_dir"
    echo "✅ Results saved to $new_run_dir"
  else
    echo "⚠️  No run directory found to rename for FM2 = $FM2"
  fi
done

# Completion message
echo "🎉 All simulations completed. Cross-sections saved to $CROSS_SECTION_FILE"
