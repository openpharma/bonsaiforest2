#!/bin/bash

# --- Securely get the password once ---
echo "Please enter your 'ace' platform password:"
read -s PASSWORD

if [ -z "$PASSWORD" ]; then
    echo "No password entered. Aborting."
    exit 1
fi

# Export the password so all child scripts (run_one.sh) can see it
export ACE_PASSWORD=$PASSWORD
# Clear the variable from this script's memory for security
unset PASSWORD

echo "Password received. Starting parallel batch (2 jobs at a time)..."
echo ""

# --- Configuration: Define the simulation list ---
data_types=("Count")
models=("Horseshoe_strong" "Horseshoe_mid" "Horseshoe_low" "R2D2_strong" "R2D2_mid" "R2D2_low")

# --- Main Logic: Generate list and pipe to xargs ---

# This function generates the list of all jobs, one per line
generate_job_list() {
  for dtype in "${data_types[@]}"; do
    for model in "${models[@]}"; do

      script_name="${dtype}/${model}.R"
      # Just print the script path to standard output
      echo "${script_name}"

    done
  done
}

# The core of the parallel processing:
# 1. generate_job_list: Prints the full list of R scripts.
# 2. | (pipe): Sends that list as input to xargs.
# 3. xargs -P 3: Reads the list and runs up to 3 commands in parallel.
# 3. -I {}: Takes one line from the list and puts it where {} is.
# 5. bash -c './run_one.sh "{}"': The command to run. We use bash -c
#    to ensure './run_one.sh "TTE/Horseshoe_strong.R"' is called correctly.

generate_job_list | xargs -P 1 -I {} bash -c './run_one.sh "{}"'

echo ""
echo "--------------------------------------------------"
echo "All simulation jobs completed."
