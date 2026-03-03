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
data_types=("TTE" "Continuous")
models_TTE=(
  "HN_global_phi1"
  "HN_global_phi_delta_plan"
  "HN_global_phi_delta_plan_half"
  "OVAT_1_oneway_HN_phi_1"
  "OVAT_1_oneway_HN_phi_delta_plan"
  "OVAT_1_oneway_HN_phi_delta_plan_half"
  "OVAT_2_oneway_HN_phi_1"
  "OVAT_2_oneway_HN_phi_delta_plan"
  "OVAT_2_oneway_HN_phi_delta_plan_half"
  "OVAT_3_oneway_HN_phi_1"
  "OVAT_3_oneway_HN_phi_delta_plan"
  "OVAT_3_oneway_HN_phi_delta_plan_half"
  "OVAT_4_oneway_HN_phi_1"
  "OVAT_4_oneway_HN_phi_delta_plan"
  "OVAT_4_oneway_HN_phi_delta_plan_half"
  "RHS_theta0_1_s_2"
  "RHS_theta0_delta_plan_10_s_2"
  "RHS_theta0_delta_plan_half_s_2"
  "RHS_theta0_delta_plan_s_2"
)

models_Continuous=(
  "HN_global_phi_delta_half"
  "HN_global_phi_delta_plan"
  "HN_global_phi_sigma_plan"
  "Horseshoe_global_delta_10"
  "Horseshoe_global_delta_half"
  "Horseshoe_global_delta_plan"
  "Horseshoe_global_sigma_plan"
  "OVAT_oneway_HN_phi_delta_half"
  "OVAT_oneway_HN_phi_delta_plan"
  "OVAT_oneway_HN_phi_sigma_plan"
)

# --- Main Logic: Generate list and pipe to xargs ---

# This function generates the list of all jobs, one per line
generate_job_list() {
  for dtype in "${data_types[@]}"; do
    # Dynamically reference the models array for this data type
    models_var="models_${dtype}"

    # Get the array values using indirect expansion
    eval "local models_list=(\"\${${models_var}[@]}\")"

    for model in "${models_list[@]}"; do
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

generate_job_list | xargs -P 2 -I {} bash -c './run_one.sh "{}"'

echo ""
echo "--------------------------------------------------"
echo "All simulation jobs completed."
