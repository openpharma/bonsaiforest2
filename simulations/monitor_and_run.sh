#!/bin/bash

# Monitor running ACE jobs and run next batch after they complete

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

echo ""
echo "============================================="
echo "ACE Job Monitor & Batch Queue"
echo "============================================="
echo ""

# --- Configuration: Job IDs to monitor ---
JOB_1="d9debaf5-bce5-4068-8c89-9cf40b5acd96"
JOB_2="1e7e9b00-816d-4812-b50b-cf9007a3f590"

# Track how many times we've checked
check_count=0

# Monitor until both jobs are no longer RUNNING
while true; do
    check_count=$((check_count + 1))
    
    # Get status of both jobs
    status_1=$(ace status -j "$JOB_1" 2>/dev/null | jq -r '.Status // "UNKNOWN"')
    status_2=$(ace status -j "$JOB_2" 2>/dev/null | jq -r '.Status // "UNKNOWN"')
    
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] Check #$check_count - Job 1: $status_1 | Job 2: $status_2"
    
    # Check if both jobs are done (either SUCCEEDED or FAILED)
    if [[ "$status_1" != "RUNNING" ]] && [[ "$status_2" != "RUNNING" ]]; then
        echo ""
        echo "✓ Both jobs have completed!"
        echo "  Job 1 ($JOB_1): $status_1"
        echo "  Job 2 ($JOB_2): $status_2"
        break
    fi
    
    # Wait 30 seconds before checking again
    sleep 30
done

echo ""
echo "============================================="
echo "Starting next batch of simulations..."
echo "============================================="
echo ""

# Define the next batch of jobs
next_batch=(
  "TTE/RHS_theta0_1_s_2_scenarios_1_3.R"
  "TTE/RHS_theta0_1_s_2_scenarios_4_6.R"
  "TTE/OVAT_1_oneway.R"
  "TTE/OVAT_2_oneway.R"
  "TTE/OVAT_3_oneway.R"
  "TTE/OVAT_4_oneway.R"
  "TTE/OVAT_5_oneway.R"
  "TTE/OVAT_6_oneway.R"
)

# Generate and run jobs 2 at a time
run_jobs() {
  for job in "${next_batch[@]}"; do
    echo "${job}"
  done
}

echo "Submitting ${#next_batch[@]} jobs (2 in parallel)..."
echo ""

run_jobs | xargs -P 2 -I {} bash -c './run_one.sh "{}"'

echo ""
echo "============================================="
echo "All batch jobs completed!"
echo "============================================="
