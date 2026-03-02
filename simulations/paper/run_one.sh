#!/bin/bash

# (Make sure these are correct!)
STATUS_COMMAND="ace status -j"
COMPLETED_TEXT="SUCCEEDED"
FAILED_TEXT="FAILED"
POLLING_INTERVAL=60
# --- End of Configuration ---

# 1. Get arguments
R_SCRIPT_PATH=$1
if [ -z "$R_SCRIPT_PATH" ]; then
    echo "Error: No R script path provided to run_one.sh. Aborting."
    exit 1
fi

# 2. Get password from environment
if [ -z "$ACE_PASSWORD" ]; then
    echo "Error: ACE_PASSWORD environment variable is not set. Aborting."
    exit 1
fi

# 3. Set up job variables
script_command="Rscript ${R_SCRIPT_PATH}"
description="Simulation job for ${R_SCRIPT_PATH}"
# Create a unique temp file name based on the R script name
safe_name=$(echo ${R_SCRIPT_PATH} | tr '/' '_')
temp_yaml="current_job_${safe_name}.yaml"

echo "--- [Worker] Starting job: ${R_SCRIPT_PATH} ---"

# 4. Create the job.yaml from the template
sed -e "s|PLACEHOLDER_COMMAND|${script_command}|" \
    -e "s|PLACEHOLDER_DESCRIPTION|${description}|" \
    job_template.yaml > ${temp_yaml}

if [ $? -ne 0 ]; then
  echo "Error: Failed to create ${temp_yaml}."
  exit 1
fi

# 5. Submit the job
submit_output=$(ace submit-yaml -p "$ACE_PASSWORD" -f ${temp_yaml})
if [ $? -ne 0 ]; then
  echo "Error: 'ace submit-yaml' failed for ${R_SCRIPT_PATH}."
  rm ${temp_yaml}
  exit 1
fi

# 6. Extract the Job ID (using the robust 'sed' method)
JOB_ID=$(echo ${submit_output} | sed -n "s/.*'Job ID': '\([^']*\)'.*/\1/p")
if [ -z "${JOB_ID}" ] || [ "${JOB_ID}" == "null" ]; then
    echo "Error: Could not parse Job ID for ${R_SCRIPT_PATH}. Output was:"
    echo ${submit_output}
    rm ${temp_yaml}
    exit 1
fi

echo "[Worker] Job ${R_SCRIPT_PATH} submitted. Job ID: ${JOB_ID}. Polling..."

# 7. Polling loop: Wait for the job to finish
while true; do
    # Get the current status text
    current_status=$(${STATUS_COMMAND} ${JOB_ID} | jq -r '.Status')

    # Check if status is COMPLETED
    if [[ "${current_status}" == "${COMPLETED_TEXT}" ]]; then
        echo "--- [Worker] Job ${R_SCRIPT_PATH} (ID: ${JOB_ID}) COMPLETED. ---"
        break # Exit the polling loop

    # Check if status is FAILED
    elif [[ "${current_status}" == "${FAILED_TEXT}" ]]; then
        echo "--- [Worker] Error: Job ${R_SCRIPT_PATH} (ID: ${JOB_ID}) FAILED. ---"
        rm ${temp_yaml}
        exit 1 # Exit this worker script with an error

    # Otherwise, wait and check again
    else
        # This line is removed to avoid spamming the console
        # echo "Status: ${current_status} (waiting...)"
        sleep ${POLLING_INTERVAL}
    fi
done # End of polling loop

# 8. Clean up
rm ${temp_yaml}
exit 0
