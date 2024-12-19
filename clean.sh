#!/bin/bash

CURRENT_DIR=$(basename "$PWD")

if [ "$CURRENT_DIR" != "MARFA" ]; then
  echo "Error: This script must be run from the MARFA directory."
  exit 1
fi

prompt_user() {
  local action_desc="$1"
  local command="$2"
  
  while true; do
    read -p "Do you want to ${action_desc}? (y/n): " yn
    case $yn in
      [Yy]* ) 
        echo "Executing: ${command}"
        eval "$command"
        break
        ;;
      [Nn]* ) 
        echo "Skipping: ${action_desc}"
        break
        ;;
      * ) 
        echo "Please answer y or n."
        ;;
    esac
  done
}

prompt_user "delete dependencies and .mod files (new build will be required)" \
  'find ./src -type f -name "*.mod" -exec rm -f {} + && find ./build -mindepth 1 -exec rm -rf {} +'

prompt_user "delete all previous PT-tables calculated" \
  'rm -rf ./output/ptTables/*'

prompt_user "delete all human-readable absorption data" \
  'rm -rf ./output/processedData/*'

prompt_user "delete all plots" \
  'rm -rf ./output/plots/*'

prompt_user "remove old .par files from all directories" \
  'find . -type f -name "*.par" -exec rm -f {} +'

prompt_user "remove old .dat files ONLY from the root directory" \
  'find . -maxdepth 1 -type f -name "*.dat" -exec rm -f {} +'

