#!/bin/bash

CURRENT_DIR=$(basename "$PWD")

if [ "$CURRENT_DIR" != "MARFA" ]; then
  echo "Error: This script must be run from the MARFA directory."
  exit 1
fi

# Find and delete all .par files in the root directory and its subdirectories
# Find and delete all .dat fiels in the root directory only
find . -type f -name "*.par" -exec rm -f {} +
find . -maxdepth 1 -type f -name "*.dat" -exec rm -f {} +

echo "All .par files in the root directory and its subdirectories have been deleted."
echo "All .dat files in the root directory are deleted"
