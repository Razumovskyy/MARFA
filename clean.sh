#!/bin/bash

# Find and delete all .mod files in the src directory
# Removes all dependencies (new build is needed after that)
# Find and delete all .par files in the root directory and its subdirectories
# Find and delete all .dat fiels in the root directory only

CURRENT_DIR=$(basename "$PWD")

if [ "$CURRENT_DIR" != "MARFA" ]; then
  echo "Error: This script must be run from the MARFA directory."
  exit 1
fi

find ./src -type f -name "*.mod" -exec rm -f {} +
find ./build -mindepth 1 -exec rm -rf {} +
find . -type f -name "*.par" -exec rm -f {} +
find . -maxdepth 1 -type f -name "*.dat" -exec rm -f {} +

