#!/bin/bash

# Check if a folder name is provided as an argument
if [ -z "$1" ]; then
    echo "Usage: $0 <folder_name>"
    exit 1
fi

# Create the folder
mkdir -p "$1"

# Check if the folder was created successfully
if [ $? -eq 0 ]; then
    echo "Folder '$1' created successfully."
else
    echo "Failed to create folder '$1'."
    exit 1
fi