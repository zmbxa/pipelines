#!/bin/bash
### This script is for auto-detect all Report file in txt format in current directory and sub-directories. 
### It will generate a summary table to make statistics and visualization easier.
### Usage:
###	bash /storage/zhangyanxiaoLab/niuyuxiao/pipelines/Paired-Tag-Pipeline/pipelines/ReportSummary.sh [path/to/workingDir]
###  If not specifiied, it will work on the current directory.
###

help() {
    sed -rn 's/^### ?//;T;p;' "$0"
}

if [[ "$1" == "-h" ]]|| [[ "$1" == "--help" ]]; then
  help
  exit 1
fi

if [[ -z $1 ]]; then 
  folder_path=$(pwd)
else
   folder_path=$1
fi

output_file="reportSummary.csv"


indicators=("Raw" "Fully ligated" "Useful reads" "Duplication" "library")

file_paths=()

# add *Reports.txt path to array
find_files() {
    for file in "$1"/*; do
        if [[ -d "$file" ]]; then
            find_files "$file"
        elif [[ "$file" == *"Report.txt" ]]; then
            file_paths+=("$file")
        fi
    done
}

# find txt file
find_files "$folder_path"

# save values
declare -A values

# read content for each file
for file_path in "${file_paths[@]}"; do
    while IFS= read -r line; do
        line=$(echo "$line" | tr -d '\r')  # delete enter

        # extract indicators and values
        for indicator in "${indicators[@]}"; do
            if [[ "$line" == *"$indicator"* ]]; then
                value=$(echo "$line" | awk -F ':' '{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}')
                values["$file_path,$indicator"]=$value
            fi
        done
    done < "$file_path"
done

# save to summary file
echo "File,Raw,Fully ligated,Useful reads,Duplication,library,Sequencing Depth(G)" > "$output_file"
for file_path in "${file_paths[@]}"; do
    line="$file_path"
    for indicator in "${indicators[@]}"; do
        line+=","${values["$file_path,$indicator"]}
    done
    raw_value=${values["$file_path,Raw"]}
    raw_divided=$(awk "BEGIN { printf \"%.6f\", ${raw_value% (*}/3300000 }")
    line+=",$raw_divided"
   
    echo "$line" >> "$output_file"
done