#!/usr/bin/env bash
set -euo pipefail

min_model=2047

for file in IniFile/IC_Fornax*.ini; do
    [ -f "$file" ] || continue
    name=${file##*/}
    number=${name#IC_Fornax}
    number=${number%.ini}
    [[ "$number" =~ ^[0-9]+$ ]] || continue

    if (( 10#$number >= min_model )); then
        git add -f -- "$file"
    fi
done

echo "Staged Fornax${min_model}+ model files. Review with: git diff --cached --stat"
