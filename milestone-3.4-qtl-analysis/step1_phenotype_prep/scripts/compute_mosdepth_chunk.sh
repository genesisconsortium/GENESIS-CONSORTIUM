#!/bin/bash
set -euo pipefail

chunk_file="$1"
out_file="$2"

echo -e "sample_id\tmean_autosome_depth\tmean_chrX_depth\tmean_chrY_depth\tn_autosome\tn_chrX\tn_chrY\tchrX_autosome_ratio\tchrY_autosome_ratio" > "$out_file"

awk -F'\t' '{print $1"\t"$6}' "$chunk_file" | while IFS=$'\t' read -r sample_id input_path
do
  if [[ -z "$sample_id" || -z "$input_path" ]]; then
    continue
  fi

  if [[ ! -f "$input_path" ]]; then
    echo "Missing: $sample_id $input_path" >&2
    continue
  fi

  zcat "$input_path" | awk -v sample="$sample_id" '
  BEGIN { OFS="\t" }
  {
    chr=$1
    depth=$5

    if (chr == "X" || chr == "chrX") group="chrX"
    else if (chr == "Y" || chr == "chrY") group="chrY"
    else if (chr ~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/) group="autosome"
    else next

    sum[group] += depth
    n[group]++
  }
  END {
    auto = (n["autosome"] > 0 ? sum["autosome"] / n["autosome"] : "NA")
    x    = (n["chrX"] > 0 ? sum["chrX"] / n["chrX"] : "NA")
    y    = (n["chrY"] > 0 ? sum["chrY"] / n["chrY"] : "NA")

    xr = (auto != "NA" && x != "NA" && auto > 0 ? x / auto : "NA")
    yr = (auto != "NA" && y != "NA" && auto > 0 ? y / auto : "NA")

    print sample, auto, x, y, n["autosome"]+0, n["chrX"]+0, n["chrY"]+0, xr, yr
  }' >> "$out_file"
done
