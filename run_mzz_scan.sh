#!/usr/bin/env bash

LHE="/home/hamzeh-khanpour/MG5_aMC_v3_6_6/AA_ZZ_madspin/Events/run_02_decayed_1/decayed_events.lhe"
SCRIPT="plot_dimuon_pt_mZZ_EFTchecks_from_lhe.py"

for CUT in 1000 1250 1500 1750 2000; do
  OUTDIR="scan_mzz_${CUT}GeV"
  echo "Running mZZ cut = ${CUT} GeV ..."
  python3.10 "$SCRIPT" "$LHE" \
    --mzz-cut "$CUT" \
    --truncate-on mzz \
    --outdir "$OUTDIR"
done

echo
printf "%-12s %-16s %-20s\n" "mZZcut[GeV]" "pass_fraction" "sigma_trunc[pb]"
printf "%-12s %-16s %-20s\n" "-----------" "-------------" "-------------------"

for CUT in 1000 1250 1500 1750 2000; do
  SUM="scan_mzz_${CUT}GeV/eft_truncation_summary.txt"

  FRAC=$(awk -F'= ' '/Passing events:/ {print $2}' "$SUM")
  SIG=$(awk -F': ' '/Estimated truncated cross section:/ {gsub(/ pb/,"",$2); print $2}' "$SUM")

  printf "%-12s %-16s %-20s\n" "$CUT" "$FRAC" "$SIG"
done
