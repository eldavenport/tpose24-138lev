#!/usr/bin/env bash
#
# mitgcm_max_uvw.sh — scan MITgcm STDOUT log(s) from the monitor package and
# report, for u/v/w and for the advective CFL number in each direction:
#   - the max value (and the timestep/file it occurred at)
#   - the min value (and the timestep/file it occurred at)
#   - the mean value
#
# u/v/w block:
#   Uses dynstat_uvel_max/min/mean, dynstat_vvel_max/min/mean,
#   dynstat_wvel_max/min/mean (requires monitorSelect >= 2 in "data" so
#   per-field statistics get printed). max/min are the true domain
#   extrema (over the whole run); mean is the time-average of the
#   per-timestep domain-mean value.
#
# CFL block:
#   MITgcm's monitor only ever reports the domain MAX advective CFL each
#   timestep (advcfl_uvel_max, advcfl_vvel_max, advcfl_wvel_max) — there
#   is no spatial min/mean CFL in STDOUT. So here "max/min/mean" are
#   computed across the time series of that per-step max value, i.e.:
#     max  = worst-case CFL anywhere/anytime in the run
#     min  = the smallest the per-step worst-case CFL ever got
#     mean = the time-average of the per-step worst-case CFL
#
# Usage:
#   ./mitgcm_max_uvw.sh [STDOUT_file_or_glob ...]
#
# Examples:
#   ./mitgcm_max_uvw.sh                     # uses STDOUT.* in cwd
#   ./mitgcm_max_uvw.sh STDOUT.0000
#   ./mitgcm_max_uvw.sh run1/STDOUT.* run2/STDOUT.*

set -euo pipefail

if [ "$#" -gt 0 ]; then
    files=("$@")
else
    files=(STDOUT.*)
fi

# Expand globs safely / fail with a clear message if nothing matches
existing=()
for f in "${files[@]}"; do
    [ -e "$f" ] && existing+=("$f")
done
if [ "${#existing[@]}" -eq 0 ]; then
    echo "No STDOUT files found (looked for: ${files[*]})" >&2
    exit 1
fi

awk '
    /%MON time_tsnumber/ { ts = $NF }

    # track running max/min (with location) and running sum/count (for mean)
    function upd(comp, val,    v) {
        v = val + 0
        if (!(comp in seen)) {
            umax[comp] = v; umaxts[comp] = ts; umaxfile[comp] = FILENAME
            umin[comp] = v; umints[comp] = ts; uminfile[comp] = FILENAME
            seen[comp] = 1
        } else {
            if (v > umax[comp]) { umax[comp] = v; umaxts[comp] = ts; umaxfile[comp] = FILENAME }
            if (v < umin[comp]) { umin[comp] = v; umints[comp] = ts; uminfile[comp] = FILENAME }
        }
    }

    # accumulate for a straight mean of whatever values are fed to it
    function accum(comp, val) {
        meanSum[comp] += (val + 0)
        meanCnt[comp] += 1
    }

    /%MON dynstat_uvel_max/  { upd("u", $NF);              haveDynstat = 1 }
    /%MON dynstat_uvel_min/  { upd("u", $NF);              haveDynstat = 1 }
    /%MON dynstat_uvel_mean/ { accum("u", $NF);            haveDynstat = 1 }
    /%MON dynstat_vvel_max/  { upd("v", $NF);              haveDynstat = 1 }
    /%MON dynstat_vvel_min/  { upd("v", $NF);              haveDynstat = 1 }
    /%MON dynstat_vvel_mean/ { accum("v", $NF);            haveDynstat = 1 }
    /%MON dynstat_wvel_max/  { upd("w", $NF);              haveDynstat = 1 }
    /%MON dynstat_wvel_min/  { upd("w", $NF);              haveDynstat = 1 }
    /%MON dynstat_wvel_mean/ { accum("w", $NF);            haveDynstat = 1 }

    # advective CFL: only a per-step max is ever printed by MITgcm, so
    # max/min/mean here are computed across the time series of that value
    /%MON advcfl_uvel_max/ { upd("cfl_u", $NF); accum("cfl_u", $NF); haveCfl = 1 }
    /%MON advcfl_vvel_max/ { upd("cfl_v", $NF); accum("cfl_v", $NF); haveCfl = 1 }
    /%MON advcfl_wvel_max/ { upd("cfl_w", $NF); accum("cfl_w", $NF); haveCfl = 1 }

    END {
        if (haveDynstat) {
            print "=== velocity (dynstat) ==="
            printf "%-3s %16s  %-10s %-24s   %16s  %-10s %-24s   %16s\n", \
                "var", "max", "@ tstep", "(file)", "min", "@ tstep", "(file)", "mean"
            n = split("u v w", comps, " ")
            for (i = 1; i <= n; i++) {
                c = comps[i]
                if (c in seen) {
                    meanStr = (meanCnt[c] > 0) ? sprintf("%16.6e", meanSum[c] / meanCnt[c]) : sprintf("%16s", "NA")
                    printf "%-3s %16.6e  %-10s %-24s   %16.6e  %-10s %-24s   %s\n", \
                        c, umax[c], umaxts[c], umaxfile[c], umin[c], umints[c], uminfile[c], meanStr
                }
            }
        } else {
            print "No dynstat_{u,v,w}vel lines found (per-field stats not enabled)."
            print "Set monitorSelect >= 2 in \"data\" to get dynstat_*vel_max/min/mean."
        }

        print ""

        if (haveCfl) {
            print "=== advective CFL (per-timestep domain max; min/mean are over time, not space) ==="
            printf "%-6s %16s  %-10s %-24s   %16s  %-10s %-24s   %16s\n", \
                "var", "max", "@ tstep", "(file)", "min", "@ tstep", "(file)", "mean"
            n = split("cfl_u cfl_v cfl_w", comps, " ")
            for (i = 1; i <= n; i++) {
                c = comps[i]
                if (c in seen) {
                    meanStr = (meanCnt[c] > 0) ? sprintf("%16.6e", meanSum[c] / meanCnt[c]) : sprintf("%16s", "NA")
                    printf "%-6s %16.6e  %-10s %-24s   %16.6e  %-10s %-24s   %s\n", \
                        c, umax[c], umaxts[c], umaxfile[c], umin[c], umints[c], uminfile[c], meanStr
                }
            }
        } else {
            print "No advcfl_*vel_max lines found. Check that pkg/monitor is enabled and monitorFreq is set in data."
        }

        if (!haveDynstat && !haveCfl) exit 1
    }
' "${existing[@]}"
