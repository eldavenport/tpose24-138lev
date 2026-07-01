#!/usr/bin/env bash
#
# mitgcm_max_uvw.sh — scan MITgcm STDOUT log(s) from the monitor package and
# report, for u/v/w/eta and for the advective CFL number in each direction:
#   - the max value (and the timestep/file it occurred at)
#   - the min value (and the timestep/file it occurred at)
#   - the mean value
# Also prints a "top CFL_w events vs. eta" table to help spot whether large
# vertical-CFL events line up with free-surface (eta) excursions.
#
# velocity/eta (dynstat) block:
#   Uses dynstat_{uvel,vvel,wvel,eta}_max/min/mean (requires monitorSelect
#   >= 2 in "data" so per-field statistics get printed; note dynstat_eta_*
#   only needs monitorSelect >= 1). max/min are the true domain extrema
#   (over the whole run); mean is the time-average of the per-timestep
#   domain-mean value. eta is sea-surface height in meters, not a velocity.
#
# CFL block:
#   MITgcm's monitor only ever reports the domain MAX advective CFL each
#   timestep (advcfl_uvel_max, advcfl_vvel_max, advcfl_wvel_max,
#   advcfl_W_hf_max) — there is no spatial min/mean CFL in STDOUT. So here
#   "max/min/mean" are computed across the time series of that per-step
#   max value, i.e.:
#     max  = worst-case CFL anywhere/anytime in the run
#     min  = the smallest the per-step worst-case CFL ever got
#     mean = the time-average of the per-step worst-case CFL
#   advcfl_W_hf_max is the partial-cell/hFac-corrected vertical CFL (see
#   earlier discussion) — it will differ from advcfl_wvel_max most near
#   topography, but also near the surface if using a nonlinear free
#   surface (since hFacC at the top layer varies with eta).
#
# Top CFL_w-vs-eta table:
#   For each timestep, remembers the most recently reported dynstat_eta_max
#   /dynstat_eta_min, then whenever advcfl_wvel_max / advcfl_W_hf_max is
#   read, records (tstep, cfl value, eta_max, eta_min) at that same
#   timestep. Prints the top N timesteps (by cfl value) so you can eyeball
#   whether the worst vertical-CFL events coincide with eta extremes
#   (pointing at the free surface / thin top cell) or not (pointing
#   elsewhere, e.g. a numerical instability in the momentum advection
#   scheme rather than a free-surface artifact).
#
# Usage:
#   ./mitgcm_max_uvw.sh [-n TOPN] [STDOUT_file_or_glob ...]
#
# Examples:
#   ./mitgcm_max_uvw.sh                     # uses STDOUT.* in cwd
#   ./mitgcm_max_uvw.sh STDOUT.0000
#   ./mitgcm_max_uvw.sh -n 20 STDOUT.0000
#   ./mitgcm_max_uvw.sh run1/STDOUT.* run2/STDOUT.*

set -euo pipefail

topn=10
if [ "${1:-}" = "-n" ]; then
    topn="$2"
    shift 2
fi

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

awk -v TOPK="$topn" '
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

    # keep the top TOPK (ts, val, etaMax, etaMin) entries for a given
    # prefix, sorted descending by val (simple bounded insertion sort)
    function insert_top(pfx, tsVal, val, eMax, eMin,   idx, cnt, t) {
        cnt = topCnt[pfx] + 0
        idx = cnt + 1
        topTs[pfx, idx] = tsVal
        topVal[pfx, idx] = val + 0
        topEMax[pfx, idx] = eMax
        topEMin[pfx, idx] = eMin
        cnt++
        while (idx > 1 && topVal[pfx, idx-1] < topVal[pfx, idx]) {
            t = topTs[pfx, idx-1];   topTs[pfx, idx-1]   = topTs[pfx, idx];   topTs[pfx, idx]   = t
            t = topVal[pfx, idx-1];  topVal[pfx, idx-1]  = topVal[pfx, idx];  topVal[pfx, idx]  = t
            t = topEMax[pfx, idx-1]; topEMax[pfx, idx-1] = topEMax[pfx, idx]; topEMax[pfx, idx] = t
            t = topEMin[pfx, idx-1]; topEMin[pfx, idx-1] = topEMin[pfx, idx]; topEMin[pfx, idx] = t
            idx--
        }
        if (cnt > TOPK) {
            delete topTs[pfx, cnt]; delete topVal[pfx, cnt]
            delete topEMax[pfx, cnt]; delete topEMin[pfx, cnt]
            cnt = TOPK
        }
        topCnt[pfx] = cnt
    }

    /%MON dynstat_uvel_max/  { upd("u", $NF);   haveDynstat = 1 }
    /%MON dynstat_uvel_min/  { upd("u", $NF);   haveDynstat = 1 }
    /%MON dynstat_uvel_mean/ { accum("u", $NF); haveDynstat = 1 }
    /%MON dynstat_vvel_max/  { upd("v", $NF);   haveDynstat = 1 }
    /%MON dynstat_vvel_min/  { upd("v", $NF);   haveDynstat = 1 }
    /%MON dynstat_vvel_mean/ { accum("v", $NF); haveDynstat = 1 }
    /%MON dynstat_wvel_max/  { upd("w", $NF);   haveDynstat = 1 }
    /%MON dynstat_wvel_min/  { upd("w", $NF);   haveDynstat = 1 }
    /%MON dynstat_wvel_mean/ { accum("w", $NF); haveDynstat = 1 }
    /%MON dynstat_eta_max/   { upd("eta", $NF);   haveDynstat = 1; curEtaMax = $NF }
    /%MON dynstat_eta_min/   { upd("eta", $NF);   haveDynstat = 1; curEtaMin = $NF }
    /%MON dynstat_eta_mean/  { accum("eta", $NF); haveDynstat = 1 }

    # advective CFL: only a per-step max is ever printed by MITgcm, so
    # max/min/mean here are computed across the time series of that value
    /%MON advcfl_uvel_max/  { upd("cfl_u", $NF);   accum("cfl_u", $NF);   haveCfl = 1 }
    /%MON advcfl_vvel_max/  { upd("cfl_v", $NF);   accum("cfl_v", $NF);   haveCfl = 1 }
    /%MON advcfl_wvel_max/  {
        upd("cfl_w", $NF); accum("cfl_w", $NF); haveCfl = 1
        insert_top("wvel", ts, $NF, curEtaMax, curEtaMin)
    }
    /%MON advcfl_W_hf_max/  {
        upd("cfl_whf", $NF); accum("cfl_whf", $NF); haveCfl = 1; haveWhf = 1
        insert_top("whf", ts, $NF, curEtaMax, curEtaMin)
    }

    END {
        if (haveDynstat) {
            print "=== velocity / eta (dynstat) ==="
            print "(u,v,w in m/s; eta = sea-surface height in m)"
            printf "%-3s %16s  %-10s %-24s   %16s  %-10s %-24s   %16s\n", \
                "var", "max", "@ tstep", "(file)", "min", "@ tstep", "(file)", "mean"
            n = split("u v w eta", comps, " ")
            for (i = 1; i <= n; i++) {
                c = comps[i]
                if (c in seen) {
                    meanStr = (meanCnt[c] > 0) ? sprintf("%16.6e", meanSum[c] / meanCnt[c]) : sprintf("%16s", "NA")
                    printf "%-3s %16.6e  %-10s %-24s   %16.6e  %-10s %-24s   %s\n", \
                        c, umax[c], umaxts[c], umaxfile[c], umin[c], umints[c], uminfile[c], meanStr
                }
            }
        } else {
            print "No dynstat_{u,v,w,eta} lines found (per-field stats not enabled)."
            print "Set monitorSelect >= 1 (eta) or >= 2 (u/v/w) in \"data\"."
        }

        print ""

        if (haveCfl) {
            print "=== advective CFL (per-timestep domain max; min/mean are over time, not space) ==="
            printf "%-8s %16s  %-10s %-24s   %16s  %-10s %-24s   %16s\n", \
                "var", "max", "@ tstep", "(file)", "min", "@ tstep", "(file)", "mean"
            n = split("cfl_u cfl_v cfl_w cfl_whf", comps, " ")
            for (i = 1; i <= n; i++) {
                c = comps[i]
                if (c in seen) {
                    meanStr = (meanCnt[c] > 0) ? sprintf("%16.6e", meanSum[c] / meanCnt[c]) : sprintf("%16s", "NA")
                    printf "%-8s %16.6e  %-10s %-24s   %16.6e  %-10s %-24s   %s\n", \
                        c, umax[c], umaxts[c], umaxfile[c], umin[c], umints[c], uminfile[c], meanStr
                }
            }
        } else {
            print "No advcfl_*vel_max lines found. Check that pkg/monitor is enabled and monitorFreq is set in data."
        }

        print ""

        # top CFL_w / CFL_W_hf events alongside eta at the same timestep
        pfxList = "wvel"
        if (haveWhf) pfxList = pfxList " whf"
        nPfx = split(pfxList, pfxArr, " ")
        for (p = 1; p <= nPfx; p++) {
            pfx = pfxArr[p]
            label = (pfx == "wvel") ? "advcfl_wvel_max" : "advcfl_W_hf_max"
            cnt = topCnt[pfx] + 0
            if (cnt == 0) continue
            print "=== top " cnt " " label " events vs. eta at same timestep ==="
            printf "%-10s %16s   %16s %16s\n", "@ tstep", label, "eta_max(t)", "eta_min(t)"
            for (i = 1; i <= cnt; i++) {
                printf "%-10s %16.6e   %16.6e %16.6e\n", \
                    topTs[pfx, i], topVal[pfx, i], topEMax[pfx, i]+0, topEMin[pfx, i]+0
            }
            print ""
        }

        if (!haveDynstat && !haveCfl) exit 1
    }
' "${existing[@]}"
