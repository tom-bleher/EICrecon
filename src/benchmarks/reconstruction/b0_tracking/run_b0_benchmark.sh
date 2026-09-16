#!/usr/bin/env bash
# SPDX-License-Identifier: LGPL-3.0-or-later
# One-command B0 tracking benchmark.
#
# Analysis only:
#   run_b0_benchmark.sh --reco reco.edm4eic.root --out results/b0
#
# Reconstruct then analyse (must run inside eic-shell, with EICrecon sourced):
#   run_b0_benchmark.sh --sim sim.edm4hep.root --out results/b0
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_b0_benchmark.sh --reco FILE.edm4eic.root [--counters FILE.json] [--out DIR]
  run_b0_benchmark.sh --sim FILE.edm4hep.root [--out DIR] [--nevents N]
  run_b0_benchmark.sh --self-test

Writes DIR/report.md, DIR/summary.json, DIR/provenance.json, and
DIR/b0_tracking_validation.pdf. Reconstruction options are the B0 defaults;
this script must not change seeder or CKF parameters.

By default reconstruction uses the ACTS material map declared by the loaded
DD4hep detector configuration and enables the strict runtime contract:
  * selected map content must match the geometry-declared map;
  * all B0 approach surfaces must have ACTS material.
Set B0_STRICT_MATERIAL_MAP=0 only for an intentional diagnostic run. Set
ACTS_MATERIAL_MAP only for an intentional map override; in strict mode an old or
different override is rejected. ACTS_MATERIAL_MAP_SHA256 can additionally pin
an expected 64-hex content digest. All settings are recorded in provenance.
EOF
}

SELF_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "${SELF_DIR}/../../../.." && pwd)
OUT=b0_tracking_benchmark
RECO=""
SIM=""
COUNTERS=""
NEVENTS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --reco) RECO=$2; shift 2 ;;
    --sim) SIM=$2; shift 2 ;;
    --counters) COUNTERS=$2; shift 2 ;;
    --out) OUT=$2; shift 2 ;;
    --nevents) NEVENTS=$2; shift 2 ;;
    --self-test) python3 "${SELF_DIR}/analyze_b0_tracking.py" --self-test; exit $? ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

mkdir -p "${OUT}"
OUT=$(cd "${OUT}" && pwd)

sha_file() {
  local file=$1
  if [[ -n ${file} && -f ${file} ]]; then
    sha256sum "${file}" | awk '{print $1}'
  else
    echo "unresolved"
  fi
}

resolve_detector_asset() {
  local file=$1
  if [[ -z ${file} ]]; then
    echo ""
  elif [[ -f ${file} ]]; then
    realpath "${file}"
  elif [[ -n ${DETECTOR_PATH:-} && -f ${DETECTOR_PATH}/${file} ]]; then
    realpath "${DETECTOR_PATH}/${file}"
  else
    echo ""
  fi
}

DETECTOR_CONFIG_NAME=${DETECTOR_CONFIG:-epic_ip6_extended}
COMPACT_FILE=${DETECTOR_PATH:-}/${DETECTOR_CONFIG_NAME}.xml
MATERIAL_MAP_OVERRIDE=${ACTS_MATERIAL_MAP:-}
EXPECTED_MATERIAL_MAP_SHA256=${ACTS_MATERIAL_MAP_SHA256:-}
STRICT_MATERIAL_MAP=${B0_STRICT_MATERIAL_MAP:-1}
GEOMETRY_DECLARED_MATERIAL_MAP=""

if [[ ${STRICT_MATERIAL_MAP} != 0 && ${STRICT_MATERIAL_MAP} != 1 ]]; then
  echo "B0_STRICT_MATERIAL_MAP must be 0 or 1" >&2
  exit 2
fi
if [[ -n ${EXPECTED_MATERIAL_MAP_SHA256} && ! ${EXPECTED_MATERIAL_MAP_SHA256} =~ ^[0-9A-Fa-f]{64}$ ]]; then
  echo "ACTS_MATERIAL_MAP_SHA256 must be exactly 64 hexadecimal digits" >&2
  exit 2
fi

# The generated ePIC XML contains the material-map constant when the selected
# detector configuration declares a validated map. Read it for provenance only;
# EICrecon itself remains the authority that loads the DD4hep geometry and map.
if [[ -n ${DETECTOR_PATH:-} && -f ${COMPACT_FILE} ]]; then
  GEOMETRY_DECLARED_MATERIAL_MAP=$(python3 - "${COMPACT_FILE}" <<'PY'
import re
import sys
from pathlib import Path

text = Path(sys.argv[1]).read_text(errors="replace")
match = re.search(r'<constant\s+name=["\']material-map["\']\s+value=["\']([^"\']+)["\']', text)
print(match.group(1) if match else "")
PY
)
fi

if [[ -n ${MATERIAL_MAP_OVERRIDE} ]]; then
  MATERIAL_MAP_POLICY="explicit ACTS_MATERIAL_MAP override"
  MATERIAL_MAP_SELECTED=${MATERIAL_MAP_OVERRIDE}
else
  MATERIAL_MAP_POLICY="geometry-declared map"
  MATERIAL_MAP_SELECTED=${GEOMETRY_DECLARED_MATERIAL_MAP:-geometry-declared}
fi
MATERIAL_MAP_SELECTED_RESOLVED=$(resolve_detector_asset "${MATERIAL_MAP_SELECTED}")
GEOMETRY_MATERIAL_MAP_RESOLVED=$(resolve_detector_asset "${GEOMETRY_DECLARED_MATERIAL_MAP}")

{
  echo "{"
  echo "  \"eicrecon_git_sha\": \"$(git -C "${REPO_ROOT}" rev-parse HEAD 2>/dev/null || echo unknown)\","
  echo "  \"eicrecon_git_describe\": \"$(git -C "${REPO_ROOT}" describe --always --dirty 2>/dev/null || echo unknown)\","
  echo "  \"detector_git_sha\": \"$(git -C "${DETECTOR_PATH:-${HOME}/eic/epic}" rev-parse HEAD 2>/dev/null || echo unknown)\","
  echo "  \"DETECTOR_PATH\": \"${DETECTOR_PATH:-}\","
  echo "  \"DETECTOR_CONFIG\": \"${DETECTOR_CONFIG_NAME}\","
  echo "  \"compact_file\": \"${COMPACT_FILE}\","
  echo "  \"beamline_field_file\": \"compact/fields/beamline_18x275.xml (included by generated epic_ip6_extended.xml; B0PF_Bmax matches beamline_5x41.xml)\","
  echo "  \"material_map_policy\": \"${MATERIAL_MAP_POLICY}\","
  echo "  \"strict_material_map_contract\": $([[ ${STRICT_MATERIAL_MAP} == 1 ]] && echo true || echo false),"
  echo "  \"strict_required_detector_constant\": \"B0Tracker_Station_1_ID\","
  echo "  \"geometry_declared_material_map\": \"${GEOMETRY_DECLARED_MATERIAL_MAP:-unknown}\","
  echo "  \"geometry_declared_material_map_resolved\": \"${GEOMETRY_MATERIAL_MAP_RESOLVED}\","
  echo "  \"geometry_declared_material_map_sha256\": \"$(sha_file "${GEOMETRY_MATERIAL_MAP_RESOLVED}")\","
  echo "  \"material_map_override\": \"${MATERIAL_MAP_OVERRIDE}\","
  echo "  \"material_map_selected\": \"${MATERIAL_MAP_SELECTED}\","
  echo "  \"material_map_selected_resolved\": \"${MATERIAL_MAP_SELECTED_RESOLVED}\","
  echo "  \"material_map_selected_sha256\": \"$(sha_file "${MATERIAL_MAP_SELECTED_RESOLVED}")\","
  echo "  \"expected_material_map_sha256\": \"${EXPECTED_MATERIAL_MAP_SHA256}\","
  echo "  \"eicrecon\": \"$(command -v eicrecon 2>/dev/null || echo missing)\","
  echo "  \"eicrecon_version\": \"$(eicrecon --version 2>/dev/null | tr '\n' ' ' || echo unknown)\","
  echo "  \"station_definition\": \"ion-frame z single-linkage clustering, gap 50 mm, min 3 physical stations\","
  echo "  \"matching\": \"MCRecoTrackParticleAssociation weight >= 0.5; seed match via raw-hit associations\","
  echo "  \"reference_surface\": \"origin perigee (0,0,0)\","
  echo "  \"sample\": \"prompt forward protons, 5x41, p=8-41 GeV, ion-frame theta=4-22 mrad\","
  echo "  \"historical_yields\": \"config comments quote 698/765 stub and 723/765 truth on seeded four-station events; not the reconstructible denominator\""
  echo "}"
} > "${OUT}/provenance.json"

if [[ -n ${SIM} ]]; then
  if ! command -v eicrecon >/dev/null 2>&1; then
    echo "eicrecon not in PATH; source install/bin/eicrecon-this.sh inside eic-shell" >&2
    exit 1
  fi
  RECO=${OUT}/b0_benchmark.edm4eic.root
  COUNTERS=${OUT}/counters.json
  export B0_TRACKING_COUNTERS_FILE=${COUNTERS}
  eicrecon_args=(
    -Ppodio:output_file="${RECO}"
    -Ppodio:output_collections=MCParticles,B0TrackerHits,B0TrackerRecHits,B0TrackerRawHitAssociations,B0TrackerSeeds,B0TrackerSeedParameters,B0TrackerCKFTracks,B0TrackerCKFTrackParameters,B0TrackerCKFTrackAssociations,B0TrackerCKFTrajectories,B0TrackerCKFTruthSeededTracks,B0TrackerCKFTruthSeededTrackParameters,B0TrackerCKFTruthSeededTrackAssociations
  )
  if [[ ${STRICT_MATERIAL_MAP} == 1 ]]; then
    eicrecon_args+=(
      -Pacts:RequireGeometryMaterialMap=true
      -Pacts:RequireMaterialCoverage=true
      -Pacts:RequiredMaterialDetectorConstants=B0Tracker_Station_1_ID
    )
  fi
  if [[ -n ${EXPECTED_MATERIAL_MAP_SHA256} ]]; then
    eicrecon_args+=(-Pacts:ExpectedMaterialMapSHA256="${EXPECTED_MATERIAL_MAP_SHA256}")
  fi
  if [[ -n ${MATERIAL_MAP_OVERRIDE} ]]; then
    eicrecon_args+=(-Pacts:MaterialMap="${MATERIAL_MAP_OVERRIDE}")
  fi
  if [[ ${NEVENTS} -gt 0 ]]; then
    eicrecon_args+=(-Pjana:nevents="${NEVENTS}")
  fi
  eicrecon "${eicrecon_args[@]}" "${SIM}" | tee "${OUT}/eicrecon.log"
fi

if [[ -z ${RECO} ]]; then
  echo "need --reco or --sim" >&2
  usage
  exit 2
fi

analyze_args=(--reco "${RECO}" --out "${OUT}" --provenance "${OUT}/provenance.json")
if [[ -n ${COUNTERS} ]]; then
  analyze_args+=(--counters "${COUNTERS}")
elif [[ -f ${OUT}/counters.json ]]; then
  analyze_args+=(--counters "${OUT}/counters.json")
fi

python3 "${SELF_DIR}/analyze_b0_tracking.py" "${analyze_args[@]}"
echo "benchmark report: ${OUT}/report.md"
