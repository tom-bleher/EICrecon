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
  if [[ -f ${file} ]]; then
    sha256sum "${file}" | awk '{print $1}'
  else
    echo "missing"
  fi
}

material_map=${ACTS_MATERIAL_MAP:-${REPO_ROOT}/calibrations/materials-map-ip6-extended-cbc0605e892b.cbor}
if [[ ! -e ${material_map} ]]; then
  material_map=${REPO_ROOT}/calibrations/materials-map.cbor
fi

{
  echo "{"
  echo "  \"eicrecon_git_sha\": \"$(git -C "${REPO_ROOT}" rev-parse HEAD 2>/dev/null || echo unknown)\","
  echo "  \"eicrecon_git_describe\": \"$(git -C "${REPO_ROOT}" describe --always --dirty 2>/dev/null || echo unknown)\","
  echo "  \"detector_git_sha\": \"$(git -C "${DETECTOR_PATH:-${HOME}/eic/epic}" rev-parse HEAD 2>/dev/null || echo unknown)\","
  echo "  \"DETECTOR_PATH\": \"${DETECTOR_PATH:-}\","
  echo "  \"DETECTOR_CONFIG\": \"${DETECTOR_CONFIG:-epic_ip6_extended}\","
  echo "  \"compact_file\": \"${DETECTOR_PATH:-}/${DETECTOR_CONFIG:-epic_ip6_extended}.xml\","
  echo "  \"beamline_field_file\": \"compact/fields/beamline_18x275.xml (included by generated epic_ip6_extended.xml; B0PF_Bmax matches beamline_5x41.xml)\","
  echo "  \"material_map\": \"${material_map}\","
  echo "  \"material_map_sha256\": \"$(sha_file "${material_map}")\","
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
    -Pacts:MaterialMap="${material_map}"
    -Ppodio:output_collections=MCParticles,B0TrackerHits,B0TrackerRecHits,B0TrackerRawHitAssociations,B0TrackerSeeds,B0TrackerSeedParameters,B0TrackerCKFTracks,B0TrackerCKFTrackParameters,B0TrackerCKFTrackAssociations,B0TrackerCKFTrajectories,B0TrackerCKFTruthSeededTracks,B0TrackerCKFTruthSeededTrackParameters,B0TrackerCKFTruthSeededTrackAssociations
  )
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
