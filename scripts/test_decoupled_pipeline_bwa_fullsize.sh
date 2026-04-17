#!/usr/bin/env bash
# =============================================================================
# test_decoupled_pipeline_bwa_fullsize.sh
#
# Runs the original joint (-W 1) BWA-only CIRI3 pipeline AND the four-stage
# decoupled pipeline (SCAN1 → BUILD_UNIVERSE → SCAN2 → FINALIZE) on full-size
# BWA data, then verifies that both pipelines produce identical BSJ and FSJ
# matrices.
#
# Uses the same bwa.sam files that were produced during the STAR re-alignment
# step, but runs only the BWA path (no -Ma 1).
#
# Each stage is skipped if its outputs already exist, so the script is safe
# to re-run after a partial failure.
#
# Usage:
#   bash scripts/test_decoupled_pipeline_bwa_fullsize.sh [options]
#
#   --threads N      threads for each per-sample stage (default: 8)
#   --output-dir D   output directory (default: DATA_DIR/decoupled_bwa_comparison)
#   --keep           keep output directory after the run
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
THREADS=8
KEEP=0
OUT_ROOT=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)     THREADS="$2";  shift 2 ;;
        --output-dir)  OUT_ROOT="$2"; shift 2 ;;
        --keep)        KEEP=1;        shift   ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Paths — edit to match your environment
# ---------------------------------------------------------------------------
DATA_DIR=/pastel/tools/circRNA_tools/test_data
REF_FA=${DATA_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa
GTF_FILE=${DATA_DIR}/gencode.v32.primary_assembly.annotation.gtf

SAMPLES=(
    "Div_100_S91"
    "Div_101_S92"
    "PARDOS_1_S1"
    "PARDOS_2_S2"
)

# ---------------------------------------------------------------------------
# Derived paths
# ---------------------------------------------------------------------------
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN_DIR="${REPO_ROOT}/bin"
LIB_DIR="${REPO_ROOT}/lib"
CLASSPATH="${BIN_DIR}:${LIB_DIR}/htsjdk-3.0.4.jar"

[[ -z "$OUT_ROOT" ]] && OUT_ROOT="${DATA_DIR}/decoupled_bwa_comparison"

ORIG_DIR="${OUT_ROOT}/original"
DECOUPLED_DIR="${OUT_ROOT}/decoupled"
SCAN1_DIR="${DECOUPLED_DIR}/scan1"
UNIVERSE_DIR="${DECOUPLED_DIR}/universe"
SCAN2_DIR="${DECOUPLED_DIR}/scan2"
FINALIZE_DIR="${DECOUPLED_DIR}/finalize"

SCAN1_META_TSV="${UNIVERSE_DIR}/samples_scan1.tsv"
FINALIZE_TSV="${FINALIZE_DIR}/finalize_samples.tsv"

PASS=0
FAIL=0

# ---------------------------------------------------------------------------
# Java commands
# ---------------------------------------------------------------------------
if [[ -n "${CONDA_PREFIX:-}" ]] && [[ -x "${CONDA_PREFIX}/bin/java" ]]; then
    JAVA_BIN="${CONDA_PREFIX}/bin/java"
else
    JAVA_BIN="$(command -v java)"
fi
JAVAC_FLAGS="-source 8 -target 8"

JAVA_SRC="${JAVA_BIN} -cp ${CLASSPATH} com.zx.test.TestParameters"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
info()  { echo "[INFO]  $*"; }
ok()    { echo "[PASS]  $*"; PASS=$((PASS+1)); }
fail()  { echo "[FAIL]  $*"; FAIL=$((FAIL+1)); }
die()   { echo "[ERROR] $*" >&2; exit 1; }

check_exists() {
    local label="$1" f="$2"
    if [[ -s "$f" ]]; then ok "$label exists and is non-empty"
    else fail "$label missing or empty: $f"; fi
}

cleanup() {
    if [[ $KEEP -eq 0 ]]; then
        rm -rf "$OUT_ROOT"
        info "Output directory removed."
    else
        info "Outputs kept at: $OUT_ROOT"
    fi
}
trap cleanup EXIT

# ---------------------------------------------------------------------------
# 0. Pre-flight checks
# ---------------------------------------------------------------------------
info "=== Pre-flight checks ==="
[[ -f "$REF_FA"   ]] || die "Reference FASTA not found: $REF_FA"
[[ -f "$GTF_FILE" ]] || die "GTF not found: $GTF_FILE"
[[ -x "$JAVA_BIN" ]] || die "Java not found at: $JAVA_BIN (is the CIRI3 conda env active?)"

for S in "${SAMPLES[@]}"; do
    # STAR_DIR="${DATA_DIR}/STAR_output_${S}"
    [[ -f "${DATA_DIR}/${S}.sam" ]] || die "Missing: ${S}.sam"
done
info "All input files found for ${#SAMPLES[@]} samples."
info "Output directory: ${OUT_ROOT}"

mkdir -p "$ORIG_DIR" "$SCAN1_DIR" "$UNIVERSE_DIR" "$SCAN2_DIR" "$FINALIZE_DIR"

# ---------------------------------------------------------------------------
# 1. Compile from source
# ---------------------------------------------------------------------------
info "=== Compiling from source (target: Java 8) ==="
find "${REPO_ROOT}/src" -name "*.java" > /tmp/ciri3_sources_$$.txt
# shellcheck disable=SC2086
javac ${JAVAC_FLAGS} -cp "${CLASSPATH}" -d "${BIN_DIR}" \
    @/tmp/ciri3_sources_$$.txt 2>&1 \
    | grep -v "^\(warning\|Note\)" || true
rm -f /tmp/ciri3_sources_$$.txt
info "Compilation complete."

# ---------------------------------------------------------------------------
# 2. ORIGINAL pipeline (-W 1, BWA-only, compiled source)
# ---------------------------------------------------------------------------
info "=== Stage 0: ORIGINAL pipeline (-W 1, BWA-only, source) ==="

ORIG_BSJ="${ORIG_DIR}/result.BSJ_Matrix"
if [[ -s "$ORIG_BSJ" ]]; then
    info "  [SKIP] Original results exist at ${ORIG_DIR}/result"
else
    ORIG_TSV="${ORIG_DIR}/samples.tsv"
    > "$ORIG_TSV"
    for S in "${SAMPLES[@]}"; do
        # STAR_DIR="${DATA_DIR}/STAR_output_${S}"
        echo "${DATA_DIR}/${S}.sam" >> "$ORIG_TSV"
    done

    ${JAVA_SRC} \
        -I  "${ORIG_TSV}" \
        -O  "${ORIG_DIR}/result" \
        -F  "${REF_FA}" \
        -A  "${GTF_FILE}" \
        -W  1 -T "${THREADS}" -S 0 \
        2>&1 | tee "${ORIG_DIR}/run.log" \
        | grep -E "CIRI3|scan|completed|circRNA|Mapped|time" || true
fi

check_exists "Original BSJ_Matrix" "${ORIG_DIR}/result.BSJ_Matrix"
check_exists "Original FSJ_Matrix" "${ORIG_DIR}/result.FSJ_Matrix"
ORIG_CIRCS=$(tail -n +2 "${ORIG_DIR}/result.BSJ_Matrix" | wc -l)
info "Original pipeline: ${ORIG_CIRCS} circRNAs."

# ---------------------------------------------------------------------------
# 3. DECOUPLED pipeline (compiled source, BWA-only)
# ---------------------------------------------------------------------------
info "=== Decoupled pipeline (compiled source, BWA-only) ==="

# Reset TSV files each run
> "$SCAN1_META_TSV"
> "$FINALIZE_TSV"

# --- Stage 1: SCAN1 (per sample) ---
info "--- Stage 1: SCAN1 ---"
for S in "${SAMPLES[@]}"; do
    # STAR_DIR="${DATA_DIR}/STAR_output_${S}"
    BWA_SAM="${DATA_DIR}/${S}.sam"
    OUT_PREFIX="${SCAN1_DIR}/${S}"
    META="${OUT_PREFIX}.scan1_meta"

    if [[ -s "$META" ]]; then
        info "  [SKIP] SCAN1 already done for $S"
    else
        info "  SCAN1: $S"
        ${JAVA_SRC} SCAN1 \
            -I  "${BWA_SAM}" \
            -O  "${OUT_PREFIX}" \
            -F  "${REF_FA}" \
            -A  "${GTF_FILE}" \
            -T  "${THREADS}" -S 0 \
            2>&1 | grep -E "scan|meta|time|Mapped" || true
    fi

    check_exists "SCAN1 meta ($S)" "${META}"

    SPLIT_NUM=$(grep "^fileSplitNum=" "${META}" | cut -d= -f2)
    local_fail=0
    for i in $(seq 1 "${SPLIT_NUM}"); do
        [[ -f "${BWA_SAM}BSJ${i}" ]] || { fail "Missing BSJ file: ${BWA_SAM}BSJ${i}"; local_fail=1; }
    done
    [[ $local_fail -eq 0 ]] && ok "SCAN1 BSJ files present for $S (${SPLIT_NUM} splits)"

    echo -e "${BWA_SAM}\t${META}" >> "$SCAN1_META_TSV"
done

# --- Stage 2a: BUILD_UNIVERSE ---
info "--- Stage 2a: BUILD_UNIVERSE ---"
UNIVERSE_FILE="${UNIVERSE_DIR}/cohort.universe"
if [[ -s "$UNIVERSE_FILE" ]]; then
    info "  [SKIP] Universe already exists"
else
    ${JAVA_SRC} BUILD_UNIVERSE \
        -I  "${SCAN1_META_TSV}" \
        -F  "${REF_FA}" \
        -O  "${UNIVERSE_DIR}/cohort" \
        2>&1 | grep -E "Universe|circRNA|time" || true
fi
check_exists "Universe file" "${UNIVERSE_FILE}"
UNIVERSE_CIRCS=$(grep -c "^chr" "${UNIVERSE_FILE}" || true)
info "Universe: ${UNIVERSE_CIRCS} circRNA candidates."

# --- Stage 3: SCAN2 (per sample) ---
info "--- Stage 3: SCAN2 ---"
for S in "${SAMPLES[@]}"; do
    # STAR_DIR="${DATA_DIR}/STAR_output_${S}"
    BWA_SAM="${DATA_DIR}/${S}.sam"
    META="${SCAN1_DIR}/${S}.scan1_meta"
    SPLIT_NUM=$(grep "^fileSplitNum=" "${META}" | cut -d= -f2)
    OUT_PREFIX="${SCAN2_DIR}/${S}"
    FSJ_COUNTS="${OUT_PREFIX}.fsj_counts"

    if [[ -s "$FSJ_COUNTS" ]]; then
        info "  [SKIP] SCAN2 already done for $S"
    else
        info "  SCAN2: $S"
        ${JAVA_SRC} SCAN2 \
            -I  "${BWA_SAM}" \
            -CU "${UNIVERSE_FILE}" \
            -O  "${OUT_PREFIX}" \
            -F  "${REF_FA}" \
            -T  "${THREADS}" \
            2>&1 | grep -E "scan|FSJ|BSJ|time" || true
    fi
    check_exists "SCAN2 FSJ counts ($S)" "${FSJ_COUNTS}"

    echo -e "${BWA_SAM}\t${FSJ_COUNTS}\t${SPLIT_NUM}\t${S}" >> "$FINALIZE_TSV"
done

# --- Stage 4: FINALIZE ---
info "--- Stage 4: FINALIZE ---"
FINAL_BSJ="${FINALIZE_DIR}/result.BSJ_Matrix"
if [[ -s "$FINAL_BSJ" ]]; then
    info "  [SKIP] FINALIZE already done"
else
    ${JAVA_SRC} FINALIZE \
        -I  "${FINALIZE_TSV}" \
        -F  "${REF_FA}" \
        -O  "${FINALIZE_DIR}/result" \
        -A  "${GTF_FILE}" \
        -S  0 \
        2>&1 | grep -E "FINALIZE|Summary|Matrix|circRNA|time" || true
fi
check_exists "Decoupled BSJ_Matrix" "${FINAL_BSJ}"
check_exists "Decoupled FSJ_Matrix" "${FINALIZE_DIR}/result.FSJ_Matrix"
DECOUPLED_CIRCS=$(tail -n +2 "${FINAL_BSJ}" | wc -l)
info "Decoupled pipeline: ${DECOUPLED_CIRCS} circRNAs."

# ---------------------------------------------------------------------------
# 4. Compare outputs
# ---------------------------------------------------------------------------
info "=== Comparing outputs ==="

normalise_matrix() {
    local f="$1"
    awk 'NR==1{
        printf "circRNA_ID"
        for(i=2;i<=NF;i++) printf "\ts%d", i-1
        printf "\n"; next
    }
    { print | "sort" }' "$f"
}

A_BSJ=$(normalise_matrix "${ORIG_DIR}/result.BSJ_Matrix")
B_BSJ=$(normalise_matrix "${FINAL_BSJ}")
if [[ "$A_BSJ" == "$B_BSJ" ]]; then
    ok "BSJ_Matrix: original and decoupled are IDENTICAL (${ORIG_CIRCS} circRNAs)"
else
    fail "BSJ_Matrix: original and decoupled DIFFER"
    diff <(echo "$A_BSJ") <(echo "$B_BSJ") | head -40
fi

A_FSJ=$(normalise_matrix "${ORIG_DIR}/result.FSJ_Matrix")
B_FSJ=$(normalise_matrix "${FINALIZE_DIR}/result.FSJ_Matrix")
if [[ "$A_FSJ" == "$B_FSJ" ]]; then
    ok "FSJ_Matrix: original and decoupled are IDENTICAL (${ORIG_CIRCS} circRNAs)"
else
    fail "FSJ_Matrix: original and decoupled DIFFER"
    diff <(echo "$A_FSJ") <(echo "$B_FSJ") | head -40
fi

info "Checking universe coverage..."
ORIG_IDS=$(tail -n +2 "${ORIG_DIR}/result.BSJ_Matrix" | awk '{print $1}' | sort)
UNIV_IDS=$(grep "^chr" "${UNIVERSE_FILE}" | awk '{printf "%s:%s|%s\n", $1, $2, $3}' | sort)
MISSING=$(comm -23 <(echo "$ORIG_IDS") <(echo "$UNIV_IDS") | wc -l)
if [[ "$MISSING" -eq 0 ]]; then
    ok "Universe coverage: all original circRNAs present in universe"
else
    fail "Universe coverage: $MISSING original circRNAs missing from universe"
    comm -23 <(echo "$ORIG_IDS") <(echo "$UNIV_IDS") | head -10
fi

# ---------------------------------------------------------------------------
# 5. Summary
# ---------------------------------------------------------------------------
echo ""
echo "========================================"
echo "  TEST SUMMARY"
echo "========================================"
printf "  Original pipeline : %d circRNAs\n" "${ORIG_CIRCS}"
printf "  Decoupled pipeline: %d circRNAs\n" "${DECOUPLED_CIRCS}"
echo "  PASSED : $PASS"
echo "  FAILED : $FAIL"
echo "========================================"
if [[ $FAIL -eq 0 ]]; then
    echo "  ALL TESTS PASSED"
    exit 0
else
    echo "  SOME TESTS FAILED — see output above"
    exit 1
fi
