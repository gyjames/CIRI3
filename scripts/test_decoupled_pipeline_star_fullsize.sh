#!/usr/bin/env bash
# =============================================================================
# test_decoupled_pipeline_star_fullsize.sh
#
# Runs the original joint (-W 1) CIRI3 pipeline AND the four-stage decoupled
# pipeline (SCAN1 → BUILD_UNIVERSE → SCAN2 → FINALIZE) on full-size STAR data,
# then verifies that both pipelines produce identical BSJ and FSJ matrices.
#
# Assumes alignment is already complete (STAR + BWA-MEM re-alignment of
# unmapped reads) and STAR output directories exist under DATA_DIR.
#
# Usage:
#   bash scripts/test_decoupled_pipeline_star_fullsize.sh [--threads N] [--keep]
#
#   --threads N   threads for each per-sample stage (default: 8)
#   --keep        keep output directories after the run
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
THREADS=8
KEEP=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --keep)    KEEP=1; shift ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Paths — edit to match your environment
# ---------------------------------------------------------------------------
DATA_DIR=/pastel/tools/circRNA_tools/test_data
REF_FA=${DATA_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa
GTF_FILE=${DATA_DIR}/gencode.v32.primary_assembly.annotation.gtf

# Original pipeline: use the pre-built JAR
CIRI3_JAR=/pastel/tools/circRNA_tools/CIRI3/CIRI3_Java_18.0.1.jar
JAVA_JAR="$CONDA_PREFIX/bin/java -jar ${CIRI3_JAR}"

# Decoupled pipeline: use this repository's compiled source
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN_DIR="${REPO_ROOT}/bin"
LIB_DIR="${REPO_ROOT}/lib"
CLASSPATH="${BIN_DIR}:${LIB_DIR}/htsjdk-3.0.4.jar"
JAVA_SRC="java -cp ${CLASSPATH} com.zx.test.TestParameters"

SAMPLES=(
    "Div_100_S91"
    "Div_101_S92"
)

OUT_ROOT="${DATA_DIR}/decoupled_comparison_$$"
ORIG_DIR="${OUT_ROOT}/original"
DECOUPLED_DIR="${OUT_ROOT}/decoupled"
SCAN1_DIR="${DECOUPLED_DIR}/scan1"
UNIVERSE_DIR="${DECOUPLED_DIR}/universe"
SCAN2_DIR="${DECOUPLED_DIR}/scan2"
FINALIZE_DIR="${DECOUPLED_DIR}/finalize"

PASS=0
FAIL=0

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
[[ -f "$CIRI3_JAR" ]] || die "CIRI3 JAR not found: $CIRI3_JAR"
[[ -d "$BIN_DIR"  ]] || die "Compiled bin/ not found at $BIN_DIR — run: javac -source 11 -target 11 -cp \"bin:lib/htsjdk-3.0.4.jar\" -d bin @sources.txt"

for S in "${SAMPLES[@]}"; do
    STAR_DIR="${DATA_DIR}/STAR_output_${S}"
    [[ -f "${STAR_DIR}/Chimeric.out.junction" ]] || die "Missing: ${STAR_DIR}/Chimeric.out.junction"
    [[ -f "${STAR_DIR}/Aligned.out.sam"       ]] || die "Missing: ${STAR_DIR}/Aligned.out.sam"
    [[ -f "${STAR_DIR}/bwa.sam"               ]] || die "Missing: ${STAR_DIR}/bwa.sam"
done
info "All input files found for ${#SAMPLES[@]} samples."

mkdir -p "$ORIG_DIR" "$SCAN1_DIR" "$UNIVERSE_DIR" "$SCAN2_DIR" "$FINALIZE_DIR"

# ---------------------------------------------------------------------------
# 1. Recompile from source (ensures fixes are in bin/)
# ---------------------------------------------------------------------------
info "=== Compiling from source ==="
find "${REPO_ROOT}/src" -name "*.java" > /tmp/ciri3_sources_$$.txt
javac -source 11 -target 11 -cp "${CLASSPATH}" -d "${BIN_DIR}" \
    @/tmp/ciri3_sources_$$.txt 2>&1 \
    | grep -v "^\(warning\|Note\)" || true
rm -f /tmp/ciri3_sources_$$.txt
info "Compilation complete."

# ---------------------------------------------------------------------------
# 2. Run ORIGINAL pipeline (-W 1, STAR mode) using the JAR
# ---------------------------------------------------------------------------
info "=== Running ORIGINAL pipeline (-W 1, STAR, JAR) ==="

ORIG_TSV="${ORIG_DIR}/samples.tsv"
> "$ORIG_TSV"
for S in "${SAMPLES[@]}"; do
    STAR_DIR="${DATA_DIR}/STAR_output_${S}"
    echo "${STAR_DIR}/Chimeric.out.junction,${STAR_DIR}/Aligned.out.sam,${STAR_DIR}/bwa.sam" \
        >> "$ORIG_TSV"
done

${JAVA_JAR} \
    -I  "${ORIG_TSV}" \
    -O  "${ORIG_DIR}/result" \
    -F  "${REF_FA}" \
    -A  "${GTF_FILE}" \
    -W  1 -Ma 1 -T "${THREADS}" -S 0 \
    2>&1 | tee "${ORIG_DIR}/run.log" \
    | grep -E "CIRI3|scan|completed|circRNA|Matrix|Mapped|time" || true

check_exists "Original BSJ_Matrix" "${ORIG_DIR}/result.BSJ_Matrix"
check_exists "Original FSJ_Matrix" "${ORIG_DIR}/result.FSJ_Matrix"
ORIG_CIRCS=$(tail -n +2 "${ORIG_DIR}/result.BSJ_Matrix" | wc -l)
info "Original pipeline detected ${ORIG_CIRCS} circRNAs."

# ---------------------------------------------------------------------------
# 3. Run DECOUPLED pipeline (compiled source, STAR mode)
# ---------------------------------------------------------------------------
info "=== Running DECOUPLED pipeline (compiled source, STAR mode) ==="

SCAN1_META_TSV="${UNIVERSE_DIR}/samples_scan1.tsv"
FINALIZE_TSV="${FINALIZE_DIR}/finalize_samples.tsv"
> "$SCAN1_META_TSV"
> "$FINALIZE_TSV"

# --- Stage 1: SCAN1 (per sample) ---
info "--- Stage 1: SCAN1 ---"
for S in "${SAMPLES[@]}"; do
    STAR_DIR="${DATA_DIR}/STAR_output_${S}"
    TRIPLE="${STAR_DIR}/Chimeric.out.junction,${STAR_DIR}/Aligned.out.sam,${STAR_DIR}/bwa.sam"
    BWA_SAM="${STAR_DIR}/bwa.sam"
    OUT_PREFIX="${SCAN1_DIR}/${S}"

    info "  SCAN1: $S"
    ${JAVA_SRC} SCAN1 \
        -I  "${TRIPLE}" \
        -O  "${OUT_PREFIX}" \
        -F  "${REF_FA}" \
        -A  "${GTF_FILE}" \
        -T  "${THREADS}" -Ma 1 -S 0 \
        2>&1 | grep -E "scan|meta|time|Mapped" || true

    META="${OUT_PREFIX}.scan1_meta"
    check_exists "SCAN1 meta ($S)" "${META}"
    SPLIT_NUM=$(grep "^fileSplitNum=" "${META}" | cut -d= -f2)
    for i in $(seq 1 "${SPLIT_NUM}"); do
        [[ -f "${BWA_SAM}BSJ${i}" ]] || fail "Missing BSJ file: ${BWA_SAM}BSJ${i}"
    done
    ok "SCAN1 BSJ files present for $S (${SPLIT_NUM} splits)"

    # NOTE: for STAR mode the BSJ files live next to bwa.sam, so that path
    # is used as the samFile key in BUILD_UNIVERSE and FINALIZE
    echo -e "${BWA_SAM}\t${META}" >> "$SCAN1_META_TSV"
done

# --- Stage 2a: BUILD_UNIVERSE ---
info "--- Stage 2a: BUILD_UNIVERSE ---"
${JAVA_SRC} BUILD_UNIVERSE \
    -I  "${SCAN1_META_TSV}" \
    -F  "${REF_FA}" \
    -O  "${UNIVERSE_DIR}/cohort" \
    2>&1 | grep -E "Universe|circRNA|time" || true
check_exists "Universe file" "${UNIVERSE_DIR}/cohort.universe"
UNIVERSE_CIRCS=$(grep -c "^chr" "${UNIVERSE_DIR}/cohort.universe" || true)
info "Universe contains ${UNIVERSE_CIRCS} circRNA candidates."

# --- Stage 3: SCAN2 (per sample) ---
info "--- Stage 3: SCAN2 ---"
for S in "${SAMPLES[@]}"; do
    STAR_DIR="${DATA_DIR}/STAR_output_${S}"
    TRIPLE="${STAR_DIR}/Chimeric.out.junction,${STAR_DIR}/Aligned.out.sam,${STAR_DIR}/bwa.sam"
    BWA_SAM="${STAR_DIR}/bwa.sam"
    META="${SCAN1_DIR}/${S}.scan1_meta"
    SPLIT_NUM=$(grep "^fileSplitNum=" "${META}" | cut -d= -f2)
    OUT_PREFIX="${SCAN2_DIR}/${S}"

    info "  SCAN2: $S"
    ${JAVA_SRC} SCAN2 \
        -I  "${TRIPLE}" \
        -CU "${UNIVERSE_DIR}/cohort.universe" \
        -O  "${OUT_PREFIX}" \
        -F  "${REF_FA}" \
        -T  "${THREADS}" -Ma 1 \
        2>&1 | grep -E "scan|FSJ|time" || true
    check_exists "SCAN2 FSJ counts ($S)" "${OUT_PREFIX}.fsj_counts"

    echo -e "${BWA_SAM}\t${OUT_PREFIX}.fsj_counts\t${SPLIT_NUM}\t${S}" >> "$FINALIZE_TSV"
done

# --- Stage 4: FINALIZE ---
info "--- Stage 4: FINALIZE ---"
${JAVA_SRC} FINALIZE \
    -I  "${FINALIZE_TSV}" \
    -F  "${REF_FA}" \
    -O  "${FINALIZE_DIR}/result" \
    -A  "${GTF_FILE}" \
    -S  0 \
    2>&1 | grep -E "FINALIZE|Summary|Matrix|circRNA|time" || true
check_exists "Decoupled BSJ_Matrix" "${FINALIZE_DIR}/result.BSJ_Matrix"
check_exists "Decoupled FSJ_Matrix" "${FINALIZE_DIR}/result.FSJ_Matrix"
DECOUPLED_CIRCS=$(tail -n +2 "${FINALIZE_DIR}/result.BSJ_Matrix" | wc -l)
info "Decoupled pipeline detected ${DECOUPLED_CIRCS} circRNAs."

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
B_BSJ=$(normalise_matrix "${FINALIZE_DIR}/result.BSJ_Matrix")
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

# Universe coverage: every circRNA in original output should be in the universe
info "Checking universe coverage..."
ORIG_CIRCS_IDS=$(tail -n +2 "${ORIG_DIR}/result.BSJ_Matrix" | awk '{print $1}' | sort)
UNIVERSE_IDS=$(grep "^chr" "${UNIVERSE_DIR}/cohort.universe" | \
    awk '{printf "%s:%s|%s\n", $1, $2, $3}' | sort)
MISSING=$(comm -23 \
    <(echo "$ORIG_CIRCS_IDS") \
    <(echo "$UNIVERSE_IDS") | wc -l)
if [[ "$MISSING" -eq 0 ]]; then
    ok "Universe coverage: all original circRNAs present in universe"
else
    fail "Universe coverage: $MISSING circRNAs from original missing from universe"
    comm -23 <(echo "$ORIG_CIRCS_IDS") <(echo "$UNIVERSE_IDS") | head -10
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
