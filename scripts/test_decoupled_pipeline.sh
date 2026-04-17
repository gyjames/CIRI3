#!/usr/bin/env bash
# =============================================================================
# test_decoupled_pipeline.sh
#
# Verifies that the four-stage decoupled pipeline (SCAN1 → BUILD_UNIVERSE →
# SCAN2 → FINALIZE) produces identical BSJ and FSJ matrices compared to the
# original joint pipeline (-W 1).
#
# Usage:
#   bash scripts/test_decoupled_pipeline.sh [--threads N] [--keep]
#
#   --threads N   number of threads for each per-sample stage (default: 1)
#   --keep        do not delete the working directories after the test
#
# Requirements:
#   - JDK 11+ (javac / java on PATH)  — any version 11 or newer is supported
#   - The four sample SAM files in data/circRNA/Mutiple/
#   - Reference FASTA at data/circRNA/ref.fa
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
THREADS=1
KEEP=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --keep)    KEEP=1; shift ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Paths — all relative to the repository root
# ---------------------------------------------------------------------------
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

SRC_DIR="$REPO_ROOT/src"
BIN_DIR="$REPO_ROOT/bin"
LIB_DIR="$REPO_ROOT/lib"
DATA_DIR="$REPO_ROOT/data/circRNA/Mutiple"
REF_FA="$REPO_ROOT/data/circRNA/ref.fa"
CLASSPATH="$BIN_DIR:$LIB_DIR/htsjdk-3.0.4.jar"
MAIN_CLASS="com.zx.test.TestParameters"

WORK_DIR="$REPO_ROOT/test_decoupled_$$"
ORIG_DIR="$WORK_DIR/original"
DECOUPLED_DIR="$WORK_DIR/decoupled"
SCAN1_DIR="$DECOUPLED_DIR/scan1"
UNIVERSE_DIR="$DECOUPLED_DIR/universe"
SCAN2_DIR="$DECOUPLED_DIR/scan2"
FINALIZE_DIR="$DECOUPLED_DIR/finalize"

PASS=0
FAIL=0

# ---------------------------------------------------------------------------
# Java commands — use the same JVM for both original and decoupled pipelines
# so there is no version mismatch between the JAR and the compiled source.
# ---------------------------------------------------------------------------
if [[ -n "${CONDA_PREFIX:-}" ]] && [[ -x "${CONDA_PREFIX}/bin/java" ]]; then
    JAVA_BIN="${CONDA_PREFIX}/bin/java"
else
    JAVA_BIN="$(command -v java)"
fi
JAVAC_FLAGS="-source 8 -target 8"   # compile to class-file version 52 (Java 8)

JAVA_SRC="${JAVA_BIN} -cp ${CLASSPATH} com.zx.test.TestParameters"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
info()  { echo "[INFO]  $*"; }
ok()    { echo "[PASS]  $*"; PASS=$((PASS+1)); }
fail()  { echo "[FAIL]  $*"; FAIL=$((FAIL+1)); }
die()   { echo "[ERROR] $*" >&2; exit 1; }

# Compare two matrix files after sorting rows (skip header)
compare_matrices() {
    local label="$1" file_a="$2" file_b="$3"
    if [[ ! -f "$file_a" ]]; then fail "$label: $file_a missing"; return; fi
    if [[ ! -f "$file_b" ]]; then fail "$label: $file_b missing"; return; fi

    # Compare headers (column names may differ in order for sample columns,
    # so we compare sorted column sets)
    local head_a head_b
    head_a=$(head -1 "$file_a" | tr '\t' '\n' | sort | tr '\n' '\t')
    head_b=$(head -1 "$file_b" | tr '\t' '\n' | sort | tr '\n' '\t')
    if [[ "$head_a" != "$head_b" ]]; then
        fail "$label: column headers differ"
        echo "  A: $(head -1 "$file_a")"
        echo "  B: $(head -1 "$file_b")"
        return
    fi

    # Sort data rows by circRNA_ID (col 1) and compare
    local sorted_a sorted_b
    sorted_a=$(tail -n +2 "$file_a" | sort -k1,1)
    sorted_b=$(tail -n +2 "$file_b" | sort -k1,1)
    if [[ "$sorted_a" == "$sorted_b" ]]; then
        local nlines
        nlines=$(echo "$sorted_a" | wc -l)
        ok "$label: $nlines circRNAs match"
    else
        fail "$label: content differs"
        diff <(echo "$sorted_a") <(echo "$sorted_b") | head -20
    fi
}

# Check a file is non-empty
check_exists() {
    local label="$1" f="$2"
    if [[ -s "$f" ]]; then ok "$label exists and is non-empty"
    else fail "$label missing or empty: $f"; fi
}

cleanup() {
    if [[ $KEEP -eq 0 ]]; then
        rm -rf "$WORK_DIR"
        info "Working directory removed."
    else
        info "Working directory kept at: $WORK_DIR"
    fi
}
trap cleanup EXIT

# ---------------------------------------------------------------------------
# 0. Pre-flight checks
# ---------------------------------------------------------------------------
info "=== Pre-flight checks ==="
[[ -d "$DATA_DIR" ]] || die "Test data not found at $DATA_DIR"
[[ -f "$REF_FA"  ]] || die "Reference FASTA not found at $REF_FA"
[[ -x "$JAVA_BIN"  ]] || die "Java not found at: $JAVA_BIN (is the CIRI3 conda env active?)"

SAMPLES=()
for s in 1 2 3 4; do
    SAM="$DATA_DIR/sample${s}.sam"
    [[ -f "$SAM" ]] || die "Missing test SAM: $SAM"
    SAMPLES+=("$SAM")
done
info "Found ${#SAMPLES[@]} sample SAM files."

# ---------------------------------------------------------------------------
# 1. Compile from source
#    Compile to Java 8 class-file format (-source 8 -target 8) so the classes
#    run on the Java 8 JVM bundled with the CIRI3 conda environment.
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
# 2. Build sample list TSVs
# ---------------------------------------------------------------------------
mkdir -p "$ORIG_DIR" "$SCAN1_DIR" "$UNIVERSE_DIR" "$SCAN2_DIR" "$FINALIZE_DIR"

# Original pipeline: one path per line
ORIG_SAMPLES_TSV="$ORIG_DIR/samples.tsv"
> "$ORIG_SAMPLES_TSV"
for SAM in "${SAMPLES[@]}"; do
    echo "$SAM" >> "$ORIG_SAMPLES_TSV"
done

# BUILD_UNIVERSE input: samFile<TAB>metaFile
SCAN1_META_TSV="$UNIVERSE_DIR/samples_scan1.tsv"
> "$SCAN1_META_TSV"

# FINALIZE input: samFile<TAB>fsjCountsFile<TAB>splitNum<TAB>sampleName
FINALIZE_TSV="$FINALIZE_DIR/finalize_samples.tsv"
> "$FINALIZE_TSV"

# ---------------------------------------------------------------------------
# 3. Run ORIGINAL pipeline (-W 1)
# ---------------------------------------------------------------------------
info "=== Running ORIGINAL pipeline (-W 1) ==="
${JAVA_BIN} -cp "${CLASSPATH}" "${MAIN_CLASS}" \
    -I "$ORIG_SAMPLES_TSV" \
    -O "$ORIG_DIR/result" \
    -F "$REF_FA" \
    -W 1 \
    -T "$THREADS" \
    -S 0 \
    2>&1 | tee "$ORIG_DIR/run.log" | grep -E "CIRI3|scan|completed|circRNA|Matrix|time" || true

check_exists "Original BSJ_Matrix" "$ORIG_DIR/result.BSJ_Matrix"
check_exists "Original FSJ_Matrix" "$ORIG_DIR/result.FSJ_Matrix"

# Re-index the BSJ files that were deleted by the original pipeline
# (the original pipeline deletes BSJ files; we use the matrix output as ground truth)
ORIG_BSJ_LINES=$(tail -n +2 "$ORIG_DIR/result.BSJ_Matrix" | wc -l)
info "Original pipeline detected $ORIG_BSJ_LINES circRNAs."

# ---------------------------------------------------------------------------
# 4. Run DECOUPLED pipeline
# ---------------------------------------------------------------------------
info "=== Running DECOUPLED pipeline ==="

# --- Stage 1: SCAN1 (per sample) ---
info "--- Stage 1: SCAN1 ---"
for SAM in "${SAMPLES[@]}"; do
    SAMPLE_NAME=$(basename "$SAM" .sam)
    OUT_PREFIX="$SCAN1_DIR/${SAMPLE_NAME}"
    info "  SCAN1: $SAMPLE_NAME"
    ${JAVA_BIN} -cp "${CLASSPATH}" "${MAIN_CLASS}" SCAN1 \
        -I "$SAM" \
        -O "$OUT_PREFIX" \
        -F "$REF_FA" \
        -T "$THREADS" \
        -S 0 \
        2>&1 | grep -E "scan|completed|meta|time" || true
    check_exists "SCAN1 meta ($SAMPLE_NAME)" "${OUT_PREFIX}.scan1_meta"

    META="${OUT_PREFIX}.scan1_meta"
    SPLIT_NUM=$(grep "^fileSplitNum=" "$META" | cut -d= -f2)
    # Confirm BSJ files were written
    for i in $(seq 1 "$SPLIT_NUM"); do
        [[ -f "${SAM}BSJ${i}" ]] || fail "Missing BSJ file: ${SAM}BSJ${i}"
    done
    ok "SCAN1 BSJ files present for $SAMPLE_NAME"

    echo -e "$SAM\t$META" >> "$SCAN1_META_TSV"
done

# --- Stage 2a: BUILD_UNIVERSE ---
info "--- Stage 2a: BUILD_UNIVERSE ---"
${JAVA_BIN} -cp "${CLASSPATH}" "${MAIN_CLASS}" BUILD_UNIVERSE \
    -I "$SCAN1_META_TSV" \
    -F "$REF_FA" \
    -O "$UNIVERSE_DIR/cohort" \
    2>&1 | grep -E "Universe|circRNA|time" || true
check_exists "Universe file" "$UNIVERSE_DIR/cohort.universe"
UNIVERSE_LINES=$(grep -c "^chr" "$UNIVERSE_DIR/cohort.universe" || true)
info "Universe contains $UNIVERSE_LINES circRNA candidates."

# --- Stage 3: SCAN2 (per sample) ---
info "--- Stage 3: SCAN2 ---"
for SAM in "${SAMPLES[@]}"; do
    SAMPLE_NAME=$(basename "$SAM" .sam)
    SCAN1_META="$SCAN1_DIR/${SAMPLE_NAME}.scan1_meta"
    SPLIT_NUM=$(grep "^fileSplitNum=" "$SCAN1_META" | cut -d= -f2)
    OUT_PREFIX="$SCAN2_DIR/${SAMPLE_NAME}"
    info "  SCAN2: $SAMPLE_NAME"
    ${JAVA_BIN} -cp "${CLASSPATH}" "${MAIN_CLASS}" SCAN2 \
        -I "$SAM" \
        -CU "$UNIVERSE_DIR/cohort.universe" \
        -O "$OUT_PREFIX" \
        -F "$REF_FA" \
        -T "$THREADS" \
        2>&1 | grep -E "scan|completed|FSJ|time" || true
    check_exists "SCAN2 FSJ counts ($SAMPLE_NAME)" "${OUT_PREFIX}.fsj_counts"

    echo -e "$SAM\t${OUT_PREFIX}.fsj_counts\t${SPLIT_NUM}\t${SAMPLE_NAME}" >> "$FINALIZE_TSV"
done

# --- Stage 2b: FINALIZE ---
info "--- Stage 2b: FINALIZE ---"
${JAVA_BIN} -cp "${CLASSPATH}" "${MAIN_CLASS}" FINALIZE \
    -I "$FINALIZE_TSV" \
    -F "$REF_FA" \
    -O "$FINALIZE_DIR/result" \
    -S 0 \
    2>&1 | grep -E "FINALIZE|Summary|Matrix|circRNA|time" || true
check_exists "Decoupled BSJ_Matrix" "$FINALIZE_DIR/result.BSJ_Matrix"
check_exists "Decoupled FSJ_Matrix" "$FINALIZE_DIR/result.FSJ_Matrix"

DECOUPLED_BSJ_LINES=$(tail -n +2 "$FINALIZE_DIR/result.BSJ_Matrix" | wc -l)
info "Decoupled pipeline detected $DECOUPLED_BSJ_LINES circRNAs."

# ---------------------------------------------------------------------------
# 5. Compare outputs
# ---------------------------------------------------------------------------
info "=== Comparing outputs ==="

# Reformat decoupled matrix columns to match original's sample column order
# (original uses basename of SAM; decoupled uses sampleName from FINALIZE TSV)
# Both should use the same name — original uses file basename from the path list.
# We'll compare after normalising column names to just the sample index.

normalise_matrix() {
    local f="$1"
    # Replace sample column headers with s1..sN (order by appearance), sort rows
    awk 'NR==1{
        printf "circRNA_ID"
        for(i=2;i<=NF;i++) printf "\ts%d", i-1
        printf "\n"
        next
    }
    { print | "sort" }' "$f"
}

info "Normalising and comparing BSJ matrices..."
A_BSJ=$(normalise_matrix "$ORIG_DIR/result.BSJ_Matrix")
B_BSJ=$(normalise_matrix "$FINALIZE_DIR/result.BSJ_Matrix")

if [[ "$A_BSJ" == "$B_BSJ" ]]; then
    ok "BSJ_Matrix: original and decoupled are IDENTICAL (${ORIG_BSJ_LINES} circRNAs)"
else
    fail "BSJ_Matrix: original and decoupled DIFFER"
    diff <(echo "$A_BSJ") <(echo "$B_BSJ") | head -30
fi

info "Normalising and comparing FSJ matrices..."
A_FSJ=$(normalise_matrix "$ORIG_DIR/result.FSJ_Matrix")
B_FSJ=$(normalise_matrix "$FINALIZE_DIR/result.FSJ_Matrix")

if [[ "$A_FSJ" == "$B_FSJ" ]]; then
    ok "FSJ_Matrix: original and decoupled are IDENTICAL (${ORIG_BSJ_LINES} circRNAs)"
else
    fail "FSJ_Matrix: original and decoupled DIFFER"
    diff <(echo "$A_FSJ") <(echo "$B_FSJ") | head -30
fi

# Also check universe coverage — every circRNA in original should be in universe
info "Checking universe coverage..."
ORIG_CIRCS=$(tail -n +2 "$ORIG_DIR/result.BSJ_Matrix" | awk '{print $1}' | sort)
UNIVERSE_CIRCS=$(grep "^chr" "$UNIVERSE_DIR/cohort.universe" | \
    awk '{printf "%s:%s|%s\n", $1, $2, $3}' | sort)
MISSING=$(comm -23 \
    <(echo "$ORIG_CIRCS") \
    <(echo "$UNIVERSE_CIRCS") | wc -l)
if [[ "$MISSING" -eq 0 ]]; then
    ok "Universe coverage: all original circRNAs are present in the universe"
else
    fail "Universe coverage: $MISSING circRNAs from original output are missing from universe"
    comm -23 <(echo "$ORIG_CIRCS") <(echo "$UNIVERSE_CIRCS") | head -10
fi

# ---------------------------------------------------------------------------
# 6. Checkpoint file format checks
# ---------------------------------------------------------------------------
info "=== Checkpoint file format checks ==="

# scan1_meta
META_SAMPLE="$SCAN1_DIR/sample1.scan1_meta"
for KEY in samFile readLen readNum fileSplitNum; do
    if grep -q "^${KEY}=" "$META_SAMPLE"; then
        ok "scan1_meta has field: $KEY"
    else
        fail "scan1_meta missing field: $KEY"
    fi
done

# universe file header
UNIVERSE_HEADER=$(head -1 "$UNIVERSE_DIR/cohort.universe")
if [[ "$UNIVERSE_HEADER" =~ ^seqLen=[0-9]+ ]]; then
    ok "Universe file header format: '$UNIVERSE_HEADER'"
else
    fail "Universe file header malformed: '$UNIVERSE_HEADER'"
fi

# fsj_counts format (4 tab-separated columns)
FSJ_COUNTS_FILE="$SCAN2_DIR/sample1.fsj_counts"
BAD_LINES=$(awk 'NF!=4' "$FSJ_COUNTS_FILE" | wc -l)
FSJ_TOTAL=$(wc -l < "$FSJ_COUNTS_FILE")
if [[ "$BAD_LINES" -eq 0 ]] && [[ "$FSJ_TOTAL" -gt 0 ]]; then
    ok "fsj_counts format: $FSJ_TOTAL lines, all 4 columns"
else
    fail "fsj_counts format: $BAD_LINES malformed lines out of $FSJ_TOTAL"
fi

# ---------------------------------------------------------------------------
# 7. Summary
# ---------------------------------------------------------------------------
echo ""
echo "========================================"
echo "  TEST SUMMARY"
echo "========================================"
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
