#!/bin/bash
set -e
cd "$(dirname "$0")"

YAML_DIR=boltz2_yamls
OUT_DIR=boltz2_results

mkdir -p "$OUT_DIR"

# Dark folds (single chain) + dark+HAT + dark+catcore for top 4 binders
YAMLS=(
    pep4s2_50aa_E2pppE2_dark
    pep4s2_50aa_dark_HAT
    pep4s2_50aa_dark_catcore
    pep13s3_20aa_E2pppE2_dark
    pep13s3_20aa_dark_HAT
    pep13s3_20aa_dark_catcore
    pep30s5_20aa_E2pppE2_dark
    pep30s5_20aa_dark_HAT
    pep30s5_20aa_dark_catcore
    nb62s3_120aa_E2pppE2_dark
    nb62s3_120aa_dark_HAT
    nb62s3_120aa_dark_catcore
)

TOTAL=${#YAMLS[@]}
COUNT=0

for NAME in "${YAMLS[@]}"; do
    COUNT=$((COUNT + 1))
    echo "[$COUNT/$TOTAL] Running $NAME..."
    if [ -d "$OUT_DIR/$NAME" ] && ls "$OUT_DIR/$NAME"/*/predictions/*/confidence_*.json >/dev/null 2>&1; then
        echo "  Already completed, skipping."
        continue
    fi
    boltz predict "$YAML_DIR/${NAME}.yaml" \
        --out_dir "$OUT_DIR/$NAME" \
        --diffusion_samples 3 \
        --recycling_steps 3 \
        --model boltz2 \
        --override 2>&1 | tail -3
    echo "  Done."
done

echo ""
echo "All $TOTAL predictions complete!"
