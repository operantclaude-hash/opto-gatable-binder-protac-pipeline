#!/bin/bash
set -e
cd /home/thinkingscopeanalysis/als_therapeutics/nanobody_optogenetic_protac/opto_scaffold_pipeline/campaigns/constitutive_top10/degron_comparison

YAML_DIR=boltz2_yamls
OUT_DIR=boltz2_results

mkdir -p "$OUT_DIR"

# Run all 8 predictions sequentially (GPU can't handle simultaneous)
YAMLS=(
    degron_nterm_dark_target
    degron_nterm_lit_target
    degron_cterm_dark_target
    degron_cterm_lit_target
    degron_nterm_dark_fold
    degron_nterm_lit_fold
    degron_cterm_dark_fold
    degron_cterm_lit_fold
)

TOTAL=${#YAMLS[@]}
COUNT=0

for NAME in "${YAMLS[@]}"; do
    COUNT=$((COUNT + 1))
    echo "[$COUNT/$TOTAL] Running $NAME..."
    if [ -d "$OUT_DIR/$NAME" ] && ls "$OUT_DIR/$NAME"/predictions/*/confidence_*.json >/dev/null 2>&1; then
        echo "  Already completed, skipping."
        continue
    fi
    boltz predict "$YAML_DIR/${NAME}.yaml" \
        --out_dir "$OUT_DIR/$NAME" \
        --diffusion_samples 3 \
        --recycling_steps 3 \
        --model boltz2 \
        --override 2>&1 | tail -5
    echo "  Done."
done

echo ""
echo "All $TOTAL predictions complete!"
