"""
Shared constants for the opto-gatable binder scaffold pipeline.

Contains: LOV2 sequences, cap sequences, linker library with geometric
properties, and physical constants for sweep volume analysis.
"""

# ===== LOV2 Sequences (Avena sativa LOV2) =====

LOV2_DARK = (  # 136aa, dark state (Jalpha helix packed against core)
    "LERIEKNFVITDPRLPDNPIIFASDSFLQLTEYSREEILGRNCRFLQGPETDRATVRKIR"
    "DAIDNQTEVTVQLINYTKSGKKFWNLFHLQPMRDQKGDVQYFIGVQLDGTEHVRDAAER"
    "EGVMLIKKTAENIDEAA"
)

LOV2_LIT = (  # 121aa, lit state (Jalpha released/unfolded)
    "LERIEKNFVITDPRLPDNPIIFASDSFLQLTEYSREEILGRNCRFLQGPETDRATVRKIR"
    "DAIDNQTEVTVQLINYTKSGKKFWNLFHLQPMRDQKGDVQYFIGVQLDGTEHVRDAAER"
    "EG"
)

# ===== Default Cap =====

UBQ_SEQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"  # 76aa

# ===== Linker Library =====

EAAAK = "EAAAK"

LINKER_LIBRARY = {
    # Kinked rigid linkers
    "E2ppE2": {
        "seq": EAAAK * 2 + "PP" + EAAAK * 2,
        "length": 22, "type": "kinked",
        "pre_helix": 10, "kink_res": 2, "post_helix": 10,
        "description": "~30 deg proline kink",
    },
    "E2pppE2": {
        "seq": EAAAK * 2 + "PPP" + EAAAK * 2,
        "length": 23, "type": "kinked",
        "pre_helix": 10, "kink_res": 3, "post_helix": 10,
        "description": "Extended proline kink (TOP GATING CANDIDATE)",
    },
    "E2dpgnE2": {
        "seq": EAAAK * 2 + "DPGN" + EAAAK * 2,
        "length": 24, "type": "kinked",
        "pre_helix": 10, "kink_res": 4, "post_helix": 10,
        "description": "~90-120 deg Type I' beta turn",
    },
    "E2npgdE2": {
        "seq": EAAAK * 2 + "NPGD" + EAAAK * 2,
        "length": 24, "type": "kinked",
        "pre_helix": 10, "kink_res": 4, "post_helix": 10,
        "description": "~90-120 deg Type I beta turn",
    },
    "E2gpngE2": {
        "seq": EAAAK * 2 + "GPNG" + EAAAK * 2,
        "length": 24, "type": "kinked",
        "pre_helix": 10, "kink_res": 4, "post_helix": 10,
        "description": "~120 deg Type II beta turn",
    },
    "E3dpgnE2": {
        "seq": EAAAK * 3 + "DPGN" + EAAAK * 2,
        "length": 29, "type": "kinked",
        "pre_helix": 15, "kink_res": 4, "post_helix": 10,
        "description": "Long pre-arm + DPGN turn",
    },
    "E2dpgnE3": {
        "seq": EAAAK * 2 + "DPGN" + EAAAK * 3,
        "length": 29, "type": "kinked",
        "pre_helix": 10, "kink_res": 4, "post_helix": 15,
        "description": "DPGN turn + long post-arm",
    },
    "E2ppE3": {
        "seq": EAAAK * 2 + "PP" + EAAAK * 3,
        "length": 27, "type": "kinked",
        "pre_helix": 10, "kink_res": 2, "post_helix": 15,
        "description": "PP kink + long post-arm",
    },
    "EP4": {
        "seq": (EAAAK + "P") * 4,
        "length": 24, "type": "kinked",
        "pre_helix": 5, "kink_res": 1, "post_helix": 5,
        "multi_kink": True,
        "description": "Pro-capped anti-bundle (4 kinks)",
    },
    "P16": {
        "seq": "P" * 16,
        "length": 16, "type": "ppII",
        "description": "Pure polyproline II helix",
    },
    # Straight controls
    "EAAAK3": {
        "seq": EAAAK * 3,
        "length": 15, "type": "straight",
        "description": "Short straight alpha-helical linker",
    },
    "EAAAK5": {
        "seq": EAAAK * 5,
        "length": 25, "type": "straight",
        "description": "Standard alpha-helical linker (baseline)",
    },
}

# ===== Physical Constants for Sweep Analysis =====

HELIX_RISE_PER_RESIDUE = 1.5   # Angstroms along helix axis (alpha helix)
PPII_RISE_PER_RESIDUE = 3.12   # Angstroms for polyproline II
DEFAULT_CAP_RADIUS = 14.0      # UBQ Rg in Angstroms
DEFAULT_DARK_CONE = 25.0       # degrees, constrained by Jalpha
DEFAULT_LIT_CONE = 60.0        # degrees, free after Jalpha release
STERIC_EXCLUSION_RADIUS = 3.0  # Angstroms
SWEEP_GRID_SPACING = 2.0       # Angstroms for voxelization
N_AZIMUTHAL_SAMPLES = 36
N_POLAR_SAMPLES = 10

# ===== Power Analysis Constants (from literature) =====

# ipTM → experimental binding probability (from BindCraft, BoltzGen, RFdiffusion data)
IPTM_SUCCESS_RATES = {
    (0.9, 1.0): 0.70,   # 60-80% → use 70%
    (0.8, 0.9): 0.40,   # 30-50% → use 40%
    (0.7, 0.8): 0.22,   # 15-30% → use 22%
    (0.5, 0.7): 0.10,   # 5-15%  → use 10%
    (0.0, 0.5): 0.03,   # <5%    → use 3%
}

# Gating ratio → probability of measurable optogenetic switching
GATING_SUCCESS_RATES = {
    (10.0, float("inf")): 0.80,
    (5.0, 10.0): 0.65,
    (2.0, 5.0): 0.40,
    (1.0, 2.0): 0.15,
    (0.0, 1.0): 0.05,
}
