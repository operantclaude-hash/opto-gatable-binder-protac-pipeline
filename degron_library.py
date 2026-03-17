"""
Degron system components for the opto-gatable binder scaffold pipeline.

Contains: VVD dimerizer, E3 ubiquitin ligase sequences (TRIM21, NEDD4L),
flexible GS linkers, and AAV packaging constants.

Supports two degradation pathways:
  - Proteasomal: TRIM21_RBCC (K48 ubiquitin chains -> 26S proteasome)
    For cytosolic proteins (e.g., p300 HAT)
  - Lysosomal:   NEDD4L_HECT (K63 ubiquitin chains -> endosomal/lysosomal)
    For membrane proteins

References:
  - VVD: Zoltowski et al., PNAS 2008 (light-induced dimerization)
  - TRIM21: Clift et al., Nature Communications 2023 (bioPROTAC)
  - NEDD4L: Columbia DiVa, Nature Communications 2025 (membrane targets)
"""

# ===== VVD Light-Gated Dimerizer (Neurospora crassa) =====

VVD_SEQ = (
    "MAHTLYAPGGYDIMGWLIQIMNRPNPQVELGPVDTSCALVLCDLKQKDTPVVYASEAFLY"
    "MTGYSNAEVILGRNCRFLQSPDGMVKPKSTRKYVDSNTINTMRKAIDRNAEVQVEVVNFK"
    "KNGQRFVNFLTMIPVRDETGEYRYSMGFQCETECFKLG"
)  # 158aa
VVD_LENGTH = 158

# ===== E3 Ubiquitin Ligase Library =====

E3_LIGASE_LIBRARY = {
    "TRIM21_RBCC": {
        "seq": (
            "MASAARLTMMWEEVTCPICLDPFVEPVSIECGHSFCQECISQVGKGGGSVCPVCRQRFLK"
            "NLRPNRQLANMVNNLKEISQEAREGTQGERCAVHGERLHLFCEKDGKALCWVCAQSRKHR"
            "DHAMVPLEEAAQEYQEKLQVALGELRRKQELAEKLEVEIAIKRADWKKTVETQKSRIHAE"
            "FVQQKNFLVEEEQRQLQELEKDEREQLRILGEKEAKLAQQSQALQELISELDRRCH"
        ),  # 236aa, Human TRIM21 (UniProt P19474, aa 1-232)
        "length": 236,
        "pathway": "proteasomal",
        "ubiquitin_chain": "K48",
        "source": "Human (UniProt P19474, aa 1-232)",
        "description": "RING-Bbox-Coiled-coil E3 ligase domain",
        "domains": "RING(16-55) + B-box(92-123) + Coiled-coil(128-236)",
    },
    "NEDD4L_HECT": {
        "seq": (
            "MKVLAAHPNNHTVTVQPPPTLAGEPVQLPDGPGEEPPVPFPFKKELFLGKKFHKAEGFDS"
            "AVSVSEPVVLPVADLSESEEEKTWMKLYCKEEADKLVRNFNFDDLGRDFLLERAKYEYVM"
            "DSVFPPHYHRSYISVKNLDGPMTLHHGSDVLFTCDCDPEYEGDSPCIYRLDIVKRDAQTQ"
            "ELSLWNDLLAMRPDHEWEQHKFTCFNCMMLWAGEKTYNMPCVGRCLTLPVYRADQNTVED"
            "FWRMMWEQQATTIIGSVMVFPMQPTLQRREKLILQSFQIYPPSQFDKNFHLTFENEIDNG"
            "QWCLEIDYPHMRGKPTWCERYDIITEPYLAHYLLDGNKTVECYLPPFVPYSMLNDHWEPT"
            "CSQITAFDFRLNHMGTV"
        ),  # 377aa, Human NEDD4-2 (UniProt Q96PU5, HECT domain)
        "length": 377,
        "pathway": "lysosomal",
        "ubiquitin_chain": "K63",
        "source": "Human (UniProt Q96PU5, aa 594-975)",
        "description": "HECT-type E3 ligase domain for membrane proteins",
    },
}

# ===== Degradation Pathway -> Default E3 Ligase =====

PATHWAY_DEFAULTS = {
    "proteasomal": "TRIM21_RBCC",
    "lysosomal": "NEDD4L_HECT",
}

# ===== Flexible GS Linkers (for connecting VVD/E3 to opto-fusion) =====

GS_LINKER_LIBRARY = {
    "GS5": {"seq": "GGGGS", "length": 5},
    "GS10": {"seq": "GGGGSGGGGS", "length": 10},
    "GS15": {"seq": "GGGGSGGGGSGGGGS", "length": 15},
    "GS20": {"seq": "GGGGSGGGGSGGGGSGGGGS", "length": 20},
}

# ===== AAV Packaging Constraints =====

AAV_MAX_BP = 4700
AAV_MAX_AA = AAV_MAX_BP // 3  # ~1566 amino acids
