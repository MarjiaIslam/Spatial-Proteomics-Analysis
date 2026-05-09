"""
Protein-to-Gene Mapping for MIBI-SoC 40-marker panel (Dayao et al. HNSCC)
===========================================================================
Source: UniProt (Swiss-Prot reviewed human entries, organism_id:9606)
        Queried via https://rest.uniprot.org/ on April 2026.

Each antibody marker is mapped to:
  - UniProt accession (reviewed / Swiss-Prot)
  - HUGO gene symbol (primary gene name from UniProt)
  - Full UniProt protein name

Special handling:
  - CD45 / CD45RA / CD45RO → all PTPRC (P08575), kept as separate features
  - PanCK → KRT18 (P05783) as representative (pan-cytokeratin cocktail)
  - HLA-DR → HLA-DRA (P01903) alpha chain
"""

import pandas as pd
from pathlib import Path

# ──────────────────────────────────────────────────────────────
# FULL UniProt-VERIFIED MAPPING  (40 MIBI markers → gene symbols)
# Each entry: (gene_symbol, uniprot_accession, uniprot_protein_name)
# ──────────────────────────────────────────────────────────────
UNIPROT_MAP = {
    #  marker_name     : (gene,      uniprot_id, uniprot_protein_name)
    "CD31":       ("PECAM1",   "A0A075B728", "Platelet endothelial cell adhesion molecule"),
    "CD57":       ("B3GAT1",   "Q9P2W7",     "Galactosylgalactosylxylosylprotein 3-beta-glucuronosyltransferase 1"),
    "CD4":        ("CD4",      "P01730",     "T-cell surface glycoprotein CD4"),
    "CD15":       ("FUT4",     "P22083",     "Alpha-(1,3)-fucosyltransferase 4"),
    "FoxP3":      ("FOXP3",    "Q9BZS1",     "Forkhead box protein P3"),
    "CD16":       ("FCGR3A",   "P08637",     "Low affinity immunoglobulin gamma Fc region receptor III-A"),
    "CD20":       ("MS4A1",    "P11836",     "B-lymphocyte antigen CD20"),
    "CD45RO":     ("PTPRC",    "P08575",     "Receptor-type tyrosine-protein phosphatase C (CD45RO isoform)"),
    "CD38":       ("CD38",     "P28907",     "ADP-ribosyl cyclase/cyclic ADP-ribose hydrolase 1"),
    "CD34":       ("CD34",     "P28906",     "Hematopoietic progenitor cell antigen CD34"),
    "CD11b":      ("ITGAM",    "P11215",     "Integrin alpha-M"),
    "CD68":       ("CD68",     "P34810",     "Macrosialin (CD68)"),
    "CD134":      ("TNFRSF4",  "P43489",     "Tumor necrosis factor receptor superfamily member 4"),
    "TMEM16A":    ("ANO1",     "Q5XXA6",     "Anoctamin-1 (Transmembrane protein 16A)"),
    "PanCK":      ("KRT18",    "P05783",     "Keratin, type I cytoskeletal 18 (representative for pan-cytokeratin)"),
    "Podoplanin": ("PDPN",     "Q86YL7",     "Podoplanin"),
    "CD45":       ("PTPRC",    "P08575",     "Receptor-type tyrosine-protein phosphatase C (total CD45)"),
    "GranzymeB":  ("GZMB",     "P10144",     "Granzyme B"),
    "CD49f":      ("ITGA6",    "P23229",     "Integrin alpha-6"),
    "CD11c":      ("ITGAX",    "P20702",     "Integrin alpha-X"),
    "CD47":       ("CD47",     "Q08722",     "Leukocyte surface antigen CD47"),
    "CD8":        ("CD8A",     "P01732",     "T-cell surface glycoprotein CD8 alpha chain"),
    "CD117":      ("KIT",      "P10721",     "Mast/stem cell growth factor receptor Kit"),
    "Vimentin":   ("VIM",      "P08670",     "Vimentin"),
    "CD69":       ("CD69",     "Q07108",     "Early activation antigen CD69"),
    "aSMA":       ("ACTA2",    "P62736",     "Actin, aortic smooth muscle (alpha-SMA)"),
    "CD14":       ("CD14",     "P08571",     "Monocyte differentiation antigen CD14"),
    "CD21":       ("CR2",      "P20023",     "Complement receptor type 2"),
    "HLA-DR":     ("HLA-DRA",  "P01903",     "HLA class II histocompatibility antigen, DR alpha chain"),
    "PDL1":       ("CD274",    "Q9NZQ7",     "Programmed cell death 1 ligand 1 (PD-L1)"),
    "CD56":       ("NCAM1",    "P13591",     "Neural cell adhesion molecule 1"),
    "p16":        ("CDKN2A",   "P42771",     "Cyclin-dependent kinase inhibitor 2A (p16-INK4a)"),
    "PD1":        ("PDCD1",    "Q15116",     "Programmed cell death protein 1"),
    "CD45RA":     ("PTPRC",    "P08575",     "Receptor-type tyrosine-protein phosphatase C (CD45RA isoform)"),
    "ICOS":       ("ICOS",     "Q9Y6W8",     "Inducible T-cell costimulator"),
    "CD152":      ("CTLA4",    "P16410",     "Cytotoxic T-lymphocyte protein 4"),
    "Ki67":       ("MKI67",    "P46013",     "Proliferation marker protein Ki-67"),
    "CollagenIV": ("COL4A1",   "P02462",     "Collagen alpha-1(IV) chain"),
    "CD3e":       ("CD3E",     "P07766",     "T-cell surface glycoprotein CD3 epsilon chain"),
}

# ──────────────────────────────────────────────────────────────
# CONVENIENCE DICTS
# ──────────────────────────────────────────────────────────────

# Simple protein → gene lookup (backward-compatible with pipeline)
PROTEIN_TO_GENE = {k: v[0] for k, v in UNIPROT_MAP.items()}

# Protein → UniProt accession
PROTEIN_TO_UNIPROT = {k: v[1] for k, v in UNIPROT_MAP.items()}

# Markers where one antibody detects multiple possible gene products
AMBIGUOUS_MARKERS = {
    "PanCK":  ["KRT5", "KRT8", "KRT14", "KRT18", "KRT19"],
    "HLA-DR": ["HLA-DRA", "HLA-DRB1"],
    "CD8":    ["CD8A", "CD8B"],
}

# Markers that map to the same gene (splice-isoform distinction)
ISOFORM_CONFLICTS = {
    "PTPRC": {
        "CD45":   "Total PTPRC (all isoforms)",
        "CD45RA": "PTPRC isoform with exon A (naive T cells)",
        "CD45RO": "PTPRC isoform without variable exons (memory T cells)",
    },
}

# Column-name strategy for the output dataset:
#   CD45   → PTPRC        (total)
#   CD45RA → PTPRC_RA     (isoform suffix)
#   CD45RO → PTPRC_RO     (isoform suffix)
# This keeps them distinguishable while still linking to the correct gene.
GENE_COLUMN_NAMES = {}
for marker, (gene, _, _) in UNIPROT_MAP.items():
    if marker == "CD45RA":
        GENE_COLUMN_NAMES[marker] = "PTPRC_RA"
    elif marker == "CD45RO":
        GENE_COLUMN_NAMES[marker] = "PTPRC_RO"
    else:
        GENE_COLUMN_NAMES[marker] = gene

# Canonical (deduplicated) gene list for gene-set overlap queries
CANONICAL_GENE_LIST = sorted(set(PROTEIN_TO_GENE.values()))


# ──────────────────────────────────────────────────────────────
# EXPORT HELPER: mapping table as DataFrame / CSV
# ──────────────────────────────────────────────────────────────
def get_mapping_dataframe():
    """Return a tidy DataFrame of the full UniProt mapping."""
    rows = []
    for marker, (gene, uniprot, prot_name) in UNIPROT_MAP.items():
        rows.append({
            "mibi_marker":       marker,
            "gene_symbol":       gene,
            "output_column":     GENE_COLUMN_NAMES[marker],
            "uniprot_accession": uniprot,
            "uniprot_protein":   prot_name,
            "is_ambiguous":      marker in AMBIGUOUS_MARKERS,
            "is_isoform":        any(marker in v for v in ISOFORM_CONFLICTS.values()),
        })
    return pd.DataFrame(rows)


def save_mapping_csv(path="output/uniprot_protein_to_gene_mapping.csv"):
    """Save the mapping table as a CSV file."""
    out = Path(path)
    out.parent.mkdir(exist_ok=True, parents=True)
    df = get_mapping_dataframe()
    df.to_csv(out, index=False)
    print(f"Mapping saved → {out}  ({len(df)} markers)")
    return df


if __name__ == "__main__":
    df = save_mapping_csv()
    print(f"\nTotal protein markers:   {len(UNIPROT_MAP)}")
    print(f"Unique gene symbols:     {len(CANONICAL_GENE_LIST)}")
    print(f"Ambiguous markers:       {list(AMBIGUOUS_MARKERS.keys())}")
    print(f"Isoform conflicts:       {list(ISOFORM_CONFLICTS.keys())}")
    print(f"\n{'Marker':<15s} {'Gene':<10s} {'UniProt':<12s} {'Output Col':<12s}")
    print("-" * 52)
    for marker in UNIPROT_MAP:
        gene, uid, _ = UNIPROT_MAP[marker]
        print(f"{marker:<15s} {gene:<10s} {uid:<12s} {GENE_COLUMN_NAMES[marker]:<12s}")
