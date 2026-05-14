```
SPATIAL TRAJECTORY ANALYSIS PIPELINE - SIMPLIFIED FLOW
═══════════════════════════════════════════════════════════════

┌──────────────┐
│  Tissue Img  │  Multiplexed imaging
│  (40+ marks) │  (CODEX/IMC/MIBI)
└──────┬───────┘
       │
       ▼
┌──────────────────┐
│ Cell Segmentation│  Detect individual cells proteomics data
│  (Cellpose)      │  → Cell masks
└──────┬───────────┘
       │
       ▼
┌──────────────────┐
│  Phenotyping     │  Cluster cells by protein expression
│  (Leiden)        │  → Cell type labels
└──────┬───────────┘
       │
       ▼
┌──────────────────┐
│ Graph Construction│  Build spatial network
│  (knn+Delaunay)      │  → G = (V, E)
└──────┬───────────┘
       │
       ▼
┌──────────────────────────┐
│ Feature Engineering      │  Enrich graph with:
│                          │  • Density • Entropy
│                          │  • Centrality • Gradients
└──────┬───────────────────┘
       │
       ▼
┌──────────────────┐
       ▼
┌──────────────────┐
│ Patient Aggregate│  Mean, max, distribution
│                 │  → Patient metrics
└──────┬───────────┘
       │
       ▼
┌──────────────────┐
│ Survival         │  Cox regression
│ Analysis         │  → HR, C-index, KM curves
└──────────────────┘

OUTCOME: Disease progression scores predictive of survival
```





```
FEATURE ENGINEERING STEP - COMPREHENSIVE FLOWCHART
═══════════════════════════════════════════════════════════════

INPUT: Basic Graph G=(V,E) + Expression Matrix + Coordinates + Clusters
                                    │
                    ┌───────────────┼───────────────┐
                    │               │               │
                    ▼               ▼               ▼
            ┌──────────────┐ ┌──────────────┐ ┌──────────────┐
            │ Expression   │ │ Coordinates  │ │ Cluster      │
            │ (N × 40)     │ │ (N × 2)      │ │ Labels       │
            └──────┬───────┘ └──────┬───────┘ └──────┬───────┘
                   │                │                │
                   └────────────────┬────────────────┘
                                    │
        ┌───────────────────────────┼───────────────────────────┐
        │                           │                           │
        ▼                           ▼                           ▼
   ╔═════════════════════╗    ╔═════════════════════╗    ╔═════════════════════╗
   ║  NODE FEATURES      ║    ║  NODE FEATURES      ║    ║  NODE FEATURES      ║
   ║  (Topology)         ║    ║  (Cell Context)     ║    ║  (Neighborhood)     ║
   ╠═════════════════════╣    ╠═════════════════════╣    ╠═════════════════════╣
   ║ 1. Local Density    ║    ║ 1. Degree Central   ║    ║ 1. Entropy H(i)     ║
   ║    ρᵢ = |N(i)|/πr² ║    ║    |N(i)|/(N-1)     ║    ║    -Σpₖlog(pₖ)      ║
   ║    (crowding)       ║    ║    (hub importance) ║    ║    (cell diversity)  ║
   ║                     ║    ║                     ║    ║                     ║
   ║ 2. Expression       ║    ║ 2. Closeness Cent   ║    ║ 2. Boundary Score   ║
   ║    Gradient         ║    ║    (N-1)/Σd(i,j)   ║    ║    B(i) = diff_nbr/  ║
   ║    std(neighbor     ║    ║    (network center) ║    ║    |N(i)| (interface)║
   ║    expression)      ║    ║                     ║    ║                     ║
   ║    (progression)    ║    ║ 3. Betweenness     ║    ║ 3. [Space for more] ║
   ║                     ║    ║    Σσₛₜ(i)/σₛₜ      ║    ║                     ║
   ║                     ║    ║    (bridge cells)   ║    ║                     ║
   ╚═════════════════════╝    ╚═════════════════════╝    ╚═════════════════════╝
        │                           │                           │
        └───────────────────────────┼───────────────────────────┘
                                    │
                    ┌───────────────┴───────────────┐
                    │                               │
                    ▼                               ▼
        ┌─────────────────────────┐      ┌─────────────────────────┐
        │  CONCATENATE NODE       │      │  EDGE FEATURES          │
        │  FEATURES               │      │                         │
        │  (N × F_enriched)       │      ╠═════════════════════════╣
        │  ═══════════════════    │      ║ 1. Expression Similarity║
        │  [Expression (40) +     │      ║    sim(i,j) = cosine   ║
        │   Density (1) +         │      ║    xᵢ·xⱼ/(||xᵢ||||xⱼ||)║
        │   Entropy (1) +         │      ║                         ║
        │   Boundary (1) +        │      ║ 2. Interaction Type    ║
        │   Degree Cent (1) +     │      ║    Homotypic (0)       ║
        │   Closeness Cent (1) +  │      ║    Tumor-Immune (1)    ║
        │   Betweenness Cent (1)+║      ║    Tumor-Stroma (2)    ║
        │   Gradient (1)]         │      ║    Immune-Immune (3)   ║
        │                         │      ║                         ║
        │  Total: N × ~48 features│      ║ 3. Distance Weight     ║
        └─────────┬───────────────┘      ║    w(i,j)=exp(-d²/2σ²)║
                  │                       ║    (spatial decay)      ║
                  │                       ╚═════════════════════════╝
                  │                               │
                  └───────────────────┬───────────┘
                                      │
                                      ▼
                        ┌──────────────────────────┐
                        │  ENRICHED GRAPH G'       │
                        │  (V', E')                │
                        ├──────────────────────────┤
                        │ • Node features: N × 48  │
                        │ • Edge features: E × 3   │
                        │ • Richer context for     │
                            │   downstream models      │
                        └──────────────────────────┘
```

SUMMARY TABLE - NODE FEATURES (Total: 8 Feature Types)
═══════════════════════════════════════════════════════════════

Feature                  │ Dimension │ Range    │ Biology
─────────────────────────┼───────────┼──────────┼─────────────────────
1. Local Density         │ 1         │ [0, 1]   │ Tissue architecture
2. Neighborhood Entropy  │ 1         │ [0, 1]   │ Cell type mixing
3. Boundary Score        │ 1         │ [0, 1]   │ Interface membership
4. Degree Centrality     │ 1         │ [0, 1]   │ Hub importance
5. Closeness Centrality  │ 1         │ [0, 1]   │ Network centrality
6. Betweenness Centrality│ 1         │ [0, 1]   │ Bridge importance
7. Expression Gradient   │ 1         │ [0, ∞)   │ Local heterogeneity
8. Raw Expression       │ 40        │ [0, 1]   │ Marker levels


SUMMARY TABLE - EDGE FEATURES (Total: 3 Feature Types)
═══════════════════════════════════════════════════════════════

Feature              │ Dimension │ Range   │ Biology
─────────────────────┼───────────┼─────────┼──────────────────────
1. Similarity        │ 1         │ [-1, 1] │ Expression alignment
2. Interaction Type  │ 1 (coded) │ 0-4     │ Cell-cell relationship
3. Distance Weight   │ 1         │ [0, 1]  │ Spatial decay


```
