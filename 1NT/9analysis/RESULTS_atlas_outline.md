# Pan-cancer interaction rewiring atlas

## Section title options

1. A pan-cancer atlas of recurrent interaction rewiring in matched tumor-normal samples
2. Recurrent interaction rewiring reveals coordinated cancer programs across tissues
3. A module-centric atlas uncovers organizing principles of cancer interaction rewiring

## Results structure

### 1. Paired tumor-normal analysis defines a pan-cancer interaction rewiring landscape

We first analyzed matched tumor-normal samples from 11 cancer types to quantify link-level perturbation for each interaction. For each cancer type, we summarized tumor-normal differences across paired samples and classified links as gained, lost, or unchanged in tumors relative to matched normal tissues. This yielded a pan-cancer interaction perturbation matrix capturing the direction and magnitude of rewiring across cancers.

Recommended figure panels:

- Workflow schematic: paired samples -> perturbation score -> recurrent edges -> modules
- Pan-cancer gain/loss burden
- Recurrent edge heatmap
- Recurrence distribution

Key files:

- `atlas_all.csv`
- `2recurrent_heatmap.pdf`
- `2recurrence_distribution.pdf`
- `2pan_cancer_burden.pdf`

Suggested text anchor:

Across 11 cancer types, recurrent perturbation patterns were readily apparent at the interaction level, indicating that tumor-associated rewiring is not random but follows shared pan-cancer principles.

### 2. Recurrently gained and lost edges reveal a conserved cancer rewiring backbone

To identify robust pan-cancer events, we aggregated interaction changes across cancer types and defined recurrently gained or lost edges based on direction-consistent recurrence. Using this framework, we identified 2,836 recurrently perturbed links, including 1,991 recurrently gained edges and 845 recurrently lost edges. The most recurrent gained edges were dominated by cell-cycle-associated interactions, whereas recurrently lost edges were enriched for tissue homeostasis, metabolic, and microenvironmental communication programs.

Recommended figure panels:

- Recurrent edge burden distribution
- Top recurrently gained/lost edge examples
- Cancer membership of recurrent links

Key files:

- `pan_cancer_recurrence_summary.csv`
- `universal_recurrent_links.csv`
- `recurrently_gained_links.csv`
- `recurrently_lost_links.csv`
- `3atlas_recurrence_distribution.pdf`
- `3atlas_recurrent_burden_by_cancer.pdf`

Suggested text anchor:

These recurrent edges define a conserved interaction rewiring backbone spanning multiple tumor types, with gained and lost rewiring occupying markedly different biological spaces.

### 3. Recurrent rewiring organizes into coordinated gained and lost modules

Rather than existing as isolated pairwise events, recurrently perturbed edges assembled into higher-order rewiring programs. Community detection on the recurrent edge graph identified 27 gained modules and 68 lost modules. The largest gained module was centered on mitotic regulators including CDC20, CDK1, BUB1, and BUB1B, whereas the largest lost module captured metabolic and differentiation-associated interactions involving ADH and ALDH family members.

Recommended figure panels:

- Module size-strength landscape
- Top gained and lost modules
- Module network view

Key files:

- `module_summary_all.csv`
- `gain_module_network.pdf`
- `loss_module_network.pdf`
- `3atlas_module_size_strength.pdf`
- `3atlas_top_modules.pdf`

Suggested text anchor:

This module-centric view shows that cancer rewiring is structured into coherent systems-level programs rather than dispersed edge-level abnormalities.

### 4. Functional annotation distinguishes proliferative gain modules from homeostatic loss modules

We next annotated each rewiring module using Hallmark, KEGG, and Gene Ontology enrichment. In total, 46 modules received robust functional annotation. Gained modules were dominated by chromosome segregation, DNA replication, ubiquitin-mediated protein modification, pyrimidine metabolism, and extracellular matrix remodeling. In contrast, lost modules were dominated by retinol and xenobiotic metabolism, chemokine and cytokine communication, glutathione metabolism, fatty acid catabolism, vascular smooth muscle contraction, and other tissue homeostasis programs.

Recommended figure panels:

- Module annotation landscape
- Top annotated modules
- Annotation source balance

Key files:

- `module_annotation_table.csv`
- `module_enrichment_all_terms.csv`
- `4module_annotation_landscape.pdf`
- `4top_module_annotations.pdf`

Suggested text anchor:

Together, these results indicate that tumors recurrently gain interaction programs linked to proliferation and matrix remodeling while losing programs associated with differentiated tissue function and microenvironmental communication.

### 5. Conserved and cancer-biased rewiring programs provide a framework for biological prioritization

The atlas supports a distinction between broadly conserved rewiring programs and modules with stronger cancer-type bias. Conserved gain modules were dominated by mitotic and replication machinery, consistent with core oncogenic proliferative demands shared across tissues. By contrast, many loss modules reflected context-dependent erosion of tissue-specific metabolic, stromal, and signaling programs, providing a route to classify recurrent rewiring into shared versus lineage-biased components.

Recommended next analyses:

- Quantify module conservation across cancers
- Label modules as conserved or cancer-biased
- Rank modules by recurrence entropy or cancer coverage

Suggested text anchor:

This separation between conserved and cancer-biased rewiring programs provides an organizing principle for prioritizing mechanistically informative modules.

### 6. A deep example can highlight the added value of interaction-level analysis

To demonstrate the biological utility of the atlas, one immune-related module should be selected as a deep example. An antigen-presentation or chemokine communication module would be especially suitable because it connects rewiring, tumor-immune context, and clinical interpretation. The key point of this section should be that module-level interaction rewiring captures structure that is not fully recoverable from gene-level analysis alone.

Recommended deep-example components:

- Pan-cancer module perturbation pattern
- Link-level versus gene-level comparison
- Within-module coordination
- Association with immune infiltration or outcome

## Figure storyboard

### Figure 1. Pan-cancer interaction rewiring landscape

- Paired tumor-normal workflow
- Gain/loss burden across cancers
- Recurrent edge heatmap
- Recurrence distribution

### Figure 2. Rewiring module atlas

- Recurrent edge network -> modules
- Top gained and lost modules
- Module size-strength landscape
- Module annotation summary

### Figure 3. Deep example of an immune-related rewiring module

- Module perturbation heatmap
- Module network
- Gene-level versus link-level comparison
- Immune infiltration or survival association

## Current quantitative anchors

- Cancer types analyzed: 11
- Recurrent links: 2,836
- Recurrently gained links: 1,991
- Recurrently lost links: 845
- Gained modules: 27
- Lost modules: 68
- Functionally annotated modules: 46

## Writing guidance

- Keep "interaction rewiring" as the central phrase.
- Treat edges as the discovery layer and modules as the biological interpretation layer.
- Avoid framing the atlas as a passive database.
- Emphasize organizing principles: conserved proliferation gain, tissue-homeostasis loss, and immune/stromal communication rewiring.
