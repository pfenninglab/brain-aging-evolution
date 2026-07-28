## An Epigenetic Signature of Vulnerable Neurons is Under Selective Pressure Associated with Longevity Across Placental Mammals.

A multi-scale framework to investigate cell-type-specific aging programs and their contribution to cellular vulnerability in neurodegenerative disease.

## Overview

Age is the primary risk factor for neurodegenerative diseases, which are characterized by cell-type-specific vulnerability. Identifying the biological mechanisms underlying brain aging remains challenging because of the complex set of interacting, age-associated biochemical pathways acting across a great diversity of neural cell types and neuron subtypes.

This repository provides the code for a framework that dissects cell-type- and cell-state-specific aging gene regulatory programs and their contribution to cellular vulnerability, by integrating epigenomics, AI methodology, and the natural diversity of lifespan across placental mammals.

This framework systematically integrates evolutionary, tissue-level, and cell-intrinsic dimensions to resolve how regulatory programs of aging emerge and contribute to disease vulnerability.

The framework connects two scales of aging:

1. **Cross-species to within-species.** It links cross-species lifespan measurements to aging within the human population.
2. **Intrinsic vs. systemic.** It decomposes intrinsic factors (cell-autonomous effects) from systemic factors that more broadly influence cells across the entire tissue.

   ![Figure 1: Overview of the multi-scale framework](figures/brain_aging_evolution_framework.pdf)

## Framework components

The approach integrates these scales through seven components:

1. Resolving neuronal and non-neuronal cortical cell-type-specific longevity-associated regulatory elements (LARs) across 240 mammals.
2. Developing a transcriptional brain aging clock to decompose systemic aging from intrinsic aging in cortical neurons within the human population.
3. Resolving cell-type-specific intrinsic aging.
4. Linking intrinsic age to cellular vulnerability and aging hallmarks.
5. Identifying differentially accessible regulatory elements in vulnerable cell types.
6. Integrating cross-species and human-focused analyses to identify shared regulatory elements connecting mammalian longevity with human aging at the cell-type level — vulnerability–longevity-associated regulatory elements (VLRs).
7. Investigating the evolutionary conservation of enhancer–gene synteny for VLRs across 577 vertebrate species.

## Repository structure

```
.
├── processing/        # data processing and preparation for downstream analysis
├── analysis/          # analysis code, organized by framework component
└── model/             # transcriptional brain aging clock (systemic and intrinsic aging clock)
```
