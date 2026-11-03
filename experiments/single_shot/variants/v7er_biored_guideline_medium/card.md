# Task

Extract gene–disease relations from the paper below.

1. Find all genes and all diseases mentioned in the paper.
2. For every gene × disease pair, decide whether the paper gives evidence of an association, correlation or causal relationship between them:
   - supported: the paper gives evidence of a relationship;
   - explicit_negative: the paper reports that there is no relationship;
   - uncertain: there is some evidence, but it is not conclusive;
   - not_found: no evidence.
3. For supported pairs, give the relation type (causal, correlation or association) and whether the evidence comes from this study's results or from background knowledge.
4. Cite the supporting sentences by ID.

## What counts as a gene
- Genes, proteins, mRNAs and other gene products, of any species.
- Gene families of paralogs or orthologs named as one unit (e.g., bone morphogenetic proteins, matrix metalloproteinases) are genes.
- Functional groupings are not genes (e.g., "tumor suppressor gene", "protein kinase"); neither are general words such as "gene" or "protein".
- If a gene and a disease share a name, the gene reading is a gene ("ALD gene") and the patient or disease reading is a disease ("ALD patient").

## What counts as a disease
- Diseases, symptoms and disease-related phenotypes are all diseases: give them the type `disease`, not `symptom` or `phenotype`.
- Broad disease terms such as cancer, tumor, pain and death are diseases. Very high-level terms such as "genetic disorder" or "familial disorder", and general words such as disease, syndrome, deficiency or complications on their own, are not.
- Mental states and phenotypes that are not disease-related are not diseases.

## When a pair is supported
- The paper explicitly states the relation anywhere in the text; the gene and the disease do not need to be in the same sentence.
- Count positive and negative correlations: more or less of the gene product causes, increases, prevents or treats the disease; knockdown or knockout of the gene changes the disease.
- Count genetic associations: a variant of the gene observed in patients with the disease, a variant that increases or decreases disease risk, or a disease caused by deficiency of the gene product. A confirmed variant–disease association also links the variant's gene to the disease.
- Minor associations count, e.g., an association observed in only a few patients.
- A disease that belongs to the target disease (e.g., "autosomal recessive disease" in a paper about a specific recessive disorder) is linked to the same gene.
- Do not link a gene as a prognostic factor of the disease (e.g., expression predicting survival) unless the paper states a relation with the disease itself.
- Do not link symptoms of the target disease to the gene unless the paper clearly relates the symptom to the gene.
- A statement that the gene is not related to the disease is explicit_negative. A statement too unclear to tell is uncertain.

## Source
- this study: the relation belongs to the main point of the paper (its results or conclusions).
- background: the relation is context, such as results of previous studies.
