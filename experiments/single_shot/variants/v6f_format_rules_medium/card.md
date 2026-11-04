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

# Output consistency rules

- Only supported pairs get a relation type; explicit_negative, uncertain and not_found always have relation type none (`-`).
- A P record links a gene entity to a disease or clinical_event entity; symptom and phenotype entities are never P endpoints.
- With several supporting claims, the main relation type follows the display priority causal > correlation > association. The main source describes only the selected claim; use both only when background and this study support that same type.
- Each field of a record holds only its listed value; identity is r, p or u.
