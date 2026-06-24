# Plan for Figure 2

##Logical flow
- We wanted to know how ECFs relate to each other
- They have a nested, hierarchical structure. Most of them overlap heavily, but there is a clear tail: ECFs largely contained within larger ones (Fig 2A)
- We thus next wanted to know if this tail represensts a genuine signal emerging from the ECF comparison
- So we compared the fragment ECFs (high qcov, low scov) with the highly overlapping ECFs (high qcov and scov) and asked if the former connected more diverse ECFs.
- We found that the fragment ECFs were more likely to connect unrelated ECFs as well as a greater number of ECF families (Fig 2B)
- Therefore, these results suggest that some ECFs act as bridging nodes between distant parts of the ECF space
- We next investigated how ECFs related to well defined domains, using ECOD and Pfam.
- As expected, ECFs spanned an entire spectrum, from subdomains, to superdomains (Fig 2C).
- They rarely formed a perfect equivalence to ECODs and Pfams, and many of them were novel (Fig 2C).
- That novelty remained regardless of the significance threshold, altough the proportion was affected (some could be distantly related).
- We next compared ECFs to ECODs directly and asked what types of domains they linked.
- We found that ECFs linked non-homologous ECODs within each class (Fig 2D).
- Two examples are shown (Fig 2E).

## List of Figures

- `p.qs.node`: per-node median qcov/scov diagonal (Figure 2A)
- `p.cov.heatmap`: raw pairwise qcov/scov heatmap (Suppl); an alternative to 2A
- `p.asym.dist`: qcov-scov asymmetry (Suppl or GONE); shows that the distribution is skewed to the right (more fragments)
- `p.asym.len`: asymmetry vs length; shows that fragments tend to be shorter; pretty obvious (Suppl or GONE)
- `pB`: Proportion of unrelated neighbours vs scov (Fig 2B)
- `p.cont.cont` full distribution of the neighbour's coverage; some are zero which is interesting (Suppl)

