# Field checks for PREREG A6 (STEP8_FIELDS_SPEC.md sec. 2)

Two readers, blind to each other and to every file under `results/`, check the facts the generator cannot derive
from the keyed text. One file per source and reader: `<source>.reader_a.yml`, `<source>.reader_b.yml`; the resolver
writes `<source>.yml`. Arm labels are the `label` strings of the resolved contrast file
`../contrasts/<source>.yml`, copied exactly; one entry per distinct label (arms that share a label share the facts).

```yaml
source: xu2007
reader: a            # a | b | resolved
pdf: data/articles/xu2007.pdf
fields: [F1, F2]     # which of F1 / F2 / F3 this source was checked for
ph_control: buffered # F1: buffered | adjusted_unbuffered | unadjusted | unstated (for every pH contrast, unless overridden)
ph_control_overrides: {}   # contrast id: value, where a contrast differs from the source-level value
arms:
  "pH 5.6 (final 5.5)":     # an arm label, exactly as in the contrast file
    final_ph: 5.5            # F2: the pH after heating as printed; null if not printed
    final_ph_location: "Table 1"   # where; add ", digitised" if read off a figure
    present: [cysteine, ribose]    # F3: every reactant in the pot (amine, sulfur, carbonyl, lipid, other), as printed
notes: ""            # one line: anything that could not be checked, and why
```

- **F1 definitions.** `buffered`: a named buffer at a printed concentration sets the pH. `adjusted_unbuffered`: the
  pH was set by titration (acid or base) in water or a solution with no added buffer. `unadjusted`: nobody set the
  pH (it is whatever the mixture gives). `unstated`: the paper does not say.
- **F2.** Only a final (after-heating) pH the paper prints or plots; never an initial or nominal one.
- **F3.** All reactants in the arm's pot, including a protein, peptide, hydrolysate or Amadori product if that is
  what the arm holds; not the solvent or the buffer salt.
