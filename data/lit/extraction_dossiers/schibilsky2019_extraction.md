# Schibilsky 2019 — DOSSIER STUB (file stem `schibilsky2019`)

**Thesis not on disk; this dossier records only the bundle's keyed values and their provenance; source re-read
pending (download requested).** Nothing below was read from the primary source by the author of this stub
(2026-11-06). It exists so that a register entry about these bundles can name its source honestly instead of
borrowing another paper's dossier.

| field | value |
|---|---|
| Work | Schibilsky, S. "Bildung fluechtiger Verbindungen und Carbonsaeuren in Modellsystemen unter Maillard-Bedingungen." Dissertation, Technische Universitaet Berlin, 2019 (as the bundles cite it) |
| DOI | none registered (the bundles' `source_doi_absent_reason`: CrossRef returns no record, checked 2026-08-27) |
| Identifier | DepositOnce bitstream 969d7075-4e3d-48ec-9a26-fed34b3c7a55 (as the bundles cite it) |
| On disk | **No.** Not in the article folder on 2026-11-06 (M0 integrity batch 1, MF-28 / CK-06 found it absent) |
| Status | **STUB: second-hand.** Re-read pending; the download was requested by the M0 integrity batch-1 records |

## What the repository holds, and where it came from

The two hold-out bundles `mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019` and `..._pH80_Schibilsky2019`
(data/benchmarks/external_validation/maillard_path/) key a D-glucose + L-alanine model, 0.25 M each (0.5 M stock
diluted 1:2), unbuffered (pH set to 5.00 / 8.00 at t = 0), 130 C, 2 h. Their `content_verification` block (Wave U,
2026-08-27) says the thesis was then downloaded from DepositOnce and its title page, sec. 5.1.6 / 5.1.6.2, Tabelle 10
and Tabelle 17 were read; that PDF is not on disk now, so those quotations are **second-hand for every later reader**,
including this stub.

Keyed values (Tabelle 17 as quoted in the bundles, uM converted with Tabelle 10's molar masses; the arithmetic
re-checked by `m0_integrity_batch1_other_lanes` MF-28 and `m0_integrity_batch1_strecker_sulfur` CK-06, not the cells):

| compound | pH 5 | pH 8 | keyed (ppb, pH 5 / pH 8) |
|---|---|---|---|
| Furaneol (DMHF) | 9 uM | 46 uM | 1153.17 / 5893.98 (x 128.13) |
| 5-HMF | 454 uM | 798 uM | 57253.94 / 100635.78 (x 126.11) |
| Furfural | 17 uM | 25 uM | 1633.36 / 2402.0 (x 96.08) |

## Used by

- CF-133 (the acrylamide lane carries none of the trunk's pH terms, so these pots' DMHF and HMF are pH-flat on the
  engine). The Schibilsky pots are SPENT for deciding that question (their glycine stand-in was compared with the
  measured ratios at diagnosis).

## To do when the thesis is on disk

Re-read Tabelle 17 (the six cells above), Tabelle 10 (molar masses), sec. 5.1.6 / 5.1.6.2 (charge, pH setting,
vessel: 10 mL sealed ampoules with 2 mL or 20 mL headspace vials with 6 mL), then replace this stub with a full
dossier.
