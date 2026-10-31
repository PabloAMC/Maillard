# MFT routes in pentose + Cys (2026-10-31)

## 0. New finding: Liu 2023's "19 % fragments" is an artefact (P 0.8)
- Liu Table 3: MFT is 41/3/2/1/**12**/40, so M+4 / M+5 = 0.30.
- NIST EI [M−H]+/M+ (webbook JCAMP): MFT 0.28. Furfural 23/21 against 0.99; 2-methylfuran 20/29 against 0.58.
- The unlabelled : labelled imbalance (41 : 52) equals the unlabelled ion's uncounted M−1.
- So the claimed H-loss correction failed. Corrected recombination is ~6 %, the impurity level. No measured fragment rise with temperature; withdraw sec. 11's Liu signature.

## 1. Route mix (ribose / xylose + Cys, pH ~5)
**95-100 C, hours.**
- A (0.75): ≥ 85 % of MFT is intact C5, not via free NF. Evidence:
  - Cerny 2003 T3: fed NF gives 7 %. T2: non-(0, 5) channels 5 %, the impurity level.
  - Cerny 2004: "C2+C3 not relevant".
  - de Roos Table I: NF + Cys gives 48× less MFT than ribose.
  - The NF marker 2-mercapto-3-pentanone is absent from ribose.
- B (0.15): NF 15-40 %. C (0.05): fragments > 10 %. Other: 0.05.

**140-170 C, 20-60 min.**
- E (0.5): NF is 30-80 %, fragments ≤ 6 %. Evidence:
  - 2-mercapto-3-pentanone appears from ribose only at 140-145 C (Cerny 2003 p. 2715).
  - In-situ NF is 4.8 mol % (Hofmann T5).
  - Per mol Cys, NF/Cys gives 8.5× ribose/Cys (T10 against T1).
  - My estimate: in-situ NF could supply 30-100 % of T1's MFT.
- F (0.35): intact non-NF stays ≥ 70 %. T5 says "not exclusively via NF". T4 (ribose + H2S) lacks an amino acid, which shuts the DDP step (Nedvidek).
- G (0.10): fragments ≥ 15 %.

**Other.** Thiamine ~54 % of MFT at 145 C (Cerny 2007b). FFT: furfural + H2S at both temperatures (pot E 92 : 8; NF 0 %).

**The temperature flip is too big for one barrier.** NF goes from 48× below ribose to 8.5× above it, ~400×, which would need ΔEa ≈ 170 kJ/mol. Part is the 20-min window at 145 C (ribose route still lagging; 100 C / 6 h makes ~8× more MFT).

## 2. Why every arm fails Cerny but fits Schieberle
**The engine needs NF to make the rise.**
- Its only accumulating MFT precursor is free NF (DPO → NF 38-53 % against ~3 % measured, EX-018).
- DDP is a pass-through fed by DPO that is on from t = 0 (`r_ttca_deg`).
- `ch_mft_decay` destroys 0.39 of the 0.42 mM made, so net MFT ≈ rate × lifetime. Only an accelerating precursor gives ×35 from 30 to 360 min.
- Fed NF at 145 mM then dominates pot D, because no fitted row ≤ 100 C constrains the split.

**The data do not need an accelerating precursor.** MFT rises from 6 to 12 h, and labelled standards survive heating (CF-013): the lifetime is hours. With slow loss, constant DDP supply after a 15-20 min come-up gives (t−τ) ratios of 23-35. If DDP carries the MFT, pot D needs k_nf(95)/k_nf(145) ≤ ~1/30, i.e. Ea ≥ ~85 kJ/mol, inside the class.

**Most likely missing element: MFT's loss law.** It is too fast at 100 C, or it is co-reactant dependent (van Seeventer: non-oxidative "electrophilic coupling").
- Loss law: 0.40
- NF → MFT as a reduction-limited branch of an NF–sulfide intermediate (mercaptoketones 16× MFT in Whitfield 1999; engine 0.003-0.075): 0.25
- An intact-C5 supply gate (Feng 2022's fed ARP shows no lag): 0.20
- H2S supply or headspace: 0.15

**Dry check:** sec. 13 corner, `k_mft_decay` ×0.1-0.2, route rescaled; read 30 / 60 min, pot D.

## 3. Discriminating observations
1. **Spiked labelled MFT in Schieberle's pot (100 C, 1-12 h): the loss constant.** On disk: bounds only (CF-013; van Seeventer at 50 C; Guth 1995).
2. **Time-resolved pot D (fed NF + [13C5]ribose, 100 C, 0.25-6 h).** Gated ribose route: NF share high early, falling to 7 %; slow NF step: ≤ 10 % throughout. On disk: no. Liu: 168 C, no fed NF.
3. **DDP measured as its quinoxaline, rising before MFT.** On disk: never measured.
4. **2-mercapto-3-pentanone / MFT from ribose vs temperature.** On disk: partly (absent at 95 C, present at 145 C). A presence / absence signature now.
5. **H2S time course in a pot.** On disk: none.

## 4. Papers not on disk
DOIs checked on Crossref (title, authors, journal match); abstracts via Crossref / Semantic Scholar, else content **unverified**.
1. van den Ouweland & Peer 1975, *JAFC* 23:501. doi:10.1021/jf60199a045. Title: volatiles from norfuraneol + H2S. Abstract elided.
2. Hofmann & Schieberle 1996, "Studies on intermediates generating ... 2-methyl-3-furanthiol ...", *Flavour Science* 182-187. doi:10.1533/9781845698232.4.182. No abstract.
3. Li, Schieberle & Steinhaus 2017, *Eur Food Res Technol* 243:69. doi:10.1007/s00217-016-2723-z. Title: SIDA for hydrogen sulphide, applied to durian.
4. Huber & Ledl 1990, *Carbohydr Res* 204:215. doi:10.1016/0008-6215(90)84037-U. Title: 1-amino-1,4-dideoxy-2,3-hexodiuloses in the Maillard reaction. An amine-bound 1,4-dideoxy route (hexose); unverified.
5. Nomi et al. 2021, *BBB* 85:2042. doi:10.1093/bbb/zbab119. Abstract: pentosone, 1-DP, 3-DP and MGO form mainly from Amadori degradation (carnosine + xylose, pH 5).
6. Cerny 2008, *Ann NY Acad Sci* 1126:66. doi:10.1196/annals.1433.011. Abstract: a review of labelled-precursor model reactions; low odds it holds the promised norfuraneol follow-up.
7. Hofmann dissertation, TU München, ~1995: the source of Schieberle's Table IV. No DOI, **unverified**.

**Searched, none found** (P ~0.1): time-resolved MFT labelling, fed DDP, an H2S time course in a sugar + Cys pot.
