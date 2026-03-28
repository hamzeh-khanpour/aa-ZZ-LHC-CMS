# aa-ZZ-LHC-CMS

A small analysis framework for studying exclusive photon-fusion diboson production,

$\gamma\gamma \to ZZ \to \mu^+\mu^- \nu\bar\nu,$

with an emphasis on

- anomalous quartic gauge couplings (aQGCs) in an EFT description,
- the onset of unitarity-violating growth in the high-energy tail,
- truth-level clipping of events above a chosen energy scale, and
- the impact of clipping on kinematic distributions such as $m_{ZZ}$ and $p_T(\mu^+\mu^-)$.

This repository is intended as a lightweight, transparent study environment for testing EFT-sensitive observables and for comparing uncut and clipped signal shapes.

---

## Physics motivation

In the Standard Model, neutral quartic gauge interactions are highly constrained, and exclusive high-mass photon-fusion diboson production is a clean place to search for deviations from the SM prediction. The CMS+PPS analysis note for exclusive $\gamma\gamma\to WW$ and $\gamma\gamma\to ZZ$ with intact protons focuses precisely on this kind of topology and on the high-mass tail accessible with proton tagging. In that note, the two-photon collision energy can be determined event by event from the proton kinematics, and the signal interpretation is performed in terms of anomalous quartic gauge couplings. See **CMS AN-19-211** for the analysis context, PPS acceptance, and signal modeling choices. 

At a more formal level, EFT descriptions of genuine quartic gauge couplings are known to produce amplitudes that grow with energy. The paper **“Unitarity Constraints on Anomalous Quartic Couplings”** studies these effects systematically using coupled-channel partial-wave unitarity, including all relevant helicities and both $J=0$ and $J=1$ partial waves. In a linear realization of electroweak symmetry, the lowest genuine QGC operators arise at dimension eight.  

---

## EFT picture

The basic EFT logic is simple:

$$\mathcal{L}_{\text{eff}} = \mathcal{L}_{\text{SM}} + \sum_i \frac{f_i}{\Lambda^4}\, O_i$$

for the genuine dimension-8 quartic operators in a linear EFT description. In practice, the repository is used to study one benchmark coefficient at a time, with all other anomalous couplings set to zero.

Two important points are worth keeping in mind:

1. **The EFT is an expansion.** It is trustworthy only below some characteristic energy scale.
2. **The anomalous amplitudes grow with energy.** If one pushes the EFT too far into the ultraviolet, the predicted cross section can become dominated by events in a region where perturbative unitarity is no longer reliable. 

The CMS note uses a LEP-style anomalous-coupling parameterization for its central result, and later discusses the mapping to dimension-8 $f_{M,i}$ operators. This repository instead uses a modern UFO-based sample generation workflow and truth-level observables from the generated event record, but the physical issue is the same: the far high-mass tail can be EFT-sensitive and eventually unphysical if no validity criterion is imposed. 

---

## Why unitarity becomes a problem

If anomalous quartic couplings are large enough, the partial-wave amplitudes violate the perturbative unitarity bound at sufficiently high energy. The general condition can be expressed in terms of the partial-wave amplitudes $T^J$, with the strongest constraints obtained after diagonalizing the coupled-channel matrix in particle and helicity space. The dedicated unitarity study by Almeida, Éboli, and Gonzalez-Garcia explicitly shows that one must account for all coupled channels and that $J=1$ can matter in addition to $J=0$ in realistic multidimensional operator scans. 

The CMS analysis note also discusses this issue explicitly. In its $\gamma\gamma\to VV$ study, it computes the energy at which the expected anomalous signal would violate unitarity and then evaluates the impact of removing the unphysical region. In the note, the unitarity-violating part of the anomalous signal is removed from the analysis while the data and backgrounds are left unchanged.  

---

## What “clipping” means in this repository

In this repository, **clipping** means:

> keep only events whose **truth-level hard scale** lies below a chosen clipping scale $\Lambda_{\text{clip}}$.

For the exclusive Born-level process

$\gamma\gamma \to ZZ,$

the natural hard scale is the diboson invariant mass:

$m_{ZZ}^{\text{truth}} \simeq \sqrt{\hat s} \simeq \sqrt{s_{\gamma\gamma}}.$

So the practical clipping prescription used here is

$m_{ZZ}^{\text{truth}} < \Lambda_{\text{clip}}.$

This is conceptually the same kind of operation used in the CMS note: the clipping threshold is an energy scale motivated by unitarity, while the actual implementation removes anomalous signal events **at generator/truth level** above that scale. The reconstructed distributions are then recomputed from the surviving events.  

Because this repository works directly at LHE truth level, the clipping is especially clean: there is no detector smearing, no proton misreconstruction, and no migration from truth below the cut to reconstructed masses above the cut. In a full detector analysis, a residual reconstructed tail can survive above the nominal clipping scale due to resolution and combinatorics; that is one of the effects visible in the CMS note. 

---

## Why clipping should be applied

Clipping is not a cosmetic choice. It addresses a real physics problem.

Without a validity criterion, the anomalous signal can become dominated by the far UV region. This makes the resulting limits deceptively strong, because they are being driven by a part of phase space where the EFT may not be trustworthy. The CMS note explicitly compares unclipped and clipped limits and shows that clipping weakens the expected sensitivity, as it should once the unphysical high-energy growth is removed. It also emphasizes that the effect of clipping depends on the channel and on the acceptance.  

In this repository, clipping is used for the same reason:

- to remove the truth-level events above a chosen EFT-validity scale,
- to inspect how much of the signal rate lives in the problematic tail,
- to compare full and clipped shapes, and
- to provide a more conservative interpretation of the anomalous benchmark.

---

## Sample generation strategy

The CMS PPS note uses signal MC samples for $\gamma\gamma\to WW$ and $\gamma\gamma\to ZZ$ produced with **FPMC**, including the outgoing protons, and then combines the central detector simulation with a dedicated PPS proton fast simulation. It also states that these samples were produced with the **dimension-6 AQGC formalism**, with **no form factors and no unitarization**, which is exactly why the clipping discussion becomes necessary in that note. 

In this repository, the workflow is different but closely analogous in spirit:

1. generate $\gamma\gamma\to ZZ$ hard events with a UFO model,
2. decay the bosons with MadSpin,
3. reconstruct truth-level observables from the LHE record,
4. study the distributions before and after truth-level clipping.

A typical generation pattern is:

```text
import model <your_UFO_model>
generate a a > z z NP=1
output <process_dir>
launch
```

with the decays applied using MadSpin, for example

```text
define vl = ve vm vt
define vl~ = ve~ vm~ vt~
decay z > mu+ mu- @1
decay z > vl vl~ @1
launch
```

The `@1` syntax is the MadSpin way to enforce the mixed exclusive topology $ZZ\to(\mu^+\mu^-)(\nu\bar\nu)$.

---

## Analysis observables used here

The main truth-level observables reconstructed in this project are:

- **$m_{ZZ}$**: from the two status-2 $Z$ bosons in the LHE record,
- **$p_T(\mu^+\mu^-)$**: from the final-state muon pair,
- **$\sqrt{\hat s}$**: from the incoming partons/photons in the LHE record.

For an exclusive $2\to2$ Born-level sample, these variables are tightly related. In particular,

$p_T(\mu^+\mu^-) = p_T(Z),$

and the kinematic boundary is

$p_T^{\max} = \frac{1}{2}\sqrt{m_{ZZ}^2 - 4m_Z^2}.$

This is why the hardest $p_T$ tail is strongly correlated with the largest $m_{ZZ}$ values.

---

## Clipping workflow implemented in the scripts

The plotting scripts in this repository keep the full-sample plots and then produce clipped comparisons and clipped-only distributions.

The logic is:

1. read the post-decay LHE file,
2. reconstruct event-level truth observables,
3. choose a clipping threshold $\Lambda_{\text{clip}}$,
4. keep only events with
   $m_{ZZ}^{\text{truth}} < \Lambda_{\text{clip}},$
5. rebuild the distributions for:
   - full sample,
   - clipped sample,
   - overlay comparison,
   - clipped-only shapes.

The current scripts also write summary files and produce survival/fraction curves such as

- fraction of events with $m_{ZZ} > \Lambda_{\text{cut}}$,
- fraction of events with $\sqrt{\hat s} > \Lambda_{\text{cut}}$.

That makes it easy to quantify how much of the sample lives above a chosen EFT-validity scale.

---

## Interpreting $\Lambda_{\text{clip}}$

A key point is that **$\Lambda_{\text{clip}}$** is the **threshold**, not the variable.

In the CMS note, the clipping scale is derived from a unitarity argument and then applied to the generated anomalous signal. In this repository, the practical implementation is

$m_{ZZ}^{\text{truth}} < \Lambda_{\text{clip}},$

because for the exclusive truth-level process studied here, $m_{ZZ}^{\text{truth}}$ is the cleanest event-by-event proxy for the hard scale.

So, for example,

```bash
python3.10 plot_dimuon_pt_mZZ_EFTchecks_truthClippingLabels_from_lhe.py \
  /path/to/decayed_events.lhe \
  --mzz-cut 1400 \
  --truncate-on mzz
```

means:

> apply **truth-level clipping** at $\Lambda_{\text{clip}} = 1.4$ TeV by requiring $m_{ZZ}^{\text{truth}} < 1.4$ TeV.

---

## Why the clipped distributions matter

The unclipped and clipped plots answer slightly different physics questions.

### Unclipped distributions
These show the raw EFT-enhanced signal prediction and are useful to identify where the anomalous contribution becomes large.

### Clipped distributions
These show the same signal after removing the events above the chosen validity scale. They are the distributions that should be used for a conservative EFT interpretation.

If the clipped distributions change dramatically compared to the full sample, that is a warning sign that the analysis sensitivity is being driven by the far UV tail.

---

## Repository contents

Typical contents of the workflow include:

- truth-level LHE parsing scripts,
- MadGraph / MadSpin generation cards,
- clipping and EFT-check plotting scripts,
- overlay plots for full vs clipped distributions,
- clipped-only distributions for CMS-style truth-level clipping studies.

Representative outputs include:

- `pt_dimuon_normalized.pdf`
- `pt_dimuon_dsigma_dpT.pdf`
- `mZZ_normalized.pdf`
- `mZZ_dsigma_dmZZ.pdf`
- `mZZ_vs_pt_dimuon_2D_normalized.pdf`
- `profile_mZZ_vs_pt_dimuon.pdf`
- `fraction_mZZ_gt_Lambdacut.pdf`
- `pt_dimuon_normalized_truncated_overlay.pdf`
- `mZZ_normalized_truncated_overlay.pdf`
- `pt_dimuon_normalized_clipped_only.pdf`
- `mZZ_normalized_clipped_only.pdf`

---

## Important caveat

The CMS PPS note used a LEP-style dim-6 parameterization for its central anomalous-coupling study and later discussed the mapping to dimension-8 $f_{M,i}$ operators. This repository may use a different UFO model and different benchmark conventions, including dimension-8 coefficients such as $f_{M,0}/\Lambda^4$.

So the **numerical value** of the clipping scale is **benchmark dependent**.

The reusable physics logic is this:

1. determine or estimate a unitarity-safe scale for the benchmark of interest,
2. treat that scale as $\Lambda_{\text{clip}}$,
3. apply a truth-level cut on the hard event scale,
4. compare unclipped and clipped results.

The repository therefore follows the **CMS clipping philosophy**, while adapting the implementation to the exclusive truth-level $\gamma\gamma\to ZZ$ setup used here.

---

## References

1. **CMS AN-19-211**: *Exclusive WW and ZZ production in the fully hadronic channel with protons reconstructed in PPS*. Provides the PPS-based exclusive diboson analysis context, the signal MC setup with FPMC, the absence of built-in unitarization in those samples, and the clipping-based treatment of unitarity-violating signal tails. 

2. **Eduardo da Silva Almeida, O. J. P. Éboli, M. C. Gonzalez-Garcia**, *Unitarity Constraints on Anomalous Quartic Couplings*, arXiv:2004.05174. Provides the coupled-channel partial-wave unitarity framework with all relevant helicities and both $J=0$ and $J=1$, and is a useful formal reference when deciding how conservative a clipping scale should be in a dimension-8 EFT study. 

---

## Suggested citation for this repository

If you use this repository as part of a phenomenology study, please cite the relevant CMS PPS analysis note and the unitarity paper above, and describe clearly which anomalous-coupling benchmark and which truth-level clipping prescription were used.
