# ifo-fep: Markov blanket claims depend on the observation protocol

**English** | [繁體中文](README.zh-TW.md)

Papers, code and checks for two linked working papers by CK Hung and the Emmy Team (Amoy Technology), version 1.0.1, October 2026 (version 1.0 is in `papers/*/_old/`).

| Paper | What it is |
|---|---|
| [Interface With Whom? Markov Blanket Claims Depend on the Observation Protocol](papers/fep-blanket/20261022-fep-blanket-v1.0.1.md) ([PDF](papers/fep-blanket/20261022-fep-blanket-v1.0.1.pdf), [繁體中文](papers/fep-blanket/20261022-fep-blanket-v1.0.1-zh-tw.md)) | The main paper: the point of view, the main results, and what is not yet shown. |
| [Finite-Window Path Screening in Linear Diffusions: An Exact Condition and a Residual-Leak Formula](papers/path-blanket/20261022-path-blanket-v1.0.1.md) ([PDF](papers/path-blanket/20261022-path-blanket-v1.0.1.pdf), [繁體中文](papers/path-blanket/20261022-path-blanket-v1.0.1-zh-tw.md)) | The technical companion: proofs, numerical checks, and a disclosure of failed predictions and corrections (Appendix A). |

## In one paragraph

The free energy principle describes a Markov blanket as the interface between inside and outside. We build on that idea and ask a second question: an interface with whom? For linear stochastic systems we give exact conditions for when a blanket screens inside from outside under different observation protocols, a separation of screening from encoding, a single-file tool, and a real-data case that remains undecided. The point of view comes from a wider framework we are developing (IFO); every result holds without it. We hope experts will test it, correct it and take it further.

## Code

- `code/tools/blanket_report.py`: a single file that reports a proposed blanket for several observers. Run it without arguments for a demo. Needs NumPy and SciPy.
- `code/track-a/`: one script per check (M1–M17d and earlier). Each script states its prediction (or "post hoc") in its header, prints its result and writes one JSON file, which is included here. The companion paper's Reproducibility table maps scripts to sections.
- `code/tools/md2pdf.py`: turns a paper in Markdown into a PDF (Chrome and KaTeX).

Some script headers refer to internal planning notes (`plan/...`, `inbox/...`) that are not part of this repository; the prediction and the result in each header stand on their own.

Not included: the fMRI data (a public copy of Gilson et al., 2018, at github.com/mb-BCA/notebooks_review2019), pyMOU (github.com/mb-BCA/pyMOU, commit f0626be), and pyDMBD (github.com/bayesianempirimancer/pyDMBD, commit e058e9c). They are needed only for the fMRI and DMBD checks.

## Licence and citation

Papers: CC BY 4.0. Code: MIT. See `LICENSE` and `LICENSE-papers.md`. To cite, see `CITATION.cff`. DOI of version 1.0: [10.5281/zenodo.23524009](https://doi.org/10.5281/zenodo.23524009). DOI for all versions (always the latest): [10.5281/zenodo.23524008](https://doi.org/10.5281/zenodo.23524008). Record: https://zenodo.org/records/23524009

## Feedback

Corrections, counter-examples and pointers to prior art are very welcome: please open an issue.
