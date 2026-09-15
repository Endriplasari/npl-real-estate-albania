# Depositing the package and obtaining a DOI

A GitHub repository alone is not sufficient for a journal data availability
statement: it has no DOI, no guarantee of permanence, and its contents can be
changed after publication. The standard solution is to keep the working copy on
GitHub and mint a citable, frozen DOI through Zenodo.

## Steps

1. Create the GitHub repository and push this directory.
2. Sign in to Zenodo with the GitHub account and authorise the integration.
3. In Zenodo, under GitHub, switch the repository toggle to "on".
4. In GitHub, create a release (for example `v1.0.0`). Zenodo archives the
   release automatically and mints a DOI.
5. Copy the concept DOI, which always resolves to the latest version, into the
   data availability statement of the manuscript.

`CITATION.cff` in the repository root supplies the metadata Zenodo reads.

## Data availability statement

Suggested wording, to be inserted in the manuscript once the DOI exists:

> All data used in this study are publicly available from the Bank of Albania
> and the Institute of Statistics of Albania. The assembled dataset, the code
> that constructs it from the primary sources, and the code implementing every
> test, simulation and calibration reported in this paper are archived at
> [DOI]. The semi-annual real estate market surveys from which the house price
> index was reconstructed are published by the Bank of Albania and are not
> redistributed; the figures extracted from each issue are documented in the
> repository. No confidential or supervisory data were used.

## Alternatives

If institutional policy requires a different repository, the Open Science
Framework and Harvard Dataverse both mint DOIs and both accept this structure
unchanged. Zenodo is preferred here only because of the automatic GitHub
integration.
