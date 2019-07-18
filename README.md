# GlypNir

An prototyped automated workflow for estimation of uqniue PSM N-glycan composition.

## Requirements
* Python >= 3.5
* Python packages: Biopython, Pandas, Tornado, xlrd

## Operation

`
python3 main.web-glycoform.py
`

The backend would be start up at port `9001`. The GUI is available at `http://localhost:9001/` if access
on the computer running the backend. If accessing from outside of the computer running the backend,
it should be accessible at `http://ip-address:9001` through the web browser.

1. The initial inputs are a collection of Excel files where filename include replicate label as well as experiment conditions.
1. In the auto filename parser, you can enter the identifier for each replicate, condition, and protein delimited by a comma.
These information will be parsed into the correct files. You can also directly modify these replicate, condition and protein labels on the GUI.
However condition must be the same for those replicates you want to group together as well as the proteins. These proteins label will have to be
the same as the master accession id in the Excel files and the id of the sequence in the fasta file that you will supply later.
1. In the reference sequences, a multi-fasta file should be supplied that contains all of the targeted proteins sequences.
1. An optional alignment files might be submitted for reference purpose when working with paralogs/homologs.
1. Minimal Area Intensity input decides the minimal cut-off threshold for PSM to be consider in the analysis.
1. Maximum N-linked Sites input decides the maximum cut-off N-linked sites count where PSMs that have more than the specified amount will be discarded.

