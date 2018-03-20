# motif-mark
Tool for visualizing protein binding motifs related to splicing.

See plot.svg for example output.

`motif-mark.py` takes a fasta file in which each sequence must be of the form intronEXONintron. It also required a text file of desired motifs separated by lines.

The output is a .svg file showing each gene and the positions of any query motifs.

A jupyter-notebook file is used for development of the script and is included in the repo.
