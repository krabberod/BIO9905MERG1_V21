Input Data
===================

query_sequences.fasta
reference_alignment.fasta
reference_alignment.phylip
reference_tree.newick



PaPaRa - aligning
===================

command:
papara -t reference_tree.newick -s reference_alignment.phylip -q query_sequences.fasta -r

output:
papara_alignment.default
papara_log.default
papara_quality.default



epa-ng - placement
===================

command:
epa-ng --split reference_alignment.fasta papara_alignment.default

output:
query.fasta
reference.fasta

command:
raxml-ng --evaluate --msa reference.fasta --tree reference_tree.newick --model GTR+G

output:
reference.fasta.raxml.bestModel
reference.fasta.raxml.bestTree
reference.fasta.raxml.log
reference.fasta.raxml.rba
reference.fasta.raxml.startTree

command:
epa-ng -t reference_tree.newick -s reference.fasta -q query.fasta --model reference.fasta.raxml.bestModel

output:
epa_result.jplace



gappa - downstream analyses
===================

HEAT TREE
command:
gappa examine heat-tree --jplace-path epa_result.jplace --mass-norm absolute --write-svg-tree --write-newick-tree --write-nexus-tree

output:
tree.svg
tree.nexus
tree.newick

TAX. ASSIGNMENT
command:
gappa examine assign --jplace-path epa_result.jplace --taxon-file clades.txt --krona

output:
assign_krona.profile
assign_labelled_tree.newick
assign_profile.tsv

EXTRACTION
command:
gappa prepare extract --jplace-path epa_result.jplace --clade-list-file clades.txt --fasta-path query.fasta --color-tree-file extract_tree --samples-out-dir samples --sequences-out-dir sequences

output: 
extract.sh.txt
extract_tree.svg
results.log
samples
sequences

TREE WITH QUERY SEQS.
command:
gappa examine graft --jplace-path epa_result.jplace --fully-resolve --out-dir labelled_tree


LWR HISTOGRAM
command:
gappa examine lwr --jplace-path epa_result.jplace

output:
lwr_histogram.csv
lwr_list.csv

EDPL HISTOGRAM
command:
gappa examine edpl --jplace-path epa_result.jplace

output:
edpl_histogram.csv
edpl_list.csv



ZIP
=================

command:
zip -r results.zip *
