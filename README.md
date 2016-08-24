## What is BackToDNA?
BackToDNA is written in Go language. The script converts an amino acid alignment back to its corresponding nucleotide alignment by processing and analyzing the BLAST results obtained from blasting these amino acid sequences against their original nucleotide sequences. This means that BackToDNA requires both the amino acid alignment and the DNA database containing the original sequences.

## How does it work?
BackToDNA reads BLAST results obtained from blasting amino acid sequences against their original nucleotide sequences and selects hits that have 100% identity and the highest bit score.
Matching regions in both amino acid sequences and its top nucleotide hits are recorded to locate the positions that needed to be reverse translated for each sequences in the amino acid alignment. Specifically, for a sequence that was translated from a positive frame, every amino acid in the matching regions will be replaced by its corresponding codon from the matching region of the amino acid sequence's top hit. If a sequence was translated from a negative frame, the matching nucleotide segment will be reverse complemented before its every codon can be used to replace its corresponding amino acid in the matching region.
Since the amino acid alignment comes with dash characters, each of them will be replaced with three dash characters so that BackToDNA's final output will remain aligned. For the same purpose, every amino acid that is not in the matching region will be replaced by three dash characters.
For amino acid sequences that have multiple hits with 100% identity and highest bit score, redundancy will be removed.

## Motivation
 Converting each amino acid sequence in an alignment back to its corresponding DNA nucleotide sequence without changing the integrity of the amino acid alignment will result in an alignment of all nucleotide sequences containing the gene of interest. With this nucleotide alignment, one can look at conserved gene regions and determine whether existing degenerative primers used to find this gene are matching well with many of these sequences. Using the aligned conserved region, it is also possible to design a new primer that better detect the presence of this gene in selected metagenomes.

## How to use BackToDNA?
1) Install GO

2a) If you have one nucleotide database and one blast result, in one directory:

Save your blast result as : database.out

Your nucleotide database: database.fa

Your amino acid alignment as: query.fa

2b) If you have multiple database files that results in multiple blast outputs, in one directory:

Save each blast result corresponding to its nucleotide file name. ie: If nucleotide file is nucleotide1.fa , save your result from blasting against this database as nucleotide1.out

Save your amino sequence alignment query as: query.fa

3) Usage: If the directory where you save your amino acid alignment, nucleotide database and blast out is 'yourdirectory', then command to run is:

./backtodna -outdir $yourdirectory -query $yourdirectory/query.fa > backtodna.fa

## Test
## License
