# Helpful Bioinformatics Scripts #

## Convert two-column csv to FASTA ##
- In bash run the following:

      awk '{ printf ">%s\n%s\n",$1,$2 }' testfile.txt > outputfile.txt