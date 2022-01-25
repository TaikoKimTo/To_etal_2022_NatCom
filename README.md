# To_etal_2022_NatCom

# counting_mC_within_genes.pl is the Perl script for counting methylated C over total C within genes in Arabidopsis thaliana.
# use tair10_gene_model_startsort.txt for the reference information.
# The methylation input data (input1) needs to include "Chr position ori methylC totalC context".
# The context should be X (CG), Y(CHG), Z(CHH).
# The file should be sorted by Chr and position.
# example of INPUT1 : Chr1 109 + 17 19 X
# If you start with CX_report file from bismark, CX_report file can be converted by AWK in the command line as follows;
# awk '($1 ~"^Chr"){print $1,$2,$3,$4,$4+$5,($6=="CG"?"X":($6=="CHG"?"Y":"Z"))}' "CX_report file name" |sort -k 1,1 -k 2n > INPUT1.txt
