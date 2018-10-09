# prot_align

Here you can find the python script I did as part of my MSc in Bioinformatics
at UAB.

The script aligns two protein sequences from a Multi-FASTA file using a
simplified version of the needleman-wunsch algorithm and reports the score and
the alignment in a printable format.
The algorithm uses a dynamic programming approach considering that each
insertion has a penalty of -4, each deletion has a penalty of -2, and each
substitution has a different cost depending on the amino acid substituted using
the specified BLOSUM matrix.

Running the following code in a Windows terminal
```
py -2.7 align.py --input .\\data\\GHRs.fasta --block_size 90
```
will return
```
Alignment score is: 1165

GHR-II  MA--AAL-------------T--LLFCLYILTSSALE--SA--SEQV--H-----PQRDPHLTGCVSANMETFRCRWNVGTLQNLSKPGE  62
        ||                  |  ||  |  | || |   |   |  |  |     |   || | | |   ||||| |  |   ||| ||
GHR-I   MAVFSSSSSSSSSSSSSSSSTSNLLL-L-LLVSS-LDWLSTRGSVFVMDHMTSSAPV-GPHFTECISREQETFRCWWSPGGFHNLSSPGA  86

GHR-II  LRLFYINKLSPLDPPKEWTECPHYS-IDRPNECFFNKNHTSVWTPYKVQLRSRDESTLY-DEN-TFTVDAIVQPDPPVDLTWTTLNESLS  149
        || ||  | ||     || ||| ||   |  ||||  |||||| ||  |||     | | ||   |||  || ||||| | || || | |

GHR-I   LRVFYLKKDSP-N--SEWKECPEYSHLKR--ECFFDVNHTSVWIPYCMQLRGQNNVT-YLDEDYCFTVENIVRPDPPVSLNWTLLNISPS  170

GHR-II  GTYYDIILSWKPPQSADVAMGWMTLQYEVQY--RSASSDLWHAVE--PVTVTQRSLFGLKHNVNHEVRVRCKMLAG-KEFGEFSDSIF--  232
        |  ||    | || ||||  |||   || ||  |      | | |  |   ||    ||      ||  || | |  | ||||||| |
GHR-I   GLSYDVMVNWEPPPSADVGAGWMRIEYEIQYTERN-TTN-WEALEMQP-H-TQQTIYGLQIGKEYEVHIRCRMQAFVK-FGEFSDSVFIQ  255

GHR-II  V-HIPAKVSSFPV-VALLLFGALCLVAILMLVI-ISQQEKLMFILLPPVPGPKIRGIDPELLKKGKLRELTSIL-GGPPN-LR---PELY  314
        |  ||   | ||   ||  || |    || | | ||||  || ||||||| ||| |||||||||||| ||  || ||    |    |  |
GHR-I   VTEIPSQDSNFPFKLALI-FGVLGIL-ILILLIGISQQPRLMMILLPPVPAPKIKGIDPELLKKGKLDELNFILSGGGMGGLSTYAPDFY  343

GHR-II  NNDPWVEFIDLDIEE-QS-DKLTDL--DTDCLMH--RSLSS--N--CTPVSIGFRDDDSGRASCCDPDLPSDPEASPFHP-LIPNQTLSK  393
           ||||||  | |      |      ||  |       |   |  |      | ||||||||| ||||  |         | | |
GHR-I   QDEPWVEFIEVDAEDADAAEKEENQGSDTQRLLDPPQPVSHHMNTGCAN-AVSFPDDDSGRASCYDPDL-HDQDTLMLMATLLPGQPEDG  431

GHR-II  EVS--C-QTAS--EPSS-P-VQSPASGEPPFAALGREAMYTQVSEVRSSGKVLLSPEEQTEV-EKTTG-KD-TEKDIM-AE--KEKAKKE  470
        | |      |   | |  | ||    | |    |     | ||| |  || | |||  |    | |    |   |     |   ||  ||
GHR-I   EDSFDVVERAPVIERSERPLVQTQTGG-PQ-TWLNTD-FYAQVSNVMPSGGVVLSPGQQLRFQESTSAAEDEAQKKGKGSEDSEEKTQKE  518

GHR-II  --FQLLVVNADHG-GYTSELNAGKMS-PRLSIGDQSEPGLTG-D-LSPLP-PASPYHESDTTAVSP--LP--P-----APV--YTVVEGV  542
          ||||||    | ||| | ||   | |  |      ||  |     | |    |         ||  ||  |     |||  ||||  |
GHR-I   LQFQLLVVDPE-GSGYTTESNARQISTPP-ST-PM--PG-SGYQTIHPQPVETKPAATAENNQ-SPYILPDSPQSQFFAPVADYTVVQEV  601

GHR-II  DRQNSLLLTP--NSTPAPQLII-P-KTMPT-PGGYLTPDLLGSITP  583
        | | |||| |     | | |   | |     | || ||||||   |
GHR-I   DSQHSLLLNPPPRQSPPPCLPHHPTKALAAMPVGYVTPDLLGNLSP  647
```
