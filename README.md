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
py -2.7 align.py --input .\\data\\GHRs.fasta --block_size 65
```
Will return
```
Alignment score is: 1165

GHR-II  MA--AAL-------------T--LLFCLYILTSSALE--SA--SEQV--H-----PQRDPHLTGC  37
        ||                  |  ||  |  | || |   |   |  |  |     |   || | |
GHR-I   MAVFSSSSSSSSSSSSSSSSTSNLLL-L-LLVSS-LDWLSTRGSVFVMDHMTSSAPV-GPHFTEC  61

GHR-II  VSANMETFRCRWNVGTLQNLSKPGELRLFYINKLSPLDPPKEWTECPHYS-IDRPNECFFNKNHT  101
         |   ||||| |  |   ||| || || ||  | ||     || ||| ||   |  ||||  |||
GHR-I   ISREQETFRCWWSPGGFHNLSSPGALRVFYLKKDSP-N--SEWKECPEYSHLKR--ECFFDVNHT  121

GHR-II  SVWTPYKVQLRSRDESTLY-DEN-TFTVDAIVQPDPPVDLTWTTLNESLSGTYYDIILSWKPPQS  164
        ||| ||  |||     | | ||   |||  || ||||| | || || | ||  ||    | || |
GHR-I   SVWIPYCMQLRGQNNVT-YLDEDYCFTVENIVRPDPPVSLNWTLLNISPSGLSYDVMVNWEPPPS  185

GHR-II  ADVAMGWMTLQYEVQY--RSASSDLWHAVE--PVTVTQRSLFGLKHNVNHEVRVRCKMLAG-KEF  224
        |||  |||   || ||  |      | | |  |   ||    ||      ||  || | |  | |
GHR-I   ADVGAGWMRIEYEIQYTERN-TTN-WEALEMQP-H-TQQTIYGLQIGKEYEVHIRCRMQAFVK-F  245

GHR-II  GEFSDSIF--V-HIPAKVSSFPV-VALLLFGALCLVAILMLVI-ISQQEKLMFILLPPVPGPKIR  284
        |||||| |  |  ||   | ||   ||  || |    || | | ||||  || ||||||| |||

GHR-I   GEFSDSVFIQVTEIPSQDSNFPFKLALI-FGVLGIL-ILILLIGISQQPRLMMILLPPVPAPKIK  308

GHR-II  GIDPELLKKGKLRELTSIL-GGPPN-LR---PELYNNDPWVEFIDLDIEE-QS-DKLTDL--DTD  340
        |||||||||||| ||  || ||    |    |  |   ||||||  | |      |      ||
GHR-I   GIDPELLKKGKLDELNFILSGGGMGGLSTYAPDFYQDEPWVEFIEVDAEDADAAEKEENQGSDTQ  373

GHR-II  CLMH--RSLSS--N--CTPVSIGFRDDDSGRASCCDPDLPSDPEASPFHP-LIPNQTLSKEVS--  396
         |       |   |  |      | ||||||||| ||||  |         | | |    | |
GHR-I   RLLDPPQPVSHHMNTGCAN-AVSFPDDDSGRASCYDPDL-HDQDTLMLMATLLPGQPEDGEDSFD  436

GHR-II  C-QTAS--EPSS-P-VQSPASGEPPFAALGREAMYTQVSEVRSSGKVLLSPEEQTEV-EKTTG-K  454
            |   | |  | ||    | |    |     | ||| |  || | |||  |    | |
GHR-I   VVERAPVIERSERPLVQTQTGG-PQ-TWLNTD-FYAQVSNVMPSGGVVLSPGQQLRFQESTSAAE  498

GHR-II  D-TEKDIM-AE--KEKAKKE--FQLLVVNADHG-GYTSELNAGKMS-PRLSIGDQSEPGLTG-D-  509
        |   |     |   ||  ||  ||||||    | ||| | ||   | |  |      ||  |
GHR-I   DEAQKKGKGSEDSEEKTQKELQFQLLVVDPE-GSGYTTESNARQISTPP-ST-PM--PG-SGYQT  557

GHR-II  LSPLP-PASPYHESDTTAVSP--LP--P-----APV--YTVVEGVDRQNSLLLTP--NSTPAPQL  560
          | |    |         ||  ||  |     |||  ||||  || | |||| |     | | |
GHR-I   IHPQPVETKPAATAENNQ-SPYILPDSPQSQFFAPVADYTVVQEVDSQHSLLLNPPPRQSPPPCL  621

GHR-II  II-P-KTMPT-PGGYLTPDLLGSITP  583
           | |     | || ||||||   |
GHR-I   PHHPTKALAAMPVGYVTPDLLGNLSP  647
```
