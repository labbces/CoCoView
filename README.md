# CoCoView
## Brief explanation 
CoCoView is a single python v.3 script used to generate sequence logos using codons, which allows for a more detailed analysis of sequence conservation. 

## Before running:
CoCoView is a python script that requires **python v.3 or greater** (it has been tested on python v3.7 and v3.9). Other than that, there are two needs to run this script: **Nucleotide sequences in FASTA file format** (see its specifications below) and **the script [CoCoView.py](https://github.com/labbces/CoCoView/blob/main/CoCoView.py)**. Besides that, CoCoView relies on external libraries (**argparse, pandas, matplotlib, Logomaker, and biopython**) that should be installed in advance. All of them are available to be installed using PyPI:

```
pip install logomaker

pip install argparse
pip install pandas
pip install matplotlib
pip install biopython
```
### Specifications of the sequences.
The sequences should be in a **FASTA format file** and must contain only **equal-sized sequences** (aligned) whose **length is a multiple of three**, i.e., only complete codons are allowed. The script uses codons as three followed symbols, for this reason, the **nucleotides symbols must be as following the modern IUPAC nucleotide code nomenclature for incompletely specified bases**, that is:

| **IUPAC nucleotide code** | **Base** |
|---|---|
| A | Adenine |
| C | Cytosine |
| G | Guanine |
| T (or U) | Thymine (or Uracil)| 
| R | A or G |
| Y | C or T |
| S | G or C |
| W | A or T |
| K | G or T |
| M | A or C |
| B | C or G or T |
| D | A or T or G |
| H | A or C or T |
| V | A or C or G |
| N | A or C or T or G |
| . or - | gap | 

[Reference](https://www.bioinformatics.org/sms/iupac.html)

Lastly, not an obligation, but based on [Crooks, G. E, 2004. WebLogo: A Sequence Logo Generator](https://doi.org/10.1101/gr.849004) we recommend using at least 40 sequences to avoid underestimation of entropy.

## How to run:
Once you have installed the script, all its dependencies, and have the nucleotide sequences according to the specifications listed above, you can run CoCoView **on the command-line interface** informing the **fasta File (Mandatory) and the graph modulators (Optional)**. 

**Defaults:**
| Modulator | Abbreviation | Default |
| --- | --- | --- |
| prefixFileName | -p | Declared fastaFile name |
| alphaColor | -a | weblogo_protein |
| degreeOfUncertainity | -d | 0.0 |
| matrixLogoType | -m | bit | 
| datasetType | -t | redundant | 
| logoFormat | -l | png |

## Modulators - Details:
_! Remember: All the modulators are optional, if not informed on the command-line the default will be applied._

### --prefixFileName [-p]
  CoCoView produces two (matrix type selected as 'Probability') or three (for matrix type selected as 'Bit') output files.
  The **sequence logo** generated for both cases will be saved as **prefixFileName + datasetType + matrixLogoType + "CoCoView" + logoFormat**. Another file generated for both cases is the **probability matrix**, it will be saved as **prefixFileName + datasetType + "probability.matrix"**
  Especifically for matrix type as "bit" the **bit matrix** would be generated and saved as **prefixFileName + datasetType + "bits.matrix"**.
 
| Output files | Name Pattern |
| --- | --- |
| Sequence logo | prefixFileName + datasetType + matrixLogoType + "CoCoView" + logoFormat |
| Probability Matrix | prefixFileName + datasetType + "probability.matrix |
| Bit Matrix | prefixFileName + datasetType + "bits.matrix |

 _! (see --matrixLogoType, --datasetType, and --logoFormat)._

### imageTitle 	[-i]
This argument is a string that will appear as the title at the top of the sequence logo. If not provided by the user a title will be automatically generated from the input file name







