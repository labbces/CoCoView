# CoCoView
## Brief explanation 
CoCoView is a single python v.3 script used to generate sequence logos using codons, which allows for a more detailed analysis of sequence conservation. 

## Table of contents
* [Before running](https://github.com/labbces/CoCoView#before-running)
  * [Specifications of the sequences](https://github.com/labbces/CoCoView#specifications-of-the-sequences)
* [How to run](https://github.com/labbces/CoCoView#how-to-run)
* [Modulators - Details](https://github.com/labbces/CoCoView#modulators---details)
  * [--prefixFileName](https://github.com/labbces/CoCoView#--prefixfilename--p)
  * [--imageTitle](https://github.com/labbces/CoCoView#--imagetitle--i) 
  * [--alphaColor](https://github.com/labbces/CoCoView#--alphacolor--a)
  * [--customPaletteFile](https://github.com/labbces/CoCoView#--custompalettefile--c)
  * [--degreeOfUncertainty](https://github.com/labbces/CoCoView#--degreeofuncertainty--d)
  * [--matrixLogoType](https://github.com/labbces/CoCoView#--matrixlogotype--m)
  * [--datasetType](https://github.com/labbces/CoCoView#--datasettype--t)
  * [--logoFormat](https://github.com/labbces/CoCoView#--logoformat--l)
 * [Brief Examples](https://github.com/labbces/CoCoView#brief-examples)
   * [Brief Example 1](https://github.com/labbces/CoCoView#brief-example-1)
   * [Brief Example 2](https://github.com/labbces/CoCoView#brief-example-2)



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

### --imageTitle [-i]
This argument is a string that will appear as the title at the top of the sequence logo. If not provided by the user a title will be automatically generated from the input file name

### --alphaColor [-a]
Codons are collored accrding to the amino acid it encodes. The color pallets available are: "weblogo_protein (default)", "charge", "chemistry" and "hydrophobicity", demonstred bellow:
![alt text](https://github.com/labbces/CoCoView/blob/main/images/AlphaColors.png)
To add a custom pallet, set alphaColor as "custom" and add a json file in [-c (customPalleteFile)](https://github.com/labbces/CoCoView#--custompalettefile--c) option.

### --customPaletteFile [-c]
To use different color pattern, it is possible to add a customized one, selecting "custom" to --alphaColor [(-a)](https://github.com/labbces/CoCoView#--alphacolor--a) and the path to a JSON (.json) file to -c option. 

An option for a custom pallet can be found in [customPallet.json](https://github.com/labbces/CoCoView/blob/main/customPalette.json)


### --degreeOfUncertainty [-d]
As presented in [Specifications of the sequences](https://github.com/labbces/CoCoView#specifications-of-the-sequences), ambiguous nucleotides are allowed in the input sequence once they follow modern IUPAC nucleotide code nomenclature. However, these sequences can be filtered based on the ambiguous nucleotides present in the sequence using the degree of Uncertainty property. 
The --degreeOfUncertainty is a float number between 0 and 100, based on this number CoCoView removes the sequences which have a proportion of ambiguous nucleotides greater than it.
For example, a degreeOfUncertainty set to 30% (that is, set as 30) will exclude all sequences of length equal to 12 that have at least 4 ambiguous nucleotides.

### --matrixLogoType [-m]
--matrixLogoType has 'bit' and 'probability' as options. 
The original sequence logo proposed by [Schneider and Stephens (1990](https://dx.doi.org/10.1093%2Fnar%2F18.20.6097) utilizes the bit as the unit of conservation measurement, so it has major conceptual support. The state "bit" is the default option. 
However, to generate the 'bit' sequence logo, CoCoView constructs first a probability matrix that can be used to define the graphic representing the frequency of each codon at each position. It's important to punctuate the probability sequence logo is not supported by Schneider and Stephens's original proposition, it is a measurement of _frequency_, not conservation.

### --datasetType [-t]
 If duplicated sequences are present in the input dataset, setting this argument to ‘nonreduntant’ will remove duplicates from the analyses. This option is useful for small datasets. When very large datasets are used (thousands of sequences with hundreds/thousands of residues), users are advised to use third-party tools to generate non-redundant sequence sets, eg., cd-hit or UCLUST. Setting ‘nonreduntant’ may be of interest when the user wants to visualize rare variants. Default value ‘reduntant’.

### --logoFormat [-l]
The format in which the sequence logo will be saved. Options: 'png' and 'pdf'.

## Brief Examples
### Brief Example 1

! Files used and generated for and from example 1 can be found in [example 1 files](https://github.com/labbces/CoCoView/tree/main/test/example1)

The SARS-Cov-2 genome encodes two polyproteins that are cleaved post-translationally by proteases. CoCoView can be used to visualize the codon conservation of the protease cleavage sites. We generated a Codon conservation sequence logo from 209,219 24bp-long sequences collected from around the world, representing the four residues up- and downstream of the cleavage site between nsp6 and nsp7. 

The code ran was:

```
python3 CoCoView.py ./CoCoView/examples/example1/Bordas_World_6.nt.fasta -p "World6" -i "World 6 : nsp6/nsp7 nonredundant" -a "weblogo_protein" -d 1 -m "bit" ´t "nonredundant" -l "png" 
```
### Brief Example 2

! Files used and generated for and from example 2 can be found in [example 2 files](https://github.com/labbces/CoCoView/tree/main/examples/example2)

Transcription factors are proteins that bind DNA and regulate the expression of target genes. AP2 is a transcription factor involved in the regulation of growth and development, fruit ripening, defense response, and metabolism in plants. As transcription factors are usually represented with Sequence Logos, we generated a codon sequence logo from this AP2 gene.

```
python3 CoCoView.py examples/example2/AP2.Nta.veryShort.nt.aln.fa -p "Nicotiana_tabacumAP2" -i "Nicotiana_tabacum - AP2" -a "custom" -c customPalette.json -d 0 -m bit -t redundant -l png
```

To illustrate the benefits of a per-codon variation representation, we generated the same sequence logos using a nucleotide based sequence logo generator. At the 10th to 12th positions, which represent the 4th codon of that region of the CDS, one could incorrectly draw the conclusion that the triplet “GAT" is common at that position, based on the conservation of the individual nucleotides. However, when looking at the sequence logo based on condons (generated with CoCoView), it is clear that “GAT" is not common at all at this position. 
 
 [Comparison image](https://github.com/labbces/CoCoView/blob/main/images/Comparison%20image.png)


