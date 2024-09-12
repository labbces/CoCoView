# -*- coding: utf-8 -*-

# TODO: excluir seqs com certa qtde de codigos ambiguos (ns)
# TODO: fazer os gifs com ordem pre-definida

# Imports
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import logomaker
from logomaker import Glyph
from Bio import AlignIO
from os.path import exists
import re
import json
from pprint import pprint

# Using ArgParse to make easier use this script using command-line
# ver como colocar numeros, true e false
parser = argparse.ArgumentParser()
parser.add_argument(
    "fastaFile", help="path to fasta file, must have extension fasta, or fa")
parser.add_argument('-p', '--prefixFileName',
                    help="String to use as prefix of output file names, if not provided will use fastaFile as prefix")
parser.add_argument('-i', '--imageTitle',
                    help="String to use as title in the image")
parser.add_argument('-a', '--alphaColor', choices=['weblogo_protein', 'charge', 'chemistry',
                    'hydrophobicity', 'custom'], help='Alphabet color scheme', default='weblogo_protein')
parser.add_argument('-d', '--degreeOfUncertainty', type=float, default=0.0,
                    help='Proportion of ambiguous nucleotides allowed in the sequences to use in the logo')
parser.add_argument('-m', '--matrixLogoType', help='Type of matrix to built',
                    default='bit', choices=['bit', 'probability'])
parser.add_argument('-t', '--datasetType', help='whether to built a reduntand o non-redundant dataset before creating the logo',
                    default='redundant', choices=['redundant', 'nonredundant'])
parser.add_argument('-l', '--logoFormat', help='Format of the sequence logo',
                    default='png', choices=['png', 'pdf'])
parser.add_argument('-c', '--customPaletteFile',
                    help='JSON file with a custom palette to use in the logo')


args = parser.parse_args()

# Checking that input file exists, and has the right extension
if exists(args.fastaFile):
    if args.fastaFile.endswith('.fa') or args.fastaFile.endswith('.fasta'):
        fastaFile = args.fastaFile
    else:
        print(f'ERROR: {args.fastaFile} must have extension .fasta or .fa')
        parser.print_help()
        exit()
else:
    print(f'ERROR: {args.fastaFile} does not exist')
    parser.print_help()
    exit()

if args.alphaColor == 'custom' and args.customPaletteFile:
    if exists(args.customPaletteFile):
        True
    else:
        print(f'ERROR: When choosing --alphaColor custom you must provide an existing JSON file with the custom palette')
        parser.print_help()
        exit()

elif args.alphaColor in ['weblogo_protein', 'charge', 'chemistry', 'hydrophobicity']:
    True

else:
    print(f'ERROR: You must use both \'--alphaColor custom\' and \'--customPaletteFile file.json\'')


# Check if the user wants and can use a custom palette file
if args.customPaletteFile:  # Check if user called the argument
    if args.alphaColor != 'custom':
        print(f'ERROR: If proding a JSON file with a color pallete, you must use the option --alphaColor set to \'custom\'')
        parser.print_help()
        exit()
    else:
        if exists(args.customPaletteFile):  # Check whether the file exists
            print(args.customPaletteFile)
            try:
                with open(args.customPaletteFile, 'r') as f:
                    customPalette = json.load(f)
                    print(customPalette)
            except ValueError as e:
                print(f'ERROR: Invalid JSON in fille {args.customPaletteFile}')
                parser.print_help()
                exit()

# Did we get a prefix to create outfiles from the user?
# If not use the basename of the input file as prefix
if args.prefixFileName:
    prefixFileName = args.prefixFileName
else:
    prefixFileName = re.sub('.fa(sta)?$', '', fastaFile)

# Did we get an image title  from the user?
# If not use the basename of the input file with the string CoCoViewer
if args.imageTitle:
    imageTitle = args.imageTitle + ' ' + args.datasetType + ' CoCoViewer'
else:
    imageTitle = re.sub('.fa(sta)?', '', fastaFile) + \
        ' ' + args.datasetType + ' CoCoViewer'
    print(f'{imageTitle}')
# GET SEQUENCES
seqDict = {}
alignment = AlignIO.read(fastaFile, "fasta")
SeqLength = alignment.get_alignment_length()
AlphaColor = str(args.alphaColor).lower()

# Check that the alignment has complete codons, length should be multiple of 2

if SeqLength % 2 != 0:
    print(f'Your sequence length {args.SeqLength} is not a multiple of 2')
    exit()


for record in alignment:
    notKnownNucleotide = 0
    for nucleotide in record:
        if nucleotide.upper() not in "ATGC":
            notKnownNucleotide += 1
    degreeOfUncertainty = 100*notKnownNucleotide/SeqLength
    if degreeOfUncertainty <= args.degreeOfUncertainty:
        if len(record.seq) == int(SeqLength):
            if args.datasetType.upper().strip() == "REDUNDANT":
                if str(record.seq) in seqDict.keys():
                    seqDict[str(record.seq)] = seqDict[str(record.seq)] + 1
                else:
                    seqDict[str(record.seq)] = 1
            elif args.datasetType.upper().strip() == "NONREDUNDANT":
                if str(record.seq) not in seqDict.keys():
                    seqDict[str(record.seq)] = 1
            else:
                print(
                    '! Please check if "DataSetType" is spelled correctly - e.g "Redundant or NonRedundant"')
                exit()
        else:
            print('The sequence size is not as expected', fastaFile,
                  '\t', record.id, " e.g This sequence won't be used")
    else:
        print(f'The sequence degree of uncertainty {degreeOfUncertainty} is higher than the expected {args.degreeOfUncertainty}', fastaFile,
              ' ', record.id, " e.g This sequence won't be used")

# print(f'Dicionário das sequências:\n{seqDict}\n')

# Printing how many sequeces were used to construct the Graphs
totalSeqs = sum(seqDict.values())

if totalSeqs == 0:
    print('All sequences were filtered out, nothing remained to generate the Codon Sequence Logo')
    exit()
else:
    print(
        f'\nA total of {totalSeqs} sequence(s) was/were used to construct the codon sequence logo')

# IMPORTANT DICTS
matrixDict = {'AA': [], 'AT': [], 'AC': [], 'AG': [], 'TA': [], 'TT': [], 'TC': [], 'TG': [], 'CA': [], 'CT': [],
              'CC': [], 'CG': [], 'GA': [], 'GT': [], 'GC': [], 'GG': []
              }
Codon2Symbol = {'AA': "A", 'AT': "B", 'AC': "C", 'AG': "D", 'TA': "E", 'TT': "F", 'TC': "G", 'TG': "H", 'CA': "I", 'CT': "J",
                'CC': "K", 'CG': "L", 'GA': "M", 'GT': "N", 'GC': "O", 'GG': "P"
                }
extraSymbols = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l",
                "m", "n", "o", "p", "q", "r", "s", "u", "v", "w", "x", "y", "z"]
color_palett = {
    'weblogo_protein':
    {'GCT': 'black', 'GCC': 'black', 'GCA': 'black', 'GCG': 'black', 'TTT': 'black', 'TTC': 'black', 'ATT': 'black', 'ATC': 'black', 'ATA': 'black', 'TTA': 'black', 'TTG': 'black', 'CTT': 'black', 'CTC': 'black', 'CTA': 'black', 'CTG': 'black', 'ATG': 'black', 'CCT': 'black', 'CCC': 'black', 'CCA': 'black', 'CCG': 'black', 'TGG': 'black', 'GTT': 'black', 'GTC': 'black', 'GTA': 'black', 'GTG': 'black', 'TGT': 'lime', 'TGC': 'lime', 'GGT': 'lime', 'GGC': 'lime', 'GGA': 'lime', 'GGG': 'lime', 'TCT': 'lime',
        'TCC': 'lime', 'TCA': 'lime', 'TCG': 'lime', 'AGT': 'lime', 'AGC': 'lime', 'ACT': 'lime', 'ACC': 'lime', 'ACA': 'lime', 'ACG': 'lime', 'TAT': 'lime', 'TAC': 'lime', 'GAT': 'red', 'GAC': 'red', 'GAA': 'red', 'GAG': 'red', 'CAT': 'blue', 'CAC': 'blue', 'AAA': 'blue', 'AAG': 'blue', 'CGT': 'blue', 'CGC': 'blue', 'CGA': 'blue', 'CGG': 'blue', 'AGA': 'blue', 'AGG': 'blue', 'AAT': 'magenta', 'AAC': 'magenta', 'CAA': 'magenta', 'CAG': 'magenta', 'TAA': 'peru', 'TAG': 'gold', 'TGA': 'azure',

        'AA': "black", 'AT': "black", 'AC': "black", 'AG': "black", 'TA': "black", 'TT': "black", 'TC': "black", 'TG': "black", 'CA': "black", 'CT': "black",
              'CC': "black", 'CG': "black", 'GA': "black", 'GT': "black", 'GC': "black", 'GG': "black"},
    'charge':
    {'GCT': 'dimgrey', 'GCC': 'dimgrey', 'GCA': 'dimgrey', 'GCG': 'dimgrey', 'TTT': 'dimgrey', 'TTC': 'dimgrey', 'ATT': 'dimgrey', 'ATC': 'dimgrey', 'ATA': 'dimgrey', 'TTA': 'dimgrey', 'TTG': 'dimgrey', 'CTT': 'dimgrey', 'CTC': 'dimgrey', 'CTA': 'dimgrey', 'CTG': 'dimgrey', 'ATG': 'dimgrey', 'CCT': 'dimgrey', 'CCC': 'dimgrey', 'CCA': 'dimgrey', 'CCG': 'dimgrey', 'TGG': 'dimgrey', 'GTT': 'dimgrey', 'GTC': 'dimgrey', 'GTA': 'dimgrey', 'GTG': 'dimgrey', 'TGT': 'dimgrey', 'TGC': 'dimgrey', 'GGT': 'dimgrey', 'GGC': 'dimgrey', 'GGA': 'dimgrey', 'GGG': 'dimgrey', 'TCT': 'dimgrey',
        'TCC': 'dimgrey', 'TCA': 'dimgrey', 'TCG': 'dimgrey', 'AGT': 'dimgrey', 'AGC': 'dimgrey', 'ACT': 'dimgrey', 'ACC': 'dimgrey', 'ACA': 'dimgrey', 'ACG': 'dimgrey', 'TAT': 'dimgrey', 'TAC': 'dimgrey', 'GAT': 'red', 'GAC': 'red', 'GAA': 'red', 'GAG': 'red', 'CAT': 'mediumblue', 'CAC': 'mediumblue', 'AAA': 'mediumblue', 'AAG': 'mediumblue', 'CGT': 'mediumblue', 'CGC': 'mediumblue', 'CGA': 'mediumblue', 'CGG': 'mediumblue', 'AGA': 'mediumblue', 'AGG': 'mediumblue', 'AAT': 'dimgrey', 'AAC': 'dimgrey', 'CAA': 'dimgrey', 'CAG': 'dimgrey', 'TAA': 'peru', 'TAG': 'gold', 'TGA': 'azure'},
    'hydrophobicity':
    {'GCT': 'green', 'GCC': 'green', 'GCA': 'green', 'GCG': 'green', 'TTT': 'black', 'TTC': 'black', 'ATT': 'black', 'ATC': 'black', 'ATA': 'black', 'TTA': 'black', 'TTG': 'black', 'CTT': 'black', 'CTC': 'black', 'CTA': 'black', 'CTG': 'black', 'ATG': 'black', 'CCT': 'green', 'CCC': 'green', 'CCA': 'green', 'CCG': 'green', 'TGG': 'black', 'GTT': 'black', 'GTC': 'black', 'GTA': 'black', 'GTG': 'black', 'TGT': 'black', 'TGC': 'black', 'GGT': 'green', 'GGC': 'green', 'GGA': 'green', 'GGG': 'green', 'TCT': 'green', 'TCC': 'green', 'TCA': 'green', 'TCG': 'green',
        'AGT': 'green', 'AGC': 'green', 'ACT': 'green', 'ACC': 'green', 'ACA': 'green', 'ACG': 'green', 'TAT': 'black', 'TAC': 'black', 'GAT': 'mediumblue', 'GAC': 'mediumblue', 'GAA': 'mediumblue', 'GAG': 'mediumblue', 'CAT': 'green', 'CAC': 'green', 'AAA': 'mediumblue', 'AAG': 'mediumblue', 'CGT': 'mediumblue', 'CGC': 'mediumblue', 'CGA': 'mediumblue', 'CGG': 'mediumblue', 'AGA': 'mediumblue', 'AGG': 'mediumblue', 'AAT': 'mediumblue', 'AAC': 'mediumblue', 'CAA': 'mediumblue', 'CAG': 'mediumblue', 'TAA': 'peru', 'TAG': 'gold', 'TGA': 'azure'},
    'chemistry':
    {'GCT': 'black', 'GCC': 'black', 'GCA': 'black', 'GCG': 'black', 'TTT': 'black', 'TTC': 'black', 'ATT': 'black', 'ATC': 'black', 'ATA': 'black', 'TTA': 'black', 'TTG': 'black', 'CTT': 'black', 'CTC': 'black', 'CTA': 'black', 'CTG': 'black', 'ATG': 'black', 'CCT': 'black', 'CCC': 'black', 'CCA': 'black', 'CCG': 'black', 'TGG': 'black', 'GTT': 'black', 'GTC': 'black', 'GTA': 'black', 'GTG': 'black', 'TGT': 'lime', 'TGC': 'lime', 'GGT': 'lime', 'GGC': 'lime', 'GGA': 'lime', 'GGG': 'lime', 'TCT': 'lime', 'TCC': 'lime', 'TCA': 'lime',
        'TCG': 'lime', 'AGT': 'lime', 'AGC': 'lime', 'ACT': 'lime', 'ACC': 'lime', 'ACA': 'lime', 'ACG': 'lime', 'TAT': 'lime', 'TAC': 'lime', 'GAT': 'red', 'GAC': 'red', 'GAA': 'red', 'GAG': 'red', 'CAT': 'mediumblue', 'CAC': 'mediumblue', 'AAA': 'mediumblue', 'AAG': 'mediumblue', 'CGT': 'mediumblue', 'CGC': 'mediumblue', 'CGA': 'mediumblue', 'CGG': 'mediumblue', 'AGA': 'mediumblue', 'AGG': 'mediumblue', 'AAT': 'purple', 'AAC': 'purple', 'CAA': 'purple', 'CAG': 'purple', 'TAA': 'peru', 'TAG': 'gold', 'TGA': 'azure'}
}

if args.customPaletteFile:
    color_palett.update(customPalette)

# DEFINING THE POSITIONS OF THE CODONS IN THE SEQUENCE
posDict = {}
codons = []
ntAmount = 2
newCodons = 0
for seq, amount in seqDict.items():
    for cropPos in range(0, len(seq), ntAmount):
        codon = (seq[cropPos: cropPos + ntAmount]).upper()

        codons.append(codon)
    for codon_index in range(int((len(seq)) / ntAmount)):
        if args.degreeOfUncertainty == 100:
            if codons[codon_index].upper() not in matrixDict.keys():
                try:
                    Codon2Symbol[codons[codon_index]] = extraSymbols[newCodons]
                    matrixDict[codons[codon_index]] = []
                    color_palett[AlphaColor][codons[codon_index]] = 'silver'
                    newCodons += 1
                except:
                    print('To many new not-known codons. Not enough extra symbols.')
                    pass

        if codon_index not in posDict.keys():
            if codons[codon_index] in matrixDict.keys():
                posDict[codon_index] = {codons[codon_index]: seqDict[seq]}
        else:
            if codons[codon_index] in matrixDict.keys():
                if codons[codon_index] in posDict[codon_index].keys():
                    posDict[codon_index][codons[codon_index]
                                         ] = posDict[codon_index][codons[codon_index]] + seqDict[seq]
                else:
                    posDict[codon_index][codons[codon_index]] = seqDict[seq]
    codons = []
# print(f'Dicionário das Posições:\n{posDict}\n')

# BUILDING THE MATRICES
for pos in posDict.keys():
    total = sum(posDict[pos].values())
    for res in sorted(matrixDict.keys()):
        if res in posDict[pos].keys():
            freq = posDict[pos][res] / total
            # print(f'{pos} {res} {freq} {total} {posDict[pos][res]}')
            matrixDict[res].append(freq)
        else:
            matrixDict[res].append(0)
# print(f'Dicionário que formará a Matrix: \n{matrixDict}\n')

matrixSimb = {}
for codon in matrixDict.keys():
    matrixSimb[Codon2Symbol[codon]] = matrixDict[codon]
# print(f'Dicionário que formará a Matrix com símbolos: \n{matrixSimb}\n')

# Build Probability Matrix
MatrixProb = pd.DataFrame(matrixDict)
MatraixProb_name = prefixFileName + '.' + \
    args.datasetType + '.probability.matrix'
MatrixProb.to_csv(sep="\t", header=True,
                  path_or_buf=MatraixProb_name, index=True)
# print(f'Matrix Prob: \n{MatrixProb}\n')

if args.matrixLogoType.upper().strip() == "BIT":
    # Build probability symbol matrix
    MatrixProbSymbol = pd.DataFrame(matrixSimb)
    # matrixSymbol_name = prefixFileName + '.' + \args.datasetType + '.probability' + '.symbol'
    # MatrixProbSymbol.to_csv(sep="\t", header=True, path_or_buf=matrixSymbol_name, index=True)
    # print(f'Matrix de Símbolos: \n{MatrixProbSymbol}\n')

    # Converting probability matrix to information (bits) matrix
    matrixValid = logomaker.validate_matrix(
        MatrixProbSymbol, matrix_type='probability', allow_nan=True)
    matrixBit = logomaker.transform_matrix(
        matrixValid, from_type='probability', to_type='information')
    matrixBit_name = prefixFileName + '.' + args.datasetType + '.bits.matrix'
    matrixBit.to_csv(sep="\t", header=True,
                     path_or_buf=matrixBit_name, index=True)
    # print(f'Matrix de Bits:\n{matrixBit}\n')


# TRANSLATING SYMBOL BIT MATRICES TO CODON BIT MATRICES
# Dictionary symbols to codons
if args.matrixLogoType.upper().strip() == "BIT":
    Symbols2Codon = {}
    for symbol in matrixBit.columns:
        Codon = list(Codon2Symbol.keys())[
            list(Codon2Symbol.values()).index(symbol)]
        Symbols2Codon[symbol] = Codon
    # print(f'Symbols 2 codon: \n {Symbols2Codon}')

    bit_matrix_codon = matrixBit.rename(Symbols2Codon, axis='columns')
    matrixCodonBit_name = prefixFileName + '.' + args.datasetType + '.bits.matrix'
    bit_matrix_codon.to_csv(sep="\t", header=True,
                            path_or_buf=matrixCodonBit_name, index=True)
    # print(f'Bit matrix codon: \n{bit_matrix_codon}\n')


# TRANSFORMING MATRICES TO DICTIONARIES TO DEFINE GHYPH PARAMETERS
if args.matrixLogoType.upper().strip() == "BIT":
    matrixInfo = matrixBit.to_dict('Dict')
    dictCodon = {}
    for simb in matrixInfo.keys():
        for k, v in Codon2Symbol.items():
            dictCodon[k] = matrixInfo[v]
    # print(f'Dict dos bits dos codons:\n{dictCodon}\n')
elif args.matrixLogoType.upper().strip() == "PROBABILITY":
    dictCodon = MatrixProb.to_dict('Dict')
    # print(f"Matrix Prob info \n {dictCodon}")

nucleotide_colors = {
    'A': 'red',     # Example color for Adenine
    'T': 'blue',    # Example color for Thymine
    'G': 'green',   # Example color for Guanine
    'C': 'yellow',  # Example color for Cytosine
    # You can add more colors for any ambiguous nucleotides if necessary
}

# DEFINNG GLYPHS PARAMETERS
info2Glyph = {}

first_pos = True
for pos in range(0, int((SeqLength)/2)):
    a = True
    codonbits = {}
    floor = ceiling = 0
    for codon in dictCodon.keys():
        codonbits[codon] = dictCodon[codon][pos]
    codonBits_Order = sorted(
        codonbits.items(), key=lambda x: x[1], reverse=False)
    for TupleCodonBit in codonBits_Order:
        di_nt = TupleCodonBit[0]
        bit = TupleCodonBit[1]
        if bit != 0.0:
            # print(f'TupleCodonBit: {TupleCodonBit}')
            floor = ceiling
            ceiling = (floor + bit)
            # print(f'Floor: {floor}')
            # print(f'Ceiling: {ceiling}')
            p = (pos*2)+1
            for nt in di_nt:
                if first_pos:
                    p_processed = p
                    first_pos = False
                else:
                    p_processed = p + 1
                    first_pos = True
                if p_processed not in info2Glyph.keys():
                    info2Glyph[p_processed] = {di_nt: {nt: {
                        'bit': bit, 'floor': floor, 'ceiling': ceiling, 'p': p_processed, 'color': nucleotide_colors[nt]}}}
                else:
                    info2Glyph[p_processed][di_nt] = {nt: {
                        'bit': bit, 'floor': floor, 'ceiling': ceiling, 'p': p_processed, 'color': nucleotide_colors[nt]}}


# pprint(info2Glyph)

# removing positions and build list to fill glyphs
ListInfo2Glyph = []
for value in info2Glyph.values():
    for value2 in value.values():
        for nt, data in value2.items():
            # print(value2)
            ListInfo2Glyph.append({
                'codon': nt,
                'data': data,
            })

# PreListInfo2Glyph = list(info2Glyph.values())
x = SeqLength/2
fig, ax = plt.subplots(figsize=[x, 4])
# set bounding box
ax.set_xlim([0.5, (SeqLength + 1)])
ax.set_ylim([0, 2])


def generate_glyph(c, p, width, ceiling, ax, floor, color):
    return Glyph(c=c,
                 p=p,
                 width=width,
                 ceiling=ceiling,
                 ax=ax,
                 floor=floor,
                 color=color)


list(
    map(
        lambda item:
        generate_glyph(
            c=item['codon'],
            p=item['data']['p'],
            width=1.0,
            ceiling=item['data']['ceiling'],
            ax=ax,
            floor=item['data']['floor'],
            color=item['data']['color']),
        ListInfo2Glyph
    )
)
title = f'{imageTitle}'
ax.set_title(title)
ax.set_xlabel('length')
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

if args.matrixLogoType.upper().strip() == "BIT":
    ax.set_ylabel('information (bits)')
elif args.matrixLogoType.upper().strip() == "PROBABILITY":
    ax.set_ylabel('probability')
    ax.set_ylim([0, 1])

figsave_name = prefixFileName + '.' + \
    args.datasetType + '.' + args.matrixLogoType + \
    '.CoCoView' + '.' + args.logoFormat
fig.savefig(figsave_name, transparent=False, format=args.logoFormat)
plt.show()
