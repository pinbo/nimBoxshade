## modules
import std/tables
from std/sequtils import zip
import strutils
import std/parseopt
import os
import math

proc writeHelp() =
  let help = """
  This app takes an protein alignment file (fasta format, equal length, gaps are "-") and output a RTF file with beautiful shading.
  The methods are the same as boxshade.

  Usage: ./nimBoxshade input.fa -o=output.rtf -w=60 -c/--consensus -r/--ruler -n/--number -t/--threshold=0.5
  the input file is an algined fasta file
  -o: output file name, default is output.rtf
  -w: output alignment width, default is 60
  -c: print consensus line in the output
  -r: print ruler line on the top
  -n: print start residual number for each line
  -t: the fraction of sequences that must agree for a consensus 
  """
  quit(help)
proc writeVersion() =
  quit(getAppFilename() & " version 0.5")

# echo paramCount(), " parameters"
if paramCount() < 1:
    writeHelp()

var filename: string
var outfile = "out.rtf"
var outwidth = 60
var conflag = false # print consensus line?
var rulerflag = false # print ruler line?
var numflag = false # print start number for each line
var thrfrac = 0.5
# var p = initOptParser("--left --debug:3 -l -r:2")
var p = initOptParser( commandLineParams() )
for kind, key, val in p.getopt():
  case kind
  of cmdArgument:
    filename = key
  of cmdLongOption, cmdShortOption:
    case key
    of "help", "h": writeHelp()
    of "version", "v": writeVersion()
    of "o": outfile = val
    of "w": outwidth = parseInt(val)
    of "consensus", "c": conflag = true
    of "ruler", "r": rulerflag = true
    of "number", "n": numflag = true
    of "threshold", "t": thrfrac = parseFloat(val)
  of cmdEnd: assert(false) # cannot happen
if filename == "":
  # no filename has been given, so we show the help
  writeHelp()

## similar residuals boxshade
# FYW 1 # ILVM 2 # RKH 3 # DE 4 # GA 5 # TS 6 # NQ 7
# let
#   aas = ['F', 'Y', 'W', 'I', 'L', 'V', 'M', 'R', 'K', 'H', 'D', 'E', 'G', 'A', 'T', 'S', 'N', 'Q']
#   grps = [1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7]

# similar residuals from paper
# Kim, Y., Sidney, J., Pinilla, C. et al. Derivation of an amino acid similarity matrix for peptide:MHC binding and its application as a Bayesian prior. BMC Bioinformatics 10, 394 (2009). https://doi.org/10.1186/1471-2105-10-394
# group: mutually similar: FYW 1, ILM 2, HKR 3, DE 4, AP 5, TS 6, NQ 7
# additional similar: F: I, I:FVLM, V:A, A:TV, T:A
let
  aas = ['F', 'Y', 'W', 'I', 'L', 'M', 'R', 'K', 'H', 'D', 'E', 'A', 'P', 'T', 'S', 'N', 'Q']
  grps = [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7]

var grpDict = initTable[char, int]()
for pairs in zip(aas, grps):
  let (aa, grp) = pairs
  grpDict[aa] = grp

# var grpDict2 = {1:"FYW", 2:"ILVM", 3:"RKH", 4:"DE",5: "AG", 6:"ST", 7:"NQ"}.toTable # boxshade grp
let grpDict2 = {1:"FWY", 2:"ILM", 3:"HKR", 4:"DE",5: "AP", 6:"TS", 7:"NQ"}.toTable # Kim et al. grp
let simDict = {'F':"I", 'I':"FVLM", 'V':"AI", 'A':"TV", 'T':"A"}.toTable # additional similar residuals

## function to read a fasta file, return an ordered dictionary and keys
proc readFasta(infile: string): (OrderedTable[string, string], seq[string]) =
  var keys: seq[string] 
  let contents = readFile(infile).strip() 
  # echo contents
  let lines = contents.splitLines()
  var seqDict = initOrderedTable[string, string]()
  var seqName = ""
  for ll in lines:
    var ll2 = ll.strip()
    if contains(ll2, ">"):
      seqName = ll2.replace(">", "")
      keys.add(seqName)
      seqDict[seqName] = ""
    else:
      seqDict[seqName].add(ll2.toUpperAscii())
  return (seqDict, keys)

## read fasta
# var testFasta = {"seq1": "MRDRTHELRQGDN", "seq2": "MKDRLEQLKAKQL", "seq3":"MRDRLPDLTACR-"}.toTable  # creates a Table
let (testFasta, seqNames) = readFasta(filename)
# # test making consensus
# proc getKeys(myMap: OrderedTable): seq[string] =
#   var mykeys: seq[string]
#   for key in myMap.keys():
#     mykeys.add(key)
#   return mykeys

# # test 
# let seqNames = getKeys(testFasta)
echo "seq names are ", seqNames
var nseq = len(testFasta)
echo "nseq is ", nseq
let seqLen = len(testFasta[seqNames[0]])
echo "seqLen is ", seqLen
var conSeq = " ".repeat(seqLen) # consensus sequence init with spaces
let thr = thrfrac*float(nseq)

# background and forground color dictionary
# black = 0 = \red0\green0\blue0
# gray  = 1 = \red180\green180\blue180
# white = 2 = \red255\green255\blue255
# let colors = [r"\red0\green0\blue0", r"\red180\green180\blue180", r"\red255\green255\blue255"]
var bgColorDict = initTable[string, seq[int]]()
var fgColorDict = initTable[string, seq[int]]()
# ini: first value is junk
for k in testFasta.keys():
  fgColorDict[k] = @[0]
  bgColorDict[k] = @[0]

## start to check similarity
for i in 0..(seqLen-1):
  # echo "i is ", i
  var aai: seq[char] #  aa i, the ith AA in all sequences
  var gri: seq[int] # gr i, the group of ith AA in all sequences
  for k, v in testFasta:
    aai.add(v[i])
    gri.add(grpDict.getOrDefault(v[i]))
  # check grp occurence
  # echo "aai is", aai
  # echo "gri is", gri
  let aaFreq = toCountTable(aai)
  let grpFreq = toCountTable(gri)
  # echo "grpFreq ", grpFreq
  let (largestAA, countAA) = largest(aaFreq)
  let (largestGrp, countGrp) = largest(grpFreq)
  # echo largestAA
  # echo countAA
  var conAA = ' ' # consens letter
  # if largestGrp == 0: # non group AAs
  if countAA > int(thr): # there is a consensus letter
    if largestAA != '-': 
      conAA = largestAA
      if countAA != nseq:
        conAA = toLowerAscii(conAA)
  elif countGrp > int(thr) and largestGrp > 0: # there is a consensus letter
    let AAgrp = grpDict2[largestGrp] # a string
    var maxCount = 0
    var maxAA = ' '
    for j in AAgrp:
      if j in aaFreq and aaFreq[j] > maxCount:
        maxCount = aaFreq[j]
        maxAA = j
    conAA = toLowerAscii(maxAA)
  ## set consensus and colors
  # color: 0=black, 1=gray, 2=white
  conSeq[i] = conAA
  for k, v in testFasta:
    # echo "i is ", i, "and fgColorDict[k] is ", fgColorDict[k]
    var m = v[i]
    if conAA == ' ': # no consensus, black on white
      fgColorDict[k].add(0)
      bgColorDict[k].add(2)
    elif countAA > int(thr): # if seq[i] == consensus, then give it white on black, others white on gray
      if m == toUpperAscii(conAA): # white on black
        fgColorDict[k].add(2)
        bgColorDict[k].add(0)
      elif m in grpDict and grpDict[m] == largestGrp: # white on gray
        fgColorDict[k].add(2)
        bgColorDict[k].add(1)
      elif toUpperAscii(conAA) in simDict and m in simDict[toUpperAscii(conAA)]: # white on gray for additional similar residaual with the consensus
        fgColorDict[k].add(2)
        bgColorDict[k].add(1)
      else:
        fgColorDict[k].add(0)
        bgColorDict[k].add(2)
    else: # no major AA, all similar AA will be white on gray
      if m in grpDict and grpDict[m] == largestGrp:
        fgColorDict[k].add(2)
        bgColorDict[k].add(1)
      else: # black on white
        fgColorDict[k].add(0)
        bgColorDict[k].add(2)

# echo conSeq
# for k, v in testFasta:
#   echo v
#   echo join(fgColorDict[k][1..^1]) # 2 to end, because the first is just init
#   echo join(bgColorDict[k][1..^1])

# write to rtf
var rtfContent = """
{\rtf1\ansi\deff0
{\fonttbl{\f0\fmodern Courier New;}}
{\info{\author BOXSHADE}{\title output.rtf}}
{\colortbl
\red0\green0\blue0;\red180\green180\blue180;\red255\green255\blue255;}
\paperw11880\paperh16820\margl1000\margr500
\margt910\margb910\sectd\cols1\pard\plain
\fs20

"""

let
  dev_miny = 1.0
  dev_maxy = 15000.0
  fontSize = 10.0
  dev_ysize = fontSize * 20.0
  lines_per_page = int((dev_maxy - dev_miny) / dev_ysize)

# var minLeftSpace = 16 # names and line numbers on the left
var maxNameLen = len("consensus")
for i in seqNames:
  if len(i) > maxNameLen:
    maxNameLen = len(i)
# if minLeftSpace < maxNameLen:
var minLeftSpace = maxNameLen + 2 # at least 2 spaces between names and sequences
var minNumSpace = len(($seqLen)) # convert to seq len to string and get its len.
echo "max name length is ", maxNameLen

# var formattedSeqNames = seqNames
# for i in seqNames:
#   formattedSeqNames[i] = alignLeft(i, minLeftSpace) # make up minLeftSpace with spaces on the right
## AA number at the beginning of each line
var aaNumDict = initTable[string, int]()
for k in testFasta.keys():
  aaNumDict[k] = 1

# rulers
var ruler = '.'.repeat(seqLen)
for n in countup(4, seqLen-1, 10):
  ruler[n] =  ':' # every 5th is :
for n in countup(9, seqLen-1, 10):
  ruler[n+1-len($(n+1)) .. n] = $(n+1)
# echo ruler
# outwidth = 60
var
  lstart = 0 # line start number
  lend = outwidth - 1 # line end number

var numSpace = "" # in case of no need to print numbers
if not numflag:
  minNumSpace = 0
# rtf format
var bgc = 2
var fgc = 0
var nlBlock = nseq + 1 # number of lines per block = nseq + 1 blank line
if rulerflag:
  nlBlock += 1
if rulerflag:
  nlBlock += 1
if nlBlock > lines_per_page:
  nlBlock = 0 # if number of sequences > lines_per_page, then do not check
var lcount = 0 # line count for new page
while lend < int(ceil(seqLen / outwidth)) * outwidth:
  if lcount + nlBlock > lines_per_page:
    lcount = 0
    rtfContent.add("\\page\n")
  if lend >= seqLen:
    lend = seqLen - 1
  if rulerflag:
    var rulerLine = ' '.repeat(minLeftSpace + minNumSpace + 1) & ruler[lstart .. lend]
    # echo rulerLine
    rtfContent.add(r"\highlight2\cf0 " & rulerLine & "\n\\highlight2\\cf0 \\line\n")
  for k, v in testFasta:
    if numflag:
      numSpace = align($(aaNumDict[k]), minNumSpace)
    # echo alignLeft(k, minLeftSpace), numSpace, " ", v[lstart .. lend]
    let ngap = count(v[lstart .. lend], '-')
    aaNumDict[k] = aaNumDict[k] + outwidth - ngap
    # rtf format
    bgc = 2
    fgc = 0
    rtfContent.add(r"\highlight2\cf0 " & alignLeft(k, minLeftSpace) & numSpace & " ")
    for i in lstart .. lend:
      var newbgc = bgColorDict[k][i+1] # bg color 1st is junk
      var newfgc = fgColorDict[k][i+1] 
      if newbgc == bgc and newfgc == fgc: # same as last aa
        rtfContent.add(v[i])
      else:
        bgc = newbgc
        fgc = newfgc
        rtfContent.add("\n\\highlight" & $(bgc) & "\\cf" & $(fgc) & " " & v[i])
    rtfContent.add("\n\\highlight2\\cf0 \\line\n") # add a newline at the end
  if conflag:
    var conLine = alignLeft("consensus", minLeftSpace) & ' '.repeat(minNumSpace + 1) & conSeq[lstart .. lend]
    # echo conLine
    rtfContent.add(r"\highlight2\cf0 " & conLine & "\n\\highlight2\\cf0 \\line\n")
  # add one blank line
  rtfContent.add("\n\\highlight2\\cf0 \\line\n")
  lstart += outwidth
  lend += outwidth
  lcount += nlBlock
  # echo "" # blank line

rtfContent.add("}")
writeFile(outfile, rtfContent)
