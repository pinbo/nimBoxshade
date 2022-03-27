## Mini did not use the tables module because it is big
## instead, it uses only arrays/sequences

## modules
# import std/tables
import sequtils 
import strutils
import std/parseopt
from os import commandLineParams, paramCount
import math

proc writeHelp() =
  let help = """
  This app takes an protein alignment file (fasta format, equal length, gaps are "-") and output a RTF file with beautiful shading.
  The methods are the same as boxshade.

  Usage: ./nimBoxshade input.fa -o=output.rtf -w=60 -c/--consensus -r/--ruler -n/--number -t/--threshold=0.5 -d/--dna
  the input file is an algined fasta file
  -o: output file name, default is output.rtf
  -w: output alignment width, default is 60
  -c: print consensus line in the output
  -r: print ruler line on the top
  -n: print start residual number for each line
  -t: the fraction of sequences that must agree for a consensus 
  -d: input is a DNA sequence alignment
  """
  quit(help)
proc writeVersion() =
#  quit(getAppFilename() & " version 0.5")
  quit("nimBoxshade: version 1.1")

# echo paramCount(), " parameters"
if paramCount() < 1:
  writeHelp()

var filename: string
var outfile = "out.rtf"
var outwidth = 60
var conflag = false # print consensus line?
var rulerflag = false # print ruler line?
var numflag = false # print start number for each line
var dnaflag = false # input is an dna alignment
var thrfrac = 0.5
# var p = initOptParser("--left --debug:3 -l -r:2")
echo "commandLineParams() is ", commandLineParams()
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
    of "dna", "d": dnaflag = true
    of "threshold", "t": thrfrac = parseFloat(val)
  of cmdEnd: assert(false) # cannot happen
if filename == "":
  # no filename has been given, so we show the help
  writeHelp()

## function to read a fasta file, return an ordered dictionary and keys
proc readFasta(infile: string): (seq[string], seq[string]) =
  var keys: seq[string]
  var allSeqs: seq[string]
  var seqs = ""
  let contents = readFile(infile).strip() 
  # echo contents
  let lines = contents.splitLines()
#   var seqDict = initOrderedTable[string, string]()
  var seqName = ""
  for ll in lines:
    var ll2 = ll.strip()
    if contains(ll2, ">"):
      seqName = ll2.replace(">", "")
      keys.add(seqName)
      if len(keys) > 1: # at least the 2nd sequences
        allSeqs.add(seqs)
        seqs = ""
    else:
      seqs.add(ll2.toUpperAscii())
  # in the end
  allSeqs.add(seqs)
  return (allSeqs, keys)

## read fasta
# var testFasta = {"seq1": "MRDRTHELRQGDN", "seq2": "MKDRLEQLKAKQL", "seq3":"MRDRLPDLTACR-"}.toTable  # creates a Table
let (testFasta, seqNames) = readFasta(filename)
# echo testFasta
# # test 
echo "seq names are ", seqNames
let nseq = len(testFasta)
echo "nseq is ", nseq
let seqLen = len(testFasta[0])
echo "seqLen is ", seqLen
var conSeq = " ".repeat(seqLen) # consensus sequence init with spaces
let thr = thrfrac*float(nseq)

# background and forground color
type
  intSeq = seq[ seq[int] ] # an seq of seq

var
  bgColorDict: intSeq
  fgColorDict: intSeq
# init colors with all 0
for i in 0 .. nseq - 1:
  var bgc: seq[int]
  for j in 0 .. seqLen - 1:
    bgc.add(0)
  bgColorDict.add(bgc)
  fgColorDict.add(bgc)

## similar residuals boxshade
# FYW 1 # ILVM 2 # RKH 3 # DE 4 # GA 5 # TS 6 # NQ 7
# let
#   aas = ['F', 'Y', 'W', 'I', 'L', 'V', 'M', 'R', 'K', 'H', 'D', 'E', 'G', 'A', 'T', 'S', 'N', 'Q']
#   grpNum = [1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7]

# similar residuals from paper
# Kim, Y., Sidney, J., Pinilla, C. et al. Derivation of an amino acid similarity matrix for peptide:MHC binding and its application as a Bayesian prior. BMC Bioinformatics 10, 394 (2009). https://doi.org/10.1186/1471-2105-10-394
# group: mutually similar: FYW 1, ILM 2, HKR 3, DE 4, AP 5, TS 6
# additional similar: F: I, I:FVLM, V:A, A:TV, T:A
var
  aas =    @['F',   'Y',  'W',  'I',     'L',  'M',  'R',  'K',  'H',  'D', 'E', 'A',   'P', 'T',  'S', 'V'] # v is an outlier
  grpNum = @[ 1,     1,    1,    2,       2,    2,    3,    3,    3,    4,   4,   5,     5,   6,    6,   7 ] # amino acid groups
  simAA =  @["WYI", "WF", "FY", "FVLM",  "IM", "IL", "KH", "RH", "KR", "E", "D", "VPT", "A", "AS", "T", "IA"]
  grpAA = @["FWY", "ILM", "HKR", "DE", "AP", "TS", "V"]

if dnaflag:
  aas =    @['A',  'G',  'R',  'C',   'T',  'Y' ]
  grpNum = @[ 1,    1,    1,    2,     2,    2  ] # amino acid groups
  simAA =  @["GR", "AR", "AG", "TY",  "CY", "CT"]
  grpAA =  @["AGR", "CTY"]
  
# proc freq [T] (ss: openArray[T]): (seq[T], seq[int], T, int) =
#   var keys: openArray[T]
#   var freqs: seq[int]
#   var maxCount = 0
#   var maxKey: T
#   for v in ss:
#     if v not in keys:
#       let cc = count(ss, v)
#       keys.add(v)
#       freqs.add(cc)
#       if cc > maxCount:
#         maxCount = cc
#         maxKey = v
#   return (keys, freqs, maxCount, maxKey)

proc findMax [T] (ss: openArray[T]): (T, int) =
  var maxCount = 0
  var maxKey: T
  for v in ss:
    let cc = count(ss, v)
    if cc > maxCount:
      maxCount = cc
      maxKey = v
  return (maxKey, maxCount)

## start to check similarity
for i in 0..(seqLen-1):
  # echo "i is ", i
  var aai: seq[char] #  aa i, the ith AA in all sequences
  var gri: seq[int] # gr i, the group of ith AA in all sequences
  for k, v in testFasta: # k is index (nth seq), v is value (sequences)
    aai.add(v[i])
    let idx = find(aas, v[i])
    if idx == -1:
      gri.add(0)
    else:
      gri.add(grpNum[idx])
  # check grp occurence
  let (largestAA, countAA) = findMax(aai)
  let (largestGrp, countGrp) = findMax(gri)
  var conAA = ' ' # consens letter
  # if largestGrp == 0: # non group AAs
  if countAA > int(thr): # there is a consensus letter
    if largestAA != '-': 
      conAA = largestAA
      if countAA != nseq:
        conAA = toLowerAscii(conAA)
  elif countGrp > int(thr) and largestGrp > 0: # there is a consensus letter
    let AAgrp = grpAA[largestGrp-1] # a string
    var maxCount = 0
    var maxAA = ' '
    for j in AAgrp:
      let countj = count(aai, j)
      if j in aai and countj > maxCount:
        maxCount = countj
        maxAA = j
    conAA = toLowerAscii(maxAA)
  ## set consensus and colors
  # color: 1=black, 2=gray, 3=white; 0 is the system default (mostly white, but I will define mine here from 1 to 3)
  conSeq[i] = conAA # set consensus at ith postion
  for k, v in testFasta: # k is index, v is value
    var m = v[i]
    # echo "m is ", m, " conAA is ", conAA
    if conAA == ' ': # no consensus, black on white
      fgColorDict[k][i] = 1
      bgColorDict[k][i] = 3
    elif countAA > int(thr): # if seq[i] == consensus, then give it white on black, others white on gray
      if m == toUpperAscii(conAA): # white on black if same as consensus
        fgColorDict[k][i] = 3
        bgColorDict[k][i] = 1
      elif toUpperAscii(conAA) in aas and m in simAA[find(aas, toUpperAscii(conAA))]: # white on gray if similar to consensus
        fgColorDict[k][i] = 3
        bgColorDict[k][i] = 2
      else:
        fgColorDict[k][i] = 1
        bgColorDict[k][i] = 3
    else: # no major AA, all similar AA will be white on gray
      if m == toUpperAscii(conAA) or (toUpperAscii(conAA) in aas and m in simAA[find(aas, toUpperAscii(conAA))]): # white on gray if similar to consensus
        fgColorDict[k][i] = 3
        bgColorDict[k][i] = 2
      else:
        fgColorDict[k][i] = 1
        bgColorDict[k][i] = 3

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
{\colortbl;\red0\green0\blue0;\red150\green150\blue150;\red255\green255\blue255;}
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
var maxNameLen = 0
if conflag:
  maxNameLen = len("consensus")
for i in seqNames:
  if len(i) > maxNameLen:
    maxNameLen = len(i)

var minLeftSpace = maxNameLen + 2 # at least 2 spaces between names and sequences
var minNumSpace = len($seqLen) # convert to seq len to string and get its len.
echo "max name length is ", maxNameLen

## AA number at the beginning of each line
var aaNumList: seq[int]
for k in 0 .. nseq - 1:
  aaNumList.add(1)

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
# var n1 = seqLen div outwidth
# let n2 = seqLen mod outwidth # mod
# if n2 > 0:
#   n1 += 1
while lend < int(ceil(seqLen / outwidth)) * outwidth:
  if lcount + nlBlock > lines_per_page:
    lcount = 0
    rtfContent.add("\\page\n")
  if lend >= seqLen:
    lend = seqLen - 1
  if rulerflag:
    var rulerLine = ' '.repeat(minLeftSpace + minNumSpace + 1) & ruler[lstart .. lend]
    # echo rulerLine
    rtfContent.add(r"\highlight3\cf1 " & rulerLine & "\n\\highlight3\\cf1 \\line\n")
  for k, v in testFasta:
    if numflag:
      numSpace = align($(aaNumList[k]), minNumSpace)
    # echo alignLeft(k, minLeftSpace), numSpace, " ", v[lstart .. lend]
    let ngap = count(v[lstart .. lend], '-')
    aaNumList[k] = aaNumList[k] + outwidth - ngap
    # rtf format
    bgc = 2
    fgc = 0
    rtfContent.add(r"\highlight3\cf1 " & alignLeft(seqNames[k], minLeftSpace) & numSpace & " ")
    for i in lstart .. lend:
      var newbgc = bgColorDict[k][i] # bg color 
      var newfgc = fgColorDict[k][i] 
      if newbgc == bgc and newfgc == fgc: # same as last aa
        rtfContent.add(v[i])
      else:
        bgc = newbgc
        fgc = newfgc
        rtfContent.add("\n\\highlight" & $bgc & "\\cf" & $fgc & " " & v[i])
    rtfContent.add("\n\\highlight3\\cf1 \\line\n") # add a newline at the end
  if conflag:
    var conLine = alignLeft("consensus", minLeftSpace) & ' '.repeat(minNumSpace + 1) & conSeq[lstart .. lend]
    # echo conLine
    rtfContent.add(r"\highlight3\cf1 " & conLine & "\n\\highlight3\\cf1 \\line\n")
  # add one blank line
  rtfContent.add("\n\\highlight3\\cf1 \\line\n")
  lstart += outwidth
  lend += outwidth
  lcount += nlBlock
  # echo "" # blank line

rtfContent.add("}")
writeFile(outfile, rtfContent)
