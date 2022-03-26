# nimBoxshade
Shade sequence alignment with boxshade-like software written with Nim.

This is a pratice project since I am learning the Nim language. Right now it only takes an fasta file with aligned protein sequences and output an RTF file that can be opened in Word.

Please check the [boxshade](https://github.com/pinbo/boxshade) software or [pyBoxshade](https://github.com/mdbaron42/pyBoxshade) for more input and output formats.

The amino acid similarity matrix is from this paper:
https://doi.org/10.1186/1471-2105-10-394

# Compile
I recommend to use `nimBoxshadeMini.nim`, because it has smaller size. `nimBoxshade.nim` needs to import 'tables' module, so the size is bigger, but it is a good example of using the 'table' module ('dictionary' in other languages).
## Compile for local use
You need to install [Nim compiler](https://nim-lang.org/install.html) first. Then compile it with command:  
`nim c -d:release nimBoxshadeMini.nim`

Then check the usage by typing:  
`./nimBoxshadeMini -h`

## Compile for webassembly
I get the instruction on how to compile a Nim program to Webassembly [here](https://github.com/treeform/nim_emscripten_tutorial). You need to install [Emscripten](https://emscripten.org/docs/getting_started/downloads.html) first. Then use the command below to compile it:  
`nim c -d:emscripten -o:nimBoxshade.js ./nimBoxshadeMini.nim`

Then you can use them as the way in [biowasm](https://github.com/biowasm/biowasm). You can pass more parameters to `emcc` in the last line of file `config.nims`.

# Usage

``` sh
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
```

The current `parseopt` module does not support `-o output.rtf`, you have to write `-o=output.rtf` or `-o:output.rtf`. In the future, I will try other modules.