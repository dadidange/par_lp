# PAR - Pile Anchors on Reference

This repository contains the code for my masters project *par* in its initial form, written in literate programming style. 
For a more recent version I split this in two. 
The library for indexing the input and for matching can be found at [github.com/dadidange/esaMatcher](github.com/dadidange/esaMatcher). 
Currently, I am working on some adaptions at the code to compute the alignments.
At some point, this will be published at github.com/dadidange/par. 

Right now, this repository here should work best. 

## How to run
### noweb
To run it from scratch type `make`. 
This will overwrite most of the files. 

For this to work, some dependencies need to be installed:
 - [noweb](https://www.cs.tufts.edu/~nr/noweb/) 
    `$ sudo apt install noweb`
 - latex 
 - [go](https://go.dev/doc/install)
 - [libdivsufsort](https://github.com/y-256/libdivsufsort)
    `$ sudo apt install libdivsufsort-dev`

The working directory should now contain the executable *par* and the corresponding code in src.
The directory [./texfiles](./texfiles) contains the document `pardoc.pdf` for documentation.
The documentation is rather a monologue that needs to be reworked. 

### just the code
Running from scratch is not necessary, [code](./src) and [documentation](./texfiles/pardoc.pdf) are already provided.  
To just compile the code, type `go build` in the `src/par` directory. 
*libdivsufsort* and *go* are still required for this. 

## What to run
*Par* expects FASTA files as input. For example, it can be applied to the 4 covid genomes in [./samples/](./samples/).
```
./par samples/*.fa > covid.maf
```
Results in a multiple alignment in [MAF-format](https://genome.ucsc.edu/FAQ/FAQformat.html#format5).

[mview](https://desmid.github.io/mview/) is a great tool to display alignments a little better than plain text. 
Mview is available here: https://github.com/desmid/mview/ and could be executed like this:

```
mview -in maf -html head -css on -coloring mismatch covid.maf > covid.html
```
