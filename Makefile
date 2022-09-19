
all: tangle weave

weave:
	noweave -n -x par.nw > texfiles/par.tex
	make -C texfiles
tangle:
	notangle -ResaMatcher.go par.nw > ./src/esaMatcher/esaMatcher.go
	notangle -Rpar.go par.nw > ./src/par/par.go
	notangle -RseqUtil.go par.nw > ./src/seqUtil/seqUtil.go
	make -C src/seqUtil
	make -C src/esaMatcher
	make -C src/par


clean:
	make -C texfiles clean
	make -C src/seqUtil clean 
	make -C src/esaMatcher clean
	make -C src/par clean
	rm ./par
