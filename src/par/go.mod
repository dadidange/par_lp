module github.com/dadidange/par_lp/src/par

replace github.com/dadidange/par_lp/src/esaMatcher => ../esaMatcher

replace github.com/dadidange/par_lp/src/seqUtil => ../seqUtil

go 1.18

require (
	github.com/dadidange/par_lp/src/esaMatcher v0.0.0-00010101000000-000000000000
	github.com/dadidange/par_lp/src/seqUtil v0.0.0-20221130163224-c1340b459c0e
	github.com/evolbioinf/fasta v0.0.0-20220329100526-fa625fa59b5d
)

require github.com/evolbioinf/esa v0.0.0-20220119170951-53a316f07a9b // indirect
