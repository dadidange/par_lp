module github.com/dadidange/par/src/esaMatcher

go 1.17

replace github.com/dadidange/par/src/seqUtil => ../seqUtil

require (
	github.com/dadidange/par/src/seqUtil v0.0.0-00010101000000-000000000000
	github.com/evolbioinf/esa v0.0.0-20220119170951-53a316f07a9b
)

require github.com/evolbioinf/fasta v0.0.0-20220329100526-fa625fa59b5d // indirect
