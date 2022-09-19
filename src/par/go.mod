module github.com/dadidange/par/src/par

replace github.com/dadidange/par/src/esaMatcher => ../esaMatcher

replace github.com/dadidange/par/src/seqUtil => ../seqUtil

go 1.18

require (
	github.com/dadidange/par/src/esaMatcher v0.0.0-00010101000000-000000000000
	github.com/dadidange/par/src/seqUtil v0.0.0-00010101000000-000000000000
	github.com/evolbioinf/fasta v0.0.0-20220329100526-fa625fa59b5d
)

require github.com/evolbioinf/esa v0.0.0-20220119170951-53a316f07a9b // indirect
