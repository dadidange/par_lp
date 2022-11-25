package main
import (
  "flag"
  "github.com/evolbioinf/fasta"

  "os"
  "strings"
  "unicode"
  "sort"
  "github.com/dadidange/par/src/esaMatcher"
  "github.com/dadidange/par/src/seqUtil"
  "runtime"
  "strconv"
  "fmt"
  "math"
)

const version = "0.21.5-LP-Fin"

func main(){
  var optR = flag.String("r", "", "reference sequence.\n\t"+
                  "Only the first sequence of this file will be handled as reference")
  var optT = flag.Int("t", 0, "threshold for minimum anchor length.\n\t" +
    "If <=0 is specified it will be calculated according to the length of the sequence.")
  //add numCPU
  var optC = flag.Int("c", 0, "Number of Cores to use for alignment.\n\t" +
    "If < 1 it runs with the number of logical CPUs usable by the current process (runtime.NumCPU()).")

  var optRevComp = flag.Bool("revComp", true, "Build reference index for reverse complement.\n")

  var optV = flag.Bool("v", false, "print par version, do not run")

  //TODO add seed? -> Don't Need since we have no random elements
  flag.Parse()

  if *optV {
    fmt.Fprintf(os.Stderr, "par v.%s\n", version)
    return
  }
  seqFiles := flag.Args()
  if len(seqFiles) < 1 {
    printUsage()
    return
  }

  if *optR != ""{
    r:= *optR
    for i, v := range seqFiles {
      if v == r {
        seqFiles = append(seqFiles[:i], seqFiles[i+1:]...)
      }
    }
    seqFiles = append([]string{r}, seqFiles...)
  }
  sequences:= readSequenceFiles(seqFiles)

  //choosing the reference
  var reference fasta.Sequence
  if *optR != ""{
    reference = sequences[0]
    sequences = sequences[1:]
  }else{
    sort.Slice(sequences, func(i,j int) bool{
      return len(sequences[i].Data()) < len(sequences[j].Data())
    })
    refIdx := len(sequences)/2
    reference = sequences[refIdx]
    sequences = append(sequences[0:refIdx], sequences[refIdx+1:]...)
  }
  go Log(fmt.Sprintf("Selected %s for reference + %d queries.\n",
          reference.Header(), len(sequences)))

  threshold := *optT
  if threshold < 1 {
    threshold = AnchorThreshold(reference.Data())
  }
  go Log(fmt.Sprintf("Set threshold to %d\n", threshold))

  SetNumCpus(*optC)

  //Construct the ESA
  myesa := esaMatcher.BuildEsa(reference.Data(), *optRevComp)
  go Log(fmt.Sprintf("finished building ESA\r"))

  ch := make(chan struct{})

  n := len(sequences)
  homs := make([][]seqUtil.Homology, n)

  for i := 0; i < n; i++ {
    go func (i int){
      h := myesa.FindAnchors(sequences[i].Data(), threshold)
      sort.Slice(h, func(i,j int) bool{
       return h[i].StartsLeftOnRef(h[j])
      })
      h =  seqUtil.FilterOverlaps(h)
      homs[i] = h
      Log(fmt.Sprintf("Fin %d\r", i))
      //Report Completion
      ch <- struct{}{}
    }(i)
  }

  // Wait for completion
  for range sequences {
    <- ch
  }
  go Log(fmt.Sprintf("Finished matching\r"))
  blocks, foundAny := pileBlocks(&homs, n)
  if !foundAny {
      Log("no match found, stopping with empty alignment.\nTry a different reference sequence.\n")
      return
  }
  //MAF Header
  fmt.Printf("%s\n",
    "##maf version=1 scoring=none")
  fmt.Println(ToMafInfo(reference.Header(), threshold, *optRevComp))
  numBlocks := len(blocks)
  for i, b := range blocks{
    go Log(fmt.Sprintf("Building Block %d of %d\r", i+1, numBlocks))
    fmt.Println(b.MafString(reference, &sequences))
  }
}

func readSequenceFiles(seqFiles []string) []fasta.Sequence {
sequences := make([]fasta.Sequence, 0)
for _, file := range seqFiles {
  f := openFile(file)
  //Scan File
  sc := fasta.NewScanner(f)
  //Append to existing sequences
  for sc.ScanSequence(){
    seq := sc.Sequence()
    f := func(c rune) bool {
      return unicode.IsSpace(c) || c == '|' || c == ','
    }
    h := strings.FieldsFunc(seq.Header(), f)[0]
    data := strings.ToUpper(string(seq.Data()))
    fMap := func(c rune) rune{
      if !(c == 'A' || c == 'C' || c == 'G' || c == 'T'){
      return 'N'
      }
      return c
    }
    data = strings.Map(fMap, data)
    seq = fasta.NewSequence(h, []byte(data))
    sequences = append(sequences, *seq)
  }
  f.Close()
}
return sequences

}

func openFile(file string) *os.File {
  f, err := os.Open(file)
  if err != nil {
    Log(fmt.Sprintf("couldn’t open %q\nError:%s", file, err))
  }
  return f
}

func Log(msg string){
  fmt.Fprintf(os.Stderr, msg)
}

func printUsage(){
  fmt.Fprintf(os.Stderr, "par v.%s\nUsage: par [OPTIONS] FILES\n\nOptions:\n", version)
  flag.PrintDefaults()
}
func AnchorThreshold(refSeq []byte) int{
  l := len(refSeq)
  gc := seqUtil.NumGC(refSeq)
  p := float64(gc) / (2*float64(l))
  pp := 1 - 0.025
  x:= 0
  for seqUtil.ShustrProb(x, l, p) < pp {
    x ++
  }
  return x
}

func SetNumCpus(n int){
  max := runtime.NumCPU()
  if n > max {
    runtime.GOMAXPROCS(max)
    s := fmt.Sprintf("Could not set numCPU to %d since max is %d.", n, max)
    Log(fmt.Sprintf("Warn: %s Set to max.\n", s))
    return
  }
  if n < 1 {
    n = max
  }
  runtime.GOMAXPROCS(n)
  Log(fmt.Sprintf("Set numCPU to %d\n", n))
}

func pileBlocks(homs *[][]seqUtil.Homology, numSeqs int) ([]MafBlock, bool) {

  finElements := 0
  nextElement := make([]int, numSeqs)

  //Check not to extend homologies
  checkNext := func (next, i int) bool{
    if next >= len((*homs)[i]){
      //reached end of this sequence
      finElements ++
      return false
    }
    return true
  }

  for i, next := range nextElement{
    checkNext(next, i)
  }
  if finElements >= numSeqs {
    return nil, false
  }

  var blocks []MafBlock
  var minStart, minStartIdx int

  SetNextMin := func(){
    minStart = math.MaxUint32
    for i, next := range nextElement{
      if next < len((*homs)[i]){
        if((*homs)[i][next].StartR() < minStart){
          minStart = (*homs)[i][next].StartR()
          minStartIdx = i
        }
      }
    }
  }

  for finElements < numSeqs {
    SetNextMin()

    h := (*homs)[minStartIdx][nextElement[minStartIdx]]
    nextElement[minStartIdx]++

    b := NewMafBlock(minStart, minStartIdx, h)
    added := true
    for added {
      added = false
      finElements = 0
      for i, next := range nextElement{
        if !checkNext(next, i){
          continue
        }
        if b.Contains((*homs)[i][next]){
          b.AddItem(i, (*homs)[i][next])
          nextElement[i]++
          added = true
        }
      }
    }
    blocks = append(blocks, b)
  }

return blocks, true
}

type MafBlock struct{
  items map[int][]seqUtil.Homology
  ends []int
  start int
  maxEnd int
}

func NewMafBlock(start, idx int, h seqUtil.Homology) MafBlock {
  b := MafBlock{make(map[int][]seqUtil.Homology), []int{}, start, 0}
  b.AddItem(idx, h)
  return b
}

func (b *MafBlock) AddItem(idx int, h seqUtil.Homology){
  b.items[idx] = append(b.items[idx], h)
  e:= h.EndR()
  b.ends = append(b.ends, e)
  if e > b.maxEnd{
    b.maxEnd = e
  }
}

func (b *MafBlock) Contains(h seqUtil.Homology) bool{
  return h.StartR() >= b.start && h.StartR() < b.maxEnd
}

func (b *MafBlock) String() string{
  return fmt.Sprintf("MafBlock: Start=%d, End=%d, ends=%v",
                                        b.start, b.maxEnd, b.ends)
}

func (b *MafBlock) MafString(refSeq fasta.Sequence, queries *[]fasta.Sequence) string{

  sort.Ints(b.ends)
  b.ends = seqUtil.RemoveDuplicateInt(b.ends)

  var blockstr strings.Builder
  items := &b.items
  start := b.start
  refName := refSeq.Header()
  refData := refSeq.Data()

  //Iterate through all end points
  for _,end := range b.ends{
    //take next homology from all items
    s := fmt.Sprintf("a score=0\ns %-10s\t %d %d + %d %s\n",
          refName, start, end-start, len(refData), refData[start:end])
    blockstr.WriteString(s)
    for idx, homs := range *items{
      h := homs[0]
        if h.StartR() < end{
          s := GetMafLine(h, start, end, (*queries)[idx])
          blockstr.WriteString(s)
          //remove from items if ends
          if h.EndR() <= end {
            if len(homs) > 1{
              (*items)[idx] = (*items)[idx][1:]
            } else {
              delete(*items, idx)
            }
          }

        }
      }
      blockstr.WriteString(fmt.Sprintln())
      start = end
    }
    return blockstr.String()
}

func GetMafLine(h seqUtil.Homology, start, end int, seq fasta.Sequence) string{
  //Init fields for MAF sequence line
  src := seq.Header()
  size := end - start
  startQ := h.StartQ()
  endQ := h.EndQ()

  var seqStr []byte
  var text string
  var strand string

  srcSize := len(seq.Data())

  overlap := h.StartR() - start
  gap := ""

  handleFwd := func ()  {
    if overlap < 0{
      //skip positions -> h begins before current interval
      startQ = startQ - overlap
      endQ = startQ + size
    }else{
      //insert gap and reduce size of actual aligning region
      gap = strings.Repeat("-", overlap)
      size = size - overlap
      endQ = startQ + size
    }
    seqStr = seq.Data()[startQ: endQ]
    strand = "+"
  }

  handleRev := func ()  {
      if overlap < 0{
        //skip positions at the end this time since they get reverted
        end= h.EndQ() + overlap
        startQ = end - size
      }else{
        gap = strings.Repeat("-", overlap)
        size = size - overlap
        end = h.EndQ()
        startQ = end-size
      }
      s := seq.Data()[startQ: end]
      seqStr = seqUtil.RevCompDna(s)
      strand = "-"
      startQ = srcSize - end
    }

  if h.IsFwd(){
    handleFwd()
  }else{
    handleRev()
  }

  text = gap + string(seqStr)

  return fmt.Sprintf("s %-10s\t %d %d %s %d %s\n",
              src, startQ, size, strand, srcSize, text)
}

func ToMafInfo(ref string, thres int, includeReverse bool) string{
  s := fmt.Sprintf("#par v.%s, command: par -r %s -t %s -revComp=%t\n",
                  version, ref, strconv.Itoa(thres), includeReverse)
  return s
}
