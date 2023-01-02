package esaMatcher

import(
  "github.com/evolbioinf/esa"
  "github.com/dadidange/par_lp/src/seqUtil"
)

type MyEsa struct {
  s []byte
  sa []int
  lcp []int
  cld []int
  strandSize int
}

func (e *MyEsa) Sa() []int {return e.sa}
func (e *MyEsa) Lcp() []int {return e.lcp}
func (e *MyEsa) Cld() []int {return e.cld}
func (e *MyEsa) Sequence() []byte {return e.s}
func (e *MyEsa) StrandSize() int {return e.strandSize}

type Interval struct{
  start int
  end int
  mid int
  l int
}

func NewInterval(start, end int, e MyEsa) Interval{
  //Check for empty, invalid or singleton interval
  if (start >= end){
    //singleton
    if (start >=0){
      return Interval{start, end, start, e.lcp[end]}
    } else {
      //empty or invalid
      return EmptyInterval()
    }
  }

  m := e.cld[end] //CLD.L(m+1) = cld(m)
  for (m <= start){
    m = e.cld[m]
  }

  return Interval{start, end, m, e.lcp[m]}
}

func EmptyInterval() Interval{
  return Interval{-1, -1 , -1, -1}
}

func BuildEsa(s []byte, revComp bool) MyEsa{
  strandSize := len(s)
  if revComp {
    revStr := append([]byte{'#'}, seqUtil.RevCompDna(s)...)
    s = append(s, revStr...)
  }
  s=append(s, byte('$'))
  sa := esa.Sa(s)
  lcp := esa.Lcp(s, sa)
  // Add last element to lcp if necessary
  if lcp[len(lcp)-1] != -1{
    lcp = append(lcp, -1)
  }
  cld := BuildCld(lcp)
  return MyEsa{s, sa, lcp, cld, strandSize}
}

func BuildCld(lcp []int) []int{
  stack := []int{}
  top := func() int{
    return stack[len(stack)-1]
  }
  pop := func() int{
    t:= top()
    stack = stack[:len(stack)-1]
    return t
  }
  push := func(i int) {
    stack = append(stack, i)
  }

  n := len(lcp) - 1
  cld := make([]int, n+1)
  cld[0]=n
  push(0)
  var last int

  for k := 1; k <= n; k++ {
    for lcp[k] < lcp[top()]{
      last = pop()
      for lcp[top()]==lcp[last]{
        cld[top()] = last // CLD[k].R = CLD[k]
        last = pop()
      }
      if lcp[k] < lcp[top()] {
        cld[top()] = last // CLD[k].R = CLD[k]
      } else {
        cld[k - 1] = last // CLD[k].L = CLD[k-1].R = CLD[k-1]
      }
    }
    push(k)
  }
  return cld
}


func (e *MyEsa)GetInterval(i Interval, c byte) (Interval){
  // Check Singleton Interval
  if i.start == i.end{
    if(e.s[e.sa[i.start]] == c){
      return i
    } else {
      //Return empty interval
      return EmptyInterval()
    }
  }
  lower := i.start
  upper := i.mid
  l := i.l
  for e.lcp[upper] == l {
    if (e.s[e.sa[lower]+l] == c){
      //match found
      return NewInterval(lower, upper-1, *e)
    }
    //increment interval boundaries
    lower = upper
    //check for singleton
    if (lower == i.end){
      break
    }
    upper = e.cld[upper] //CLD.R(m) = cld(m)
  }
  if (e.s[e.sa[lower] + l] == c){
    return NewInterval(lower, i.end, *e)
  } else {
    return EmptyInterval()
  }
}

func (e *MyEsa) GetMatch(query []byte) Interval{
  in := NewInterval(0, len(e.s)-1, *e)
  cld := EmptyInterval()
  k := 0
  m := len(query)
  for k < m{
    cld = e.GetInterval(in, query[k])
    if (cld.start == -1 && cld.end == -1){
      if (k == 0){
        return cld
      }
      in.l = k
      return in
    }

    k++ //the k-th character was matched in
    in = cld
    l := in.l
    if(in.start == in.end || l > m){
      l = m
    }

    for saIdx:=e.sa[in.start]; k < l; k++ {
      if(e.s[saIdx+k] != query[k]){
        in.l = k
        return in
      }
    }
  }
  in.l = m
  return in
}


func (e *MyEsa) FindAnchors(query []byte,
                      threshold int) []seqUtil.Homology{
  var homs []seqUtil.Homology
  qLen := len(query)
  strandBorder := e.StrandSize()

  lastEndQue:= 0
  currentQ := 0

  lastEndRef := 0
  lastWasRight := false
  var matchL int
  var lastLen int
  var currentRef int

  var currentMatch Interval
  currentHom := seqUtil.NewHomology(0,0,0, true)

  advance := func() {
    if matchL > 0 {
      currentQ += matchL
      return
    }
    currentQ ++
  }

  naive := func() bool{
    //Check if we can just continue comparing our queries after some mismatch
    if currentQ - lastEndQue - lastLen > threshold{
      return false
    }

    q := currentQ
    step := q - lastEndQue
    r := lastEndRef + step
    match := 0
    for r < len(e.s) && q < qLen {
      if query[q] == e.s[r] {
          match++
          q++
          r++
      } else {
          break
      }
    }
    if match > threshold{
      matchL = match
      currentRef = r - match
      return true
    }
    return false
  }

  isAnchor := func () bool{
    //match naïve
    if naive(){
      return true
    }

    //Find Regular Match
    currentMatch = e.GetMatch(query[currentQ:])
    matchL = currentMatch.l
    //uniqueness, minimum length
    if(currentMatch.start != currentMatch.end ||
              matchL < threshold){ return false }
    currentRef = e.sa[currentMatch.start]
    return true
  }


  for currentQ < qLen{
    // Find anchor
    if isAnchor(){
          isFwd := currentRef < strandBorder
          // check for anchor pair
          if(currentQ - lastEndQue == currentRef - lastEndRef &&
            isFwd == currentHom.IsFwd()){
            //set new End to start of current query + matchlen
            currentHom.Expand(currentQ+matchL)
            lastWasRight = true
          } else {
            // We can not build anchor pair
            // Start new left one and append the old
            if (lastWasRight || lastLen >= 2*threshold){
              if !currentHom.IsFwd(){
                currentHom.ReverseRefCoords(strandBorder)
              }
              homs = append(homs, currentHom)
            }
            currentHom = seqUtil.NewHomology(currentQ, currentRef, matchL, isFwd)
            lastWasRight = false
          }
          // set variables for next iteration
          lastLen = matchL
          lastEndQue = currentQ + matchL
          lastEndRef = currentRef + matchL
    }
    //proceed in query
    advance()
  }

  // Check last
  if (lastWasRight || lastLen >= 2*threshold){
    if !currentHom.IsFwd(){
      currentHom.ReverseRefCoords(strandBorder)
    }
    homs = append(homs, currentHom)
  }
  return homs
}

