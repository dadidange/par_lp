package seqUtil

import (
  "fmt"
  "math"
  "math/big"
  "bytes"
)

type Homology struct{
  startQ int
  startR int
  homLen int
  isFwd bool
}

func (h *Homology) Len() int{ return h.homLen}
func (h *Homology) StartR() int{ return h.startR}
func (h *Homology) EndR() int{ return h.startR + h.homLen}
func (h *Homology) StartQ() int{ return h.startQ}
func (h *Homology) EndQ() int{ return h.startQ + h.homLen}
func (h *Homology) IsFwd() bool{ return h.isFwd}

func NewHomology(startQ, startR, homLen int, isFwd bool) Homology{
  return Homology{startQ, startR, homLen, isFwd}
}

func (h *Homology) Expand(newQEnd int){
  h.homLen = newQEnd - h.startQ
}

func (this *Homology) String() string{
  dir := "fwd"
  if !this.isFwd{
    dir = "rev"
  }
  s := fmt.Sprintf("Homology: startR=%d, startQ=%d, len=%d, dir=%s",
        this.startR, this.startQ, this.homLen, dir)
  return s
}

func (this *Homology) StartsLeftOnRef(other Homology) bool{
  return this.startR < other.startR
}

func (this *Homology) StartsLeftOnQ(other Homology) bool{
  return this.startQ < other.startQ
}

func (this *Homology) EndsLeftOnRef(other Homology) bool{
  return this.EndR() <= other.startR
}

func (this *Homology) Overlaps(other Homology) bool{
  if this.startR == other.startR {return true}

  if this.StartsLeftOnRef(other) {
    return !this.EndsLeftOnRef(other)
  } else {
    return other.EndsLeftOnRef(*this)
  }
}

func (h *Homology) ReverseRefCoords(refLen int){
  startR := h.startR
  h.startR = refLen - (startR - refLen) - h.homLen + 1
}

func FilterOverlaps(homs []Homology) []Homology{
  m := len(homs)
  if m == 1 {
    // nothing to filter
    return homs
  }
  // Initialize slices for scores and predecessors
  pred := make([]int, m)
  score := make([]int, m)
  predNum := make([]int, m)
  totalMax := -1
  totalMaxIdx := 0

  pred[0] = -1
  score[0] = homs[0].Len()

  for i := 1; i<m; i++ {
      maxIdx := -1 // or NaN
      maxScore := 0
      for j := 0; j<i; j++{
      if(homs[j].EndsLeftOnRef(homs[i]) &&
                          score[j] > maxScore){
        maxIdx = j
        maxScore = score[j]
      }
    }
    pred[i] = maxIdx
    sc := maxScore + homs[i].Len()
    score[i] = sc
    if (maxIdx >= 0){
      predNum[i] = predNum[maxIdx] + 1
    }
    if (sc > totalMax){
      totalMax = sc
      totalMaxIdx = i
    }
  }

  n := predNum[totalMaxIdx]
  reHoms := make([]Homology, n+1)
  h := totalMaxIdx
  for i:=n; i >=0; i--{
    reHoms[i] = homs[h]
    h = pred[h]
  }
  return reHoms
}

func RemoveDuplicateInt(intSlice []int) []int {
    allKeys := make(map[int]bool)
    sl := []int{}
    for _, item := range intSlice {
        if _, value := allKeys[item]; !value {
            allKeys[item] = true
            sl = append(sl, item)
        }
    }
    return sl
}

//Given the length l of a sequence and x ShustrProb returns probabilityt hat a
// match between two random, unrelated sequences of lenght l is below x with 2p
// as the GC-content of the sequence.
func ShustrProb(x, l int, p float64) float64{
        xFl := float64(x)
        lFl := float64(l)

        lExp2 := math.Pow(2,xFl)
        bin := func (n,k int) float64 {
                v := Binomial(int64(n), int64(k))
                return float64(v)
        }

        var sum float64 = 0
        for k:= 0; k<=x; k++{
                kFl := float64(k)
                tmp := math.Pow(p,kFl) * math.Pow(0.5 - p, xFl - kFl)
                b := bin(x,k)
                sum += lExp2 * b * (tmp * math.Pow(1 - tmp, lFl))


                if sum >= 1.0 {
                        return 1.0
                }
        }
        return sum
}

func Binomial(n, k int64) int64{
        z := new(big.Int).Binomial(n,k)
        return z.Int64()
}
//Counts the number of nucleotides that are either 'G' or 'C'
func NumGC(seq []byte) int{
        c := byte('C')
        g := byte('G')
        count := 0
        for _,char := range seq {
                if char == c || char == g {
                        count ++
                }
        }
        return count
}


func RevCompDna(seq []byte) []byte{
  n := len(seq)
        revSeq := make([]byte, n)

        //Reverse
        //for i, j := 0, n-1; i < j; i, j = i+1, j-1 {
    //revSeq[i], revSeq[j] = seq[j], seq[i]
        //}
  //slower but more secure
  for i:= 0; i<n; i++ {
    revSeq[(n-i)-1] = seq[i]
  }

        //Complement
        f := func (r rune) rune  {
          switch{
          case r == 'A':
                return 'T'
          case r == 'T':
                return 'A'
          case r == 'G':
                return 'C'
          case r == 'C':
                return 'G'
          default:
                return 'N'
          }
        }
        return  bytes.Map(f, revSeq)
}
