package main

import (
	"fmt"
	"github.com/mingzhi/gomath/stat/correlation"
	"github.com/mingzhi/gomath/stat/desc/meanvar"
	"github.com/mingzhi/ncbiftp/genomes/profiling"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"github.com/mingzhi/pileup"
	"log"
	"math"
	"os"
)

type cmdCt struct {
	pileupFile, fastaFile, gffFile, outFile string
	codonTableID                            string
	maxl, pos, minCoverage                  int
	regionStart, regionEnd, chunckSize      int
}

func (cmd *cmdCt) Run() {
	var f *os.File
	if cmd.pileupFile == "" {
		f = os.Stdin
	} else {
		f = openFile(cmd.pileupFile)
	}

	genome := readGenome(cmd.fastaFile)
	gffs := readGff(cmd.gffFile)
	codonTable := taxonomy.GeneticCodes()[cmd.codonTableID]
	profile := profiling.ProfileGenome(genome, gffs, codonTable)
	if cmd.regionEnd <= 0 {
		cmd.regionEnd = len(profile)
	}

	posType := convertPosType(cmd.pos)
	snpChan := cmd.filterSNP(readPileup(f), profile, posType)

	covsChan := cmd.calcCt(snpChan)
	meanVars := cmd.collect(covsChan, cmd.maxl)
	cmd.write(meanVars, cmd.outFile)
}

func (cmd *cmdCt) filterSNP(snpChan chan pileup.SNP, profile []profiling.Pos, posType byte) chan pileup.SNP {
	c := make(chan pileup.SNP)
	go func() {
		defer close(c)
		for s := range snpChan {
			if s.Pos >= cmd.regionStart && s.Pos < cmd.regionEnd {
				p := profile[s.Pos].Type
				if checkPosType(posType, p) {
					c <- s
				}
			}
		}
	}()
	return c
}

func (cmd *cmdCt) calcCt(snpChan chan pileup.SNP) chan []*correlation.BivariateCovariance {
	c := make(chan []*correlation.BivariateCovariance)
	go func() {
		defer close(c)
		covs := []*correlation.BivariateCovariance{}
		for i := 0; i < cmd.maxl; i++ {
			covs = append(covs, correlation.NewBivariateCovariance(false))
		}

		currentChunkEnd := cmd.chunckSize + cmd.regionStart
		snpArr := []pileup.SNP{}
		for s := range snpChan {
			if s.Pos > currentChunkEnd {
				c <- covs
				covs = []*correlation.BivariateCovariance{}
				for i := 0; i < cmd.maxl; i++ {
					covs = append(covs, correlation.NewBivariateCovariance(false))
				}
				currentChunkEnd += cmd.chunckSize
			}

			snpArr = append(snpArr, s)
			lag := s.Pos - snpArr[0].Pos
			if lag >= cmd.maxl {
				cmd.calc(snpArr, covs)
				snpArr = snpArr[1:]
			}
		}
		cmd.calc(snpArr, covs)
		c <- covs
	}()

	return c
}

func (cmd *cmdCt) calc(snpArr []pileup.SNP, covs []*correlation.BivariateCovariance) {
	s1 := snpArr[0]
	m := make(map[int]pileup.Allele)
	for _, a := range s1.Alleles {
		if isATGC(a.Base) {
			m[a.ReadID] = a
		}
	}

	for k := 0; k < len(snpArr); k++ {
		s2 := snpArr[k]
		l := s2.Pos - s1.Pos
		if l >= cmd.maxl {
			break
		}

		pairs := cmd.findPairs(m, s2.Alleles)
		if len(pairs) > cmd.minCoverage {
			for i := 0; i < len(pairs); i++ {
				for j := i + 1; j < len(pairs); j++ {
					var x, y float64
					p1 := pairs[i]
					p2 := pairs[j]
					if p1.A.Base != p2.A.Base {
						x = 1
					} else {
						x = 0
					}

					if p1.B.Base != p2.B.Base {
						y = 1
					} else {
						y = 0
					}

					covs[l].Increment(x, y)
				}
			}
		}

	}
}

type AllelePair struct {
	A, B pileup.Allele
}

func (cmd *cmdCt) findPairs(m map[int]pileup.Allele, mates []pileup.Allele) (pairs []AllelePair) {
	for _, b := range mates {
		if isATGC(b.Base) {
			a, found := m[b.ReadID]
			if found {
				pairs = append(pairs, AllelePair{A: a, B: b})
			}
		}
	}
	return
}

// collect
func (cmd *cmdCt) collect(covsChan chan []*correlation.BivariateCovariance, maxl int) (meanVars []*meanvar.MeanVar) {
	meanVars = []*meanvar.MeanVar{}
	for i := 0; i < cmd.maxl; i++ {
		meanVars = append(meanVars, meanvar.New())
	}

	for covs := range covsChan {
		for i := range covs {
			c := covs[i]
			v := c.GetResult()
			if !math.IsNaN(v) {
				meanVars[i].Increment(v)
			}
		}
	}

	return
}

// write
func (cmd *cmdCt) write(meanVars []*meanvar.MeanVar, filename string) {
	w, err := os.Create(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer w.Close()

	for i := 0; i < len(meanVars); i++ {
		m := meanVars[i].Mean.GetResult()
		v := meanVars[i].Var.GetResult()
		n := meanVars[i].Mean.GetN()
		w.WriteString(fmt.Sprintf("%d\t%g\t%g\t%d\n", i, m, v, n))
	}
}