package main

import (
	"fmt"
	"log"
	"math"
	"os"

	"runtime"

	"github.com/mingzhi/gomath/stat/correlation"
	"github.com/mingzhi/gomath/stat/desc/meanvar"
	"github.com/mingzhi/ncbiftp/genomes/profiling"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"github.com/mingzhi/pileup"
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
	defer f.Close()

	genome := readGenome(cmd.fastaFile)
	gffs := readGff(cmd.gffFile)
	codonTable := taxonomy.GeneticCodes()[cmd.codonTableID]
	profile := profiling.ProfileGenome(genome, gffs, codonTable)
	if cmd.regionEnd <= 0 {
		cmd.regionEnd = len(profile)
	}

	posType := convertPosType(cmd.pos)
	snpChan := readPileup(f, cmd.regionStart, cmd.regionEnd)
	filteredSNPChan := cmd.filterSNP(snpChan, profile, posType)
	snpChanChan := cmd.splitChuncks(filteredSNPChan)
	covsChan := cmd.calcCt(snpChanChan)
	meanVars, xMVs, yMVs := cmd.collect(covsChan, cmd.maxl)
	cmd.write(meanVars, xMVs, yMVs, cmd.outFile)
}

func (cmd *cmdCt) filterSNP(snpChan <-chan *pileup.SNP, profile []profiling.Pos, posType byte) chan *pileup.SNP {
	c := make(chan *pileup.SNP)
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

func (cmd *cmdCt) calcOne(snpChan chan *pileup.SNP) []*correlation.BivariateCovariance {
	ncpu := runtime.GOMAXPROCS(0)

	// create jobs.
	jobChan := make(chan []*pileup.SNP)
	go func() {
		defer close(jobChan)
		arr := []*pileup.SNP{}
		for s := range snpChan {
			arr = append(arr, s)
			lag := s.Pos - arr[0].Pos
			if lag >= cmd.maxl {
				jobChan <- arr
				arr = arr[1:]
			}
		}
		jobChan <- arr
	}()

	c := make(chan []*correlation.BivariateCovariance)
	for i := 0; i < ncpu; i++ {
		go func() {
			covs := createCovariances(cmd.maxl)
			for arr := range jobChan {
				cmd.calc(arr, covs)
			}
			c <- covs
		}()
	}

	var covs []*correlation.BivariateCovariance
	for i := 0; i < ncpu; i++ {
		cc := <-c
		if i == 0 {
			covs = cc
		} else {
			for k := 0; k < len(covs); k++ {
				covs[k].Append(cc[k])
			}
		}
	}

	return covs
}

func (cmd *cmdCt) splitChuncks(snpChan chan *pileup.SNP) chan chan *pileup.SNP {
	cc := make(chan chan *pileup.SNP)
	go func() {
		defer close(cc)
		currentChunkEnd := cmd.chunckSize + cmd.regionStart
		c := make(chan *pileup.SNP)
		cc <- c
		for s := range snpChan {
			if s.Pos > currentChunkEnd {
				close(c)
				c = make(chan *pileup.SNP)
				cc <- c
				currentChunkEnd += cmd.chunckSize
			}
			c <- s
		}
		close(c)
	}()
	return cc
}

func (cmd *cmdCt) calcCt(snpChanChan chan chan *pileup.SNP) chan []*correlation.BivariateCovariance {
	cc := make(chan []*correlation.BivariateCovariance)

	go func() {
		defer close(cc)
		for snpChan := range snpChanChan {
			covs := cmd.calcOne(snpChan)
			cc <- covs
		}
	}()

	return cc
}

func createCovariances(maxl int) []*correlation.BivariateCovariance {
	covs := []*correlation.BivariateCovariance{}
	for i := 0; i < maxl; i++ {
		covs = append(covs, correlation.NewBivariateCovariance(false))
	}
	return covs
}

func (cmd *cmdCt) calc(snpArr []*pileup.SNP, covs []*correlation.BivariateCovariance) {
	s1 := snpArr[0]
	m := make(map[string]pileup.Allele)
	for _, a := range s1.Alleles {
		if isATGC(a.Base) {
			m[a.QName] = a
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

func (cmd *cmdCt) findPairs(m map[string]pileup.Allele, mates []pileup.Allele) (pairs []AllelePair) {
	for _, b := range mates {
		if isATGC(b.Base) {
			a, found := m[b.QName]
			if found && isATGC(b.Base) {
				pairs = append(pairs, AllelePair{A: a, B: b})
			}
		}
	}
	return
}

// collect
func (cmd *cmdCt) collect(covsChan chan []*correlation.BivariateCovariance, maxl int) (meanVars, xMVs, yMVs []*meanvar.MeanVar) {
	for i := 0; i < cmd.maxl; i++ {
		meanVars = append(meanVars, meanvar.New())
		xMVs = append(xMVs, meanvar.New())
		yMVs = append(yMVs, meanvar.New())
	}

	for covs := range covsChan {
		for i := range covs {
			c := covs[i]
			v := c.GetResult()
			x := c.MeanX()
			y := c.MeanY()
			if !math.IsNaN(v) {
				meanVars[i].Increment(v)
				xMVs[i].Increment(x)
				yMVs[i].Increment(y)
			}
		}
	}

	return
}

// write
func (cmd *cmdCt) write(meanVars, xMVs, yMVs []*meanvar.MeanVar, filename string) {
	w, err := os.Create(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer w.Close()
	mvs := [][]*meanvar.MeanVar{meanVars, xMVs, yMVs}
	for i := 0; i < len(meanVars); i++ {
		w.WriteString(fmt.Sprintf("%d\t", i))
		for _, mv := range mvs {
			w.WriteString(fmt.Sprintf("%g\t%g\t", mv[i].Mean.GetResult(), mv[i].Var.GetResult()))
		}
		w.WriteString(fmt.Sprintf("%d\n", mvs[0][i].Mean.GetN()))
	}
}
