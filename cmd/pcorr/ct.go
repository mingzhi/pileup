package main

import (
	"fmt"
	"github.com/mingzhi/ncbiftp/genomes/profiling"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"github.com/mingzhi/pileup"
	"github.com/mingzhi/pileup/cov"
	"log"
	"os"
	"runtime"
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
	csMeanVars, crMeanVars, ctMeanVars := cmd.collect(covsChan, cmd.maxl)
	cmd.write(csMeanVars, crMeanVars, ctMeanVars, cmd.outFile)
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

func (cmd *cmdCt) calcOne(snpChan chan *pileup.SNP) *cov.Calculator {
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

	c := make(chan *cov.Calculator)
	for i := 0; i < ncpu; i++ {
		go func() {
			covs := cov.NewCalculator(cmd.maxl)
			for arr := range jobChan {
				cmd.calc(arr, covs)
			}
			c <- covs
		}()
	}

	var covs *cov.Calculator
	for i := 0; i < ncpu; i++ {
		cc := <-c
		if i == 0 {
			covs = cc
		} else {
			covs.Append(cc)
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

func (cmd *cmdCt) calcCt(snpChanChan chan chan *pileup.SNP) chan *cov.Calculator {
	cc := make(chan *cov.Calculator)

	go func() {
		defer close(cc)
		for snpChan := range snpChanChan {
			covs := cmd.calcOne(snpChan)
			cc <- covs
		}
	}()

	return cc
}

func (cmd *cmdCt) calc(snpArr []*pileup.SNP, calculator *cov.Calculator) {
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
			xArr := []float64{}
			yArr := []float64{}
			for i := 0; i < len(pairs); i++ {
				for j := i + 1; j < len(pairs); j++ {
					p1 := pairs[i]
					p2 := pairs[j]
					if p1.A.Base != p2.A.Base {
						xArr = append(xArr, 1)
					} else {
						xArr = append(xArr, 0)
					}

					if p1.B.Base != p2.B.Base {
						yArr = append(yArr, 1)
					} else {
						yArr = append(yArr, 0)
					}
				}
			}
			calculator.Calc(xArr, yArr, l)
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
func (cmd *cmdCt) collect(calculatorChan chan *cov.Calculator, maxl int) (csMeanVars, crMeanVars, ctMeanVars *cov.MeanVariances) {

	csMeanVars = cov.NewMeanVariances(cmd.maxl)
	ctMeanVars = cov.NewMeanVariances(cmd.maxl)
	crMeanVars = cov.NewMeanVariances(cmd.maxl)

	for calculator := range calculatorChan {
		for i := 0; i < calculator.MaxL; i++ {
			cs := calculator.Cs.GetMean(i)
			cr := calculator.Cr.GetResult(i)
			ct := calculator.Ct.GetResult(i)
			csMeanVars.Increment(i, cs)
			crMeanVars.Increment(i, cr)
			ctMeanVars.Increment(i, ct)
		}
	}

	return
}

// write
func (cmd *cmdCt) write(csMeanVars, crMeanVars, ctMeanVars *cov.MeanVariances, filename string) {
	w, err := os.Create(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer w.Close()
	mvs := []*cov.MeanVariances{csMeanVars, crMeanVars, ctMeanVars}
	for i := 0; i < mvs[0].Size(); i++ {
		w.WriteString(fmt.Sprintf("%d\t", i))
		for _, mv := range mvs {
			w.WriteString(fmt.Sprintf("%g\t%g\t%d\t", mv.GetMean(i), mv.GetVar(i), mv.GetN(i)))
		}
		w.WriteString("\n")
	}
}
