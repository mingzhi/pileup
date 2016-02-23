package main

import (
	"fmt"
	"log"
	"math"
	"os"
	"runtime"

	"github.com/mingzhi/ncbiftp/genomes/profiling"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"github.com/mingzhi/pileup"
	"github.com/mingzhi/pileup/calc"
)

type cmdCt struct {
	pileupFile, fastaFile, gffFile, outFile string
	pileupFormat                            string
	codonTableID                            string
	maxl, pos, minCoverage                  int
	regionStart, regionEnd, chunckSize      int
	debug                                   bool
}

// Run is the main function.
func (cmd *cmdCt) Run() {
	// The input of pileup can be from standard input,
	// or from a file.
	var f *os.File
	if cmd.pileupFile == "" {
		f = os.Stdin
	} else {
		f = openFile(cmd.pileupFile)
	}
	defer f.Close()

	// Prepare genome position profile.
	genome := readGenome(cmd.fastaFile)
	gffs := readGff(cmd.gffFile)
	codonTable := taxonomy.GeneticCodes()[cmd.codonTableID]
	profile := profiling.ProfileGenome(genome, gffs, codonTable)
	if cmd.regionEnd <= 0 {
		cmd.regionEnd = len(profile)
	}

	// Convert pos from int to byte.
	posType := convertPosType(cmd.pos)

	// Read SNP from pileup input.
	snpChan := readPileup(f, cmd.regionStart, cmd.regionEnd, cmd.pileupFormat)

	// Apply filters.
	filteredSNPChan := cmd.filterSNP(snpChan, profile, posType)

	// Split SNPs into different chuncks.
	snpChanChan := cmd.splitChuncks(filteredSNPChan)

	// For each chunck, do the calculation.
	covsChan := cmd.doCalculation(snpChanChan)

	// Collect results from each chunck.
	csMeanVars, crMeanVars, ctMeanVars := cmd.collect(covsChan, cmd.maxl)

	// And finally, write results into the output file.
	cmd.write(csMeanVars, crMeanVars, ctMeanVars, cmd.outFile)
}

// filterSNP returns SNPs in specific positions.
func (cmd *cmdCt) filterSNP(snpChan <-chan *pileup.SNP, profile []profiling.Pos, posType byte) chan *pileup.SNP {
	c := make(chan *pileup.SNP)
	go func() {
		defer close(c)
		for s := range snpChan {
			if s.Pos >= cmd.regionStart && s.Pos < cmd.regionEnd {
				p := profile[s.Pos].Type
				if checkPosType(posType, p) {
					c <- filterOverlap(s)
				}
			}
		}
	}()
	return c
}

func filterOverlap(s *pileup.SNP) *pileup.SNP {
	m := make(map[string]byte)
	for i := range s.Alleles {
		qname := s.Alleles[i].QName
		base := s.Alleles[i].Base
		b, found := m[qname]
		if found {
			if b != base {
				m[qname] = '*'
			}
		} else {
			m[qname] = base
		}
	}

	alleles := make([]pileup.Allele, len(s.Alleles))
	k := 0
	for i := range s.Alleles {
		qname := s.Alleles[i].QName
		if m[qname] != '*' {
			alleles[k] = s.Alleles[i]
			k++
			m[qname] = '*'
		}
	}
	s.Alleles = alleles
	return s
}

// panic
func (cmd *cmdCt) panic(msg string) {
	if cmd.debug {
		log.Panic(msg)
	} else {
		log.Fatalln(msg)
	}
}

// calcInChunck does calculation in a chunck of SNPs.
func (cmd *cmdCt) calcInChunck(snpChan chan *pileup.SNP) *calc.Calculator {
	// Create job channel.
	// Each job is a array of SNPs, which we will calculate
	// correlations of the first SNP with the rest ones.
	jobChan := make(chan []*pileup.SNP)
	go func() {
		defer close(jobChan)
		arr := []*pileup.SNP{}
		for s := range snpChan {
			arr = append(arr, s)
			lag := s.Pos - arr[0].Pos

			if lag < 0 {
				cmd.panic("SNPs are not in order.")
			}

			if lag >= cmd.maxl {
				jobChan <- arr
				arr = arr[1:]
			}
		}
		jobChan <- arr
	}()

	// Make ncpu workers.
	// For each worker, do the calculation,
	// and push the result into a channel.
	ncpu := runtime.GOMAXPROCS(0)
	c := make(chan *calc.Calculator)
	for i := 0; i < ncpu; i++ {
		go func() {
			covs := calc.New(cmd.maxl)
			maxN := 10000
			xArr := make([]float64, maxN)
			yArr := make([]float64, maxN)
			for arr := range jobChan {
				cmd.calcSNPArr(arr, covs, xArr, yArr)
			}
			c <- covs
		}()
	}

	// Wait for all the worker,
	// and collect their results.
	var covs *calc.Calculator
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

// splitChuncks split the genome of SNPs into several chuncks.
// Each chunck is a channel of SNP,
// and we returns a channel of channel.
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

func (cmd *cmdCt) doCalculation(snpChanChan chan chan *pileup.SNP) chan *calc.Calculator {
	cc := make(chan *calc.Calculator)

	go func() {
		defer close(cc)
		for snpChan := range snpChanChan {
			covs := cmd.calcInChunck(snpChan)
			cc <- covs
		}
	}()

	return cc
}

// calcSNPArr calculate correlation of the first SNP with the rest.
// It finds pairs of bases that are from the same read,
// compare every two pairs of bases,
// and compute several correlations, which is contained in a calculator.
// A calculator here is a black box.
// calcSNPArr only push inputs into the calculator.
func (cmd *cmdCt) calcSNPArr(snpArr []*pileup.SNP, calculator *calc.Calculator, xArr, yArr []float64) {
	s1 := snpArr[0]
	m := make(map[string]pileup.Allele)
	for _, a := range s1.Alleles {
		if isATGC(a.Base) {
			m[a.QName] = a
		}
	}

	pairs := make([]AllelePair, len(m))
	for k := 0; k < len(snpArr); k++ {
		s2 := snpArr[k]

		// double check the lag.
		// it expects an order array of SNPs.
		l := s2.Pos - s1.Pos
		if l >= cmd.maxl {
			break
		} else if l < 0 {
			cmd.panic("SNPs is not in order.")
		}

		numPair := cmd.findPairs(m, s2.Alleles, pairs)

		// check the coverage.
		// if less than min coverage, skip.
		if numPair < cmd.minCoverage {
			continue
		}

		k := 0
		for i := 0; i < numPair && k < len(xArr); i++ {
			p1 := pairs[i]
			for j := i + 1; j < numPair && k < len(xArr); j++ {
				p2 := pairs[j]
				x := cmd.diffBases(p1.A.Base, p2.A.Base)
				y := cmd.diffBases(p1.B.Base, p2.B.Base)
				xArr[k] = x
				yArr[k] = y
				k++
			}
		}
		calculator.Increment(xArr[:k], yArr[:k], l)
	}
}

func (cmd *cmdCt) diffBases(a, b byte) float64 {
	if a != b {
		return 1.0
	} else {
		return 0.0
	}
}

type AllelePair struct {
	A, B pileup.Allele
}

func (cmd *cmdCt) findPairs(m map[string]pileup.Allele, mates []pileup.Allele, pairs []AllelePair) (numPair int) {
	k := 0
	for _, b := range mates {
		if isATGC(b.Base) {
			a, found := m[b.QName]
			if found && isATGC(b.Base) {
				pairs[k] = AllelePair{A: a, B: b}
				k++
			}
		}
	}

	numPair = k

	return
}

// collect
func (cmd *cmdCt) collect(calculatorChan chan *calc.Calculator, maxl int) (csMeanVars, crMeanVars, ctMeanVars *calc.MeanVariances) {

	csMeanVars = calc.NewMeanVariances(cmd.maxl)
	ctMeanVars = calc.NewMeanVariances(cmd.maxl)
	crMeanVars = calc.NewMeanVariances(cmd.maxl)

	meanvars := []*calc.MeanVariances{csMeanVars, crMeanVars, ctMeanVars}
	for calculator := range calculatorChan {
		for i := 0; i < calculator.MaxL; i++ {
			cs := calculator.Cs.GetMean(i)
			cr := calculator.Cr.GetResult(i)
			ct := calculator.Ct.GetResult(i)
			n := calculator.Cs.GetN(i)
			if n > 10 {
				vs := []float64{cs, cr, ct}
				for j := range vs {
					if !math.IsNaN(vs[j]) {
						meanvars[j].Increment(i, vs[j])
					}
				}
			}

		}
	}

	return
}

// write
func (cmd *cmdCt) write(csMeanVars, crMeanVars, ctMeanVars *calc.MeanVariances, filename string) {
	w, err := os.Create(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer w.Close()
	mvs := []*calc.MeanVariances{csMeanVars, crMeanVars, ctMeanVars}
	for i := 0; i < mvs[0].Size(); i++ {
		w.WriteString(fmt.Sprintf("%d\t", i))
		for _, mv := range mvs {
			w.WriteString(fmt.Sprintf("%g\t%g\t%d\t", mv.GetMean(i), mv.GetVar(i), mv.GetN(i)))
		}
		w.WriteString("\n")
	}
}
