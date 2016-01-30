package main

import (
	"fmt"
	"github.com/mingzhi/gomath/stat/correlation"
	"github.com/mingzhi/gomath/stat/desc/meanvar"
	"github.com/mingzhi/ncbiftp/genomes/profiling"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"math"
)

type cmdCr struct {
	piFile, fastaFile, gffFile, outFile string
	codonTableID                        string
	maxl, pos, minCoverage              int
	regionStart, regionEnd, chunckSize  int
}

func (cmd *cmdCr) Run() {
	// Obtain codon table for identifying four-fold degenerate sites.
	codonTable := taxonomy.GeneticCodes()[cmd.codonTableID]
	// Profiling genome using reference sequence and protein feature data.
	genome := readGenome(cmd.fastaFile)
	gffs := readGff(cmd.gffFile)
	profile := profiling.ProfileGenome(genome, gffs, codonTable)
	posType := convertPosType(cmd.pos)
	// Read pi.
	piChan := readPi(cmd.piFile)
	piChunckChan := cmd.split(piChan)
	covsChan := cmd.calc(piChunckChan, profile, posType, cmd.maxl)
	covMVs, xMVs, yMVs := cmd.collect(covsChan)

	cmd.write(covMVs, xMVs, yMVs)
}

func (cmd *cmdCr) collect(covsChan chan []Covariance) (covMVs, xMVs, yMVs []*meanvar.MeanVar) {
	covMVs = make([]*meanvar.MeanVar, cmd.maxl)
	xMVs = make([]*meanvar.MeanVar, cmd.maxl)
	yMVs = make([]*meanvar.MeanVar, cmd.maxl)
	for i := range covMVs {
		covMVs[i] = meanvar.New()
		xMVs[i] = meanvar.New()
		yMVs[i] = meanvar.New()
	}
	for covs := range covsChan {
		for i := range covs {
			n := covs[i].GetN()
			v := covs[i].GetResult()
			if n > 10 && !math.IsNaN(v) {
				covMVs[i].Increment(v)
				xMVs[i].Increment(covs[i].MeanX())
				yMVs[i].Increment(covs[i].MeanY())
			}
		}
	}
	return
}

func (cmd *cmdCr) split(piChan chan Pi) chan []Pi {
	c := make(chan []Pi)
	go func() {
		defer close(c)
		chunckEnd := cmd.chunckSize
		pis := []Pi{}
		for pi := range piChan {
			if pi.Pos >= chunckEnd {
				c <- pis
				pis = []Pi{}
				chunckEnd += cmd.chunckSize
			}
			pis = append(pis, pi)
		}
		c <- pis
	}()
	return c
}

func (cmd *cmdCr) calc(piChunckChan chan []Pi, profile []profiling.Pos, posType byte, maxl int) chan []Covariance {
	c := make(chan []Covariance)
	go func() {
		defer close(c)
		for chunck := range piChunckChan {
			covs := cmd.calcCr(chunck, profile, posType, maxl)
			c <- covs
		}
	}()

	return c
}

type Covariance interface {
	GetN() int
	GetResult() float64
	MeanX() float64
	MeanY() float64
	Increment(x, y float64)
}

// Calculate covariance of rates.
func (cmd *cmdCr) calcCr(pis []Pi, profile []profiling.Pos, posType byte, maxl int) (covs []Covariance) {
	corrs := make([]Covariance, maxl)
	for i := 0; i < maxl; i++ {
		corrs[i] = correlation.NewBivariateCovariance(false)
	}

	for i := 0; i < len(pis); i++ {
		p1 := pis[i]
		pos1 := profile[p1.Pos]
		if checkPosType(posType, pos1.Type) {
			for j := i; j < len(pis); j++ {
				p2 := pis[j]
				pos2 := profile[p2.Pos]

				distance := p2.Pos - p1.Pos
				if distance < 0 {
					distance = -distance
				}
				if distance >= maxl {
					break
				}

				if checkPosType(posType, pos2.Type) {
					x, y := p1.Pi(), p2.Pi()
					corrs[distance].Increment(x, y)
				}
			}
		}

	}

	covs = corrs
	return
}

func (cmd *cmdCr) write(covMVs, xMVs, yMVs []*meanvar.MeanVar) {
	w := createFile(cmd.outFile)
	defer w.Close()

	for i := 0; i < len(covMVs); i++ {
		c := covMVs[i]
		x := xMVs[i]
		y := yMVs[i]
		w.WriteString(fmt.Sprintf("%d\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n", i, c.Mean.GetResult(), c.Var.GetResult(), x.Mean.GetResult(), x.Var.GetResult(), y.Mean.GetResult(), y.Var.GetResult(), c.Mean.GetN()))
	}
}
