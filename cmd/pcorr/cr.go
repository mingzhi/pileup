package main

import (
	"fmt"
	"github.com/mingzhi/gomath/stat/correlation"
	"github.com/mingzhi/gomath/stat/desc/meanvar"
	"github.com/mingzhi/ncbiftp/genomes/profiling"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"math"
	"path/filepath"
)

type cmdCr struct {
	prefix, genomeDir                  string
	codonTableID                       string
	maxl, pos, minCoverage             int
	regionStart, regionEnd, chunckSize int
}

func (cmd *cmdCr) Run() {
	// Read pi.
	piFile := cmd.prefix + ".pi"
	piChan := readPi(piFile)
	gPiCC := cmd.separate(piChan)
	for gPiChan := range gPiCC {
		cmd.runOne(gPiChan)
	}
}

func (cmd *cmdCr) runOne(gPiChan genomePiChan) {
	// Obtain codon table for identifying four-fold degenerate sites.
	codonTable := taxonomy.GeneticCodes()[cmd.codonTableID]
	// Profiling genome using reference sequence and protein feature data.
	ref := gPiChan.genome
	fnaFile := filepath.Join(cmd.genomeDir, ref+".fna")
	gffFile := filepath.Join(cmd.genomeDir, ref+".gff")
	gffs := readGff(gffFile)
	genome := readGenome(fnaFile)
	profile := profiling.ProfileGenome(genome, gffs, codonTable)
	posType := convertPosType(cmd.pos)
	piChunckChan := cmd.split(gPiChan.piChan)
	covsChan := cmd.calc(piChunckChan, profile, posType, cmd.maxl)
	covMVs, xMVs, yMVs := cmd.collect(covsChan)

	cmd.write(ref, covMVs, xMVs, yMVs)
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

type genomePiChan struct {
	genome string
	piChan chan Pi
}

func (cmd *cmdCr) separate(piChan chan Pi) chan genomePiChan {
	cc := make(chan genomePiChan)
	go func() {
		defer close(cc)
		c := genomePiChan{}
		for pi := range piChan {
			if c.genome == "" {
				c.genome = pi.Ref
				c.piChan = make(chan Pi)
			}

			if c.genome != pi.Ref {
				cc <- c
				c = genomePiChan{}
				c.genome = pi.Ref
				c.piChan = make(chan Pi)
			}

			c.piChan <- pi
		}
		cc <- c
	}()
	return cc
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

func (cmd *cmdCr) write(ref string, covMVs, xMVs, yMVs []*meanvar.MeanVar) {
	outFile := fmt.Sprintf("%s_%s_calc_cr_%d.txt", cmd.prefix, ref, cmd.pos)
	w := createFile(outFile)
	defer w.Close()

	for i := 0; i < len(covMVs); i++ {
		c := covMVs[i]
		x := xMVs[i]
		y := yMVs[i]
		w.WriteString(fmt.Sprintf("%d\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n", i, c.Mean.GetResult(), c.Var.GetResult(), x.Mean.GetResult(), x.Var.GetResult(), y.Mean.GetResult(), y.Var.GetResult(), c.Mean.GetN()))
	}
}
