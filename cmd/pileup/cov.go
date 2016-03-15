package main

import (
	"bytes"
	"github.com/bmatsuo/lmdb-go/lmdb"
	"github.com/mingzhi/biogo/pileup"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"gopkg.in/vmihailenco/msgpack.v2"
	"log"
	"math"
)

type cmdCr struct {
	dbfile        string
	codonID       string
	featureDbPath string
	minDepth      int
	gc            *taxonomy.GeneticCode

	env        *lmdb.Env
	sizeDB     int64
	featureEnv *lmdb.Env
}

func (c *cmdCr) run() {
	numDB := 10
	c.sizeDB = 1 << 30
	var err error
	c.env, err = createEnv(c.dbfile, numDB, c.sizeDB)
	for lmdb.IsMapFull(err) {
		c.sizeDB *= 2
		c.env, err = createEnv(c.dbfile, numDB, c.sizeDB)
	}
	raiseError(err)
	defer c.env.Close()

	// open feature db
	var sizeDB int64 = 1 << 30

	c.featureEnv, err = createNoLockEnv(c.featureDbPath, numDB, sizeDB)
	for lmdb.IsMapFull(err) {
		sizeDB *= 2
		c.featureEnv, err = createNoLockEnv(c.featureDbPath, numDB, sizeDB)
	}
	raiseError(err)
	defer c.featureEnv.Close()

	createDBI(c.env, "cr")
	snpChan := c.readSNPs()
	crChan := c.calculateCr(snpChan)
	c.load(crChan)
}

func (c *cmdCr) readSNPs() chan SNPArr {
	return getAllSNPs(c.env)
}

func (c *cmdCr) calculateCr(in chan SNPArr) chan CovRes {
	out := make(chan CovRes)
	fn := func(tx *lmdb.Txn) error {
		for snpArr := range in {
			dbi, err := tx.OpenDBI("feature", 0)
			if err != nil {
				return err
			}
			// get feature
			k := snpArr.Key
			feature := Feature{}
			v1, err := tx.Get(dbi, k)
			if err != nil {
				return err
			}
			if err := msgpack.Unmarshal(v1, &feature); err != nil {
				return err
			}

			seqLen := feature.End - feature.Start + 1

			piArr := make([]float64, seqLen)
			for i := range piArr {
				piArr[i] = math.NaN()
			}

			for _, s := range snpArr.Arr {
				if len(s.Bases) < c.minDepth {
					continue
				}
				pos := s.Position - feature.Start
				if feature.IsComplementaryStrand() {
					pos = seqLen - 1 - pos
				}

				if (pos+1)%3 == 0 {
					pi := s.Pi()
					piArr[pos] = pi
				}
			}

			cc := xcross(piArr)
			out <- CovRes{Key: k, Values: cc}
		}

		return nil
	}
	go func() {
		defer close(out)
		c.featureEnv.View(fn)
	}()
	return out
}

type CovRes struct {
	Key    []byte
	Values []float64
}

func (c *cmdCr) load(crChan chan CovRes) {
	bufferSize := 1000
	buffer := []CovRes{}
	for cr := range crChan {
		if len(buffer) > bufferSize {
			err := c.loadBuffer(buffer)
			if err != nil {
				log.Panicln(err)
			}
			buffer = []CovRes{}
		}
		buffer = append(buffer, cr)
	}
	err := c.loadBuffer(buffer)
	if err != nil {
		log.Panicln(err)
	}
}

func (c *cmdCr) loadBuffer(buffer []CovRes) error {
	fn := func(txn *lmdb.Txn) error {
		dbi, err := txn.OpenDBI("cr", 0)
		if err != nil {
			return err
		}
		for _, cr := range buffer {
			key := cr.Key
			value, err := msgpack.Marshal(cr.Values)
			if err != nil {
				return err
			}
			if err := txn.Put(dbi, key, value, 0); err != nil {
				return err
			}
		}
		return nil
	}

	err := c.env.Update(fn)
	return err

}

func xcross(a []float64) []float64 {
	res := []float64{}
	for l := 0; l < len(a); l++ {
		xs, ys := []float64{}, []float64{}
		for i := 0; i < len(a)-l; i++ {
			if !(math.IsNaN(a[i]) || math.IsNaN(a[i+l])) {
				xs = append(xs, a[i])
				ys = append(ys, a[i+l])
			}
		}
		res = append(res, cov(xs, ys))
	}

	return res
}

func cov(x, y []float64) float64 {
	xy := 0.0
	xbar := 0.0
	ybar := 0.0
	for i := range x {
		xy += x[i] * y[i]
		xbar += x[i]
		ybar += y[i]
	}
	n := float64(len(x))
	if n < 100 {
		return math.NaN()
	}

	return xy/n - xbar/n*ybar/n
}

func isSynomous(s pileup.SNP, pos int, sequence []byte, gc *taxonomy.GeneticCode) bool {
	codon := sequence[(pos/3)*3 : (pos/3+1)*3]
	p := pos % 3
	c := codon[p]
	a := gc.Table[string(codon)]

	if *debug {
		log.Printf("%c, %s\n", c, s.Bases)
	}

	for _, b := range bytes.ToUpper(s.Bases) {
		if b != c {
			codon[p] = b
			a1 := gc.Table[string(codon)]
			if a != a1 {
				return false
			}
		}
	}

	return true
}

// reverse and complement a nucleotide sequence
func revcomp(seq []byte) []byte {
	s := make([]byte, len(seq))
	// reverse
	for i, j := 0, len(seq)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = seq[j], seq[i]
	}

	// complement
	for i := 0; i < len(seq); i++ {
		s[i] = comp(s[i])
	}
	return s
}

// complement a nucleotide sequence (standard version)
func comp(a byte) byte {
	switch a {
	case 'A', 'a':
		return 'T'
	case 'T', 't':
		return 'A'
	case 'G', 'g':
		return 'C'
	case 'C', 'c':
		return 'G'
	}
	return 0
}
