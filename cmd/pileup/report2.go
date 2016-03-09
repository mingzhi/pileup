package main

import (
	"fmt"
	"github.com/bmatsuo/lmdb-go/lmdb"
	"github.com/montanaflynn/stats"
	"gopkg.in/vmihailenco/msgpack.v2"
	"log"
	"os"
	"sort"
)

type cmdReport2 struct {
	prefix    string
	featureDB *lmdb.Env
	resultsDB *lmdb.Env
}

func (c *cmdReport2) run() {
	geneSnpChan := getAllSNPs(c.resultsDB)
	c.getFeatures(geneSnpChan)
}

func (c *cmdReport2) getFeatures(geneSnpChan chan SNPArr) {
	w, err := os.Create(c.prefix + ".detectable.gene.csv")
	if err != nil {
		log.Fatalln(err)
	}
	defer w.Close()
	w.WriteString("patric_id,genome,figfam,sample,pi,depth\n")

	fn := func(txn *lmdb.Txn) error {
		dbi, err := txn.OpenDBI("feature", 0)
		if err != nil {
			return err
		}

		for gs := range geneSnpChan {
			if len(gs.Arr) < 100 {
				continue
			}

			k := gs.Key
			v, err := txn.Get(dbi, k)
			if err != nil {
				return err
			}
			f := Feature{}
			if err := msgpack.Unmarshal(v, &f); err != nil {
				return err
			}

			seqLen := f.End - f.Start + 1

			// calculate median of depth
			depthArr := []float64{}
			piArr := []float64{}
			for _, snp := range gs.Arr {
				pos := snp.Position - f.Start
				if f.IsComplementaryStrand() {
					pos = seqLen - 1 - pos
				}
				if (pos+1)%3 == 0 {
					depthArr = append(depthArr, float64(len(snp.Bases)))
					piArr = append(piArr, snp.Pi())
				}

			}
			depthMedian, _ := stats.Median(depthArr)
			sort.Float64s(piArr)
			piMean, _ := stats.Mean(piArr[10 : len(piArr)-10])

			w.WriteString(fmt.Sprintf("%s,%s,%s,%s,%g,%g\n",
				f.PatricID,
				f.TaxID,
				f.FigfamID,
				c.prefix,
				piMean,
				depthMedian))
		}
		return nil
	}

	err = c.featureDB.View(fn)
	if err != nil {
		log.Panicln(err)
	}
}
