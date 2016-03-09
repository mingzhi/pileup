package main

import (
	"fmt"
	"github.com/bmatsuo/lmdb-go/lmdb"
	"gopkg.in/vmihailenco/msgpack.v2"
	"log"
	"os"
)

type cmdReport2 struct {
	prefix    string
	featureDB *lmdb.Env
	resultsDB *lmdb.Env
}

func (c *cmdReport2) run() {
	geneSnpChan := getAllSNPs(c.resultsDB)
	geneResChan := c.getFeatures(geneSnpChan)
	c.groupGr(geneResChan)
}

func (c *cmdReport2) getFeatures(geneSnpChan chan SNPArr) chan GeneRes {
	ch := make(chan GeneRes)
	fn := func(txn *lmdb.Txn) error {
		dbi, err := txn.OpenDBI("feature", 0)
		if err != nil {
			return err
		}

		for gs := range geneSnpChan {
			k := gs.Key
			values := []float64{}
			for _, snp := range gs.Arr {
				values = append(values, float64(len(snp.Bases)))
			}
			v, err := txn.Get(dbi, k)
			if err != nil {
				return err
			}
			f := Feature{}
			if err := msgpack.Unmarshal(v, &f); err != nil {
				return err
			}
			ch <- GeneRes{Key: k, Values: values, Feature: f}
		}
		return nil
	}
	go func() {
		defer close(ch)
		err := c.featureDB.View(fn)
		if err != nil {
			log.Panicln(err)
		}
	}()

	return ch
}

func (c *cmdReport2) groupGr(grChan chan GeneRes) {
	w, err := os.Create(c.prefix + ".detectable.gene.csv")
	if err != nil {
		log.Fatalln(err)
	}
	defer w.Close()
	w.WriteString("patric_id,genome,figfam,sample,depth\n")

	for gr := range grChan {
		values := gr.Values
		total := 0.0
		for _, v := range values {
			total += v
		}
		m := total / float64(len(values))

		w.WriteString(fmt.Sprintf("%s,%s,%s,%s,%g\n",
			gr.Feature.PatricID,
			gr.Feature.TaxID,
			gr.Feature.FigfamID,
			c.prefix,
			m))
	}
}
