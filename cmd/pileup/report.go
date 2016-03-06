package main

import (
	"fmt"
	"github.com/bmatsuo/lmdb-go/lmdb"
	"github.com/mingzhi/gomath/stat/desc/meanvar"
	"gopkg.in/vmihailenco/msgpack.v2"
	"log"
	"math"
	"os"
	"strings"
)

type cmdReport2 struct {
	prefix    string
	featureDB *lmdb.Env
	resultsDB *lmdb.Env
	maxl      int
}

func (c *cmdReport2) run() {
	crChan := c.getCr()
	grChan := c.getFeature(crChan)
	c.groupGr(grChan)
}

// Iterate all cr results,
// output a channel of CovRes.
func (c *cmdReport2) getCr() chan CovRes {
	ch := make(chan CovRes)
	fn := func(txn *lmdb.Txn) error {
		dbi, err := txn.OpenDBI("cr", 0)
		if err != nil {
			return err
		}

		cur, err := txn.OpenCursor(dbi)
		if err != nil {
			return err
		}
		defer cur.Close()

		for {
			k, v, err := cur.Get(nil, nil, lmdb.Next)
			if lmdb.IsNotFound(err) {
				break
			}

			if err != nil {
				return err
			}

			values := []float64{}
			if err := msgpack.Unmarshal(v, &values); err != nil {
				return err
			}

			ch <- CovRes{Key: k, Values: values}
		}

		return nil
	}

	go func() {
		defer close(ch)
		err := c.resultsDB.View(fn)
		if err != nil {
			log.Panicln(err)
		}
	}()
	return ch
}

// GeneRes contains covariances and its feature.
type GeneRes struct {
	Key     []byte
	Feature Feature
	Values  []float64
}

func (c *cmdReport2) getFeature(crChan chan CovRes) chan GeneRes {
	ch := make(chan GeneRes)
	fn := func(txn *lmdb.Txn) error {
		dbi, err := txn.OpenDBI("feature", 0)
		if err != nil {
			return err
		}

		for cr := range crChan {
			v, err := txn.Get(dbi, cr.Key)
			if err != nil {
				return err
			}
			f := Feature{}
			if err := msgpack.Unmarshal(v, &f); err != nil {
				return err
			}
			ch <- GeneRes{Key: cr.Key, Values: cr.Values, Feature: f}
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
	gm := make(map[string][]*meanvar.MeanVar)
	pm := make(map[string][]*meanvar.MeanVar)
	for gr := range grChan {
		for i := 0; i < len(gr.Values) && i < c.maxl; i++ {
			v := gr.Values[i]
			if !math.IsNaN(v) {
				g := gr.Feature.TaxID

				_, found := gm[g]
				if !found {
					gm[g] = newMeanVars(c.maxl)
				}
				gm[g][i].Increment(v)

				pathways := gr.Feature.Pathway
				terms := strings.Split(pathways, ";")
				for _, s := range terms {
					id := strings.Split(s, "|")[0]
					_, found := pm[id]
					if !found {
						pm[id] = newMeanVars(c.maxl)
					}
					pm[id][i].Increment(v)
				}
			}
		}
	}

	wfn := func(filename string, m map[string][]*meanvar.MeanVar) {
		w, err := os.Create(filename)
		if err != nil {
			log.Panicln(err)
		}
		defer w.Close()

		w.WriteString("id,l,m,v,n\n")
		for id, mvs := range m {
			for i := range mvs {
				m := mvs[i].Mean.GetResult()
				v := mvs[i].Var.GetResult()
				n := mvs[i].Mean.GetN()
				if !math.IsNaN(v) {
					w.WriteString(fmt.Sprintf("%s,%d,%g,%g,%d\n", id, i, m, v, n))
				}
			}
		}
	}

	wfn(c.prefix+".genome.cr.csv", gm)
	wfn(c.prefix+".pathway.cr.csv", pm)
}
