package main

import (
	"fmt"
	"github.com/boltdb/bolt"
	"github.com/mingzhi/biogo/pileup"
	"github.com/mingzhi/gomath/stat/desc/meanvar"
	"gopkg.in/vmihailenco/msgpack.v2"
	"log"
	"math"
	"os"
	"strings"
)

type cmdReport struct {
	dbfile string
	prefix string

	db *bolt.DB
}

func (c *cmdReport) run() {
	c.openDB()
	// obtain gene coverage map.
	geneInforMap := c.getGeneInfors()
	log.Printf("Total gene number: %d\n", len(geneInforMap))

	c.reportTax(geneInforMap)
	c.reportPathway(geneInforMap)
}

func (c *cmdReport) reportTax(geneInforMap map[string]GeneInfor) {
	// group genes in tax level.
	taxGenes := make(map[string][]GeneInfor)
	for _, f := range geneInforMap {
		taxID := f.TaxID
		taxGenes[taxID] = append(taxGenes[taxID], f)
	}
	// report
	taxOut := c.prefix + ".tax.csv"
	w, err := os.Create(taxOut)
	if err != nil {
		log.Fatal(err)
	}
	defer w.Close()

	crOut := c.prefix + ".tax.cr.csv"
	crW, err := os.Create(crOut)
	if err != nil {
		log.Fatal(err)
	}
	defer crW.Close()

	geneOut := c.prefix + ".tax.coverage.csv"
	geneW, err := os.Create(geneOut)
	if err != nil {
		log.Fatal(err)
	}
	defer geneW.Close()

	for taxID, genes := range taxGenes {
		w.WriteString(fmt.Sprintf("%s,%d\n", taxID, len(genes)))

		if len(genes) > 100 {
			maxl := 300
			mvs := make([]*meanvar.MeanVar, maxl)
			for i := 0; i < maxl; i++ {
				mvs[i] = meanvar.New()
			}

			c.db.View(func(tx *bolt.Tx) error {
				for _, gene := range genes {
					k := []byte(gene.PatricID)
					b := tx.Bucket([]byte("cr"))
					v := b.Get(k)
					if len(v) == 0 {
						log.Println(gene.PatricID)
						continue
					}
					cc := []float64{}
					if err := msgpack.Unmarshal(v, &cc); err != nil {
						log.Panicln(err)
					}
					for i := 0; i < len(cc) && i < maxl; i++ {
						if !math.IsNaN(cc[i]) {
							mvs[i].Increment(cc[i])
						}
					}
				}
				return nil
			})

			for i := 0; i < maxl; i++ {
				m := mvs[i].Mean.GetResult()
				v := mvs[i].Var.GetResult()
				n := mvs[i].Var.GetN()
				if !math.IsNaN(v) {
					crW.WriteString(fmt.Sprintf("%s,%d,%g,%g,%d\n", taxID, i, m, v, n))
				}
			}

			for _, g := range genes {
				geneW.WriteString(fmt.Sprintf("%s,%s,%g,%g,%g,%d\n", taxID, g.PatricID, g.Coverage, g.Depth, g.Pi, (g.End - g.Start + 1)))
			}
		}
	}
}

type Pathway struct {
	ID, Name string
}

func (c *cmdReport) reportPathway(geneInforMap map[string]GeneInfor) {
	pathwayMap := make(map[string]Pathway)
	pathwayGenes := make(map[string][]GeneInfor)
	for _, f := range geneInforMap {
		terms := strings.Split(f.Pathway, ";")
		for _, t := range terms {
			fields := strings.Split(t, "|")
			if len(fields) < 2 {
				continue
			}
			pathwayID, pathwayName := fields[0], fields[1]
			_, found := pathwayMap[pathwayID]
			if !found {
				pathwayMap[pathwayID] = Pathway{ID: pathwayID, Name: pathwayName}
			}
			pathwayGenes[pathwayID] = append(pathwayGenes[pathwayID], f)
		}
	}

	// report

	out := c.prefix + ".pathway.csv"
	w, err := os.Create(out)
	if err != nil {
		log.Fatal(err)
	}
	defer w.Close()

	crOut := c.prefix + ".pathway.cr.csv"
	crW, err := os.Create(crOut)
	if err != nil {
		log.Fatal(err)
	}
	defer crW.Close()

	geneOut := c.prefix + ".pathway.coverage.csv"
	geneW, err := os.Create(geneOut)
	if err != nil {
		log.Fatal(err)
	}
	defer geneW.Close()

	for id, genes := range pathwayGenes {
		name := pathwayMap[id].Name
		w.WriteString(fmt.Sprintf("%s,\"%s\",%d\n", id, name, len(genes)))

		if len(genes) > 4 {
			maxl := 300
			mvs := make([]*meanvar.MeanVar, maxl)
			for i := 0; i < maxl; i++ {
				mvs[i] = meanvar.New()
			}

			c.db.View(func(tx *bolt.Tx) error {
				for _, gene := range genes {
					k := []byte(gene.PatricID)
					b := tx.Bucket([]byte("cr"))
					v := b.Get(k)
					if len(v) == 0 {
						log.Println(gene.PatricID)
						continue
					}
					cc := []float64{}
					if err := msgpack.Unmarshal(v, &cc); err != nil {
						log.Panicln(err)
					}
					for i := 0; i < len(cc) && i < maxl; i++ {
						if !math.IsNaN(cc[i]) {
							mvs[i].Increment(cc[i])
						}
					}
				}
				return nil
			})

			for i := 0; i < maxl; i++ {
				m := mvs[i].Mean.GetResult()
				v := mvs[i].Var.GetResult()
				n := mvs[i].Var.GetN()
				if !math.IsNaN(v) {
					crW.WriteString(fmt.Sprintf("%s,%d,%g,%g,%d\n", id, i, m, v, n))
				}
			}

			for _, g := range genes {
				geneW.WriteString(fmt.Sprintf("%s,%s,%g,%g,%g,%d\n", id, g.PatricID, g.Coverage, g.Depth, g.Pi, (g.End - g.Start + 1)))
			}
		}
	}
}

func (c *cmdReport) openDB() {
	if _, err := os.Stat(c.dbfile); os.IsNotExist(err) {
		log.Fatalf("can not open db: %s\n", c.dbfile)
	}
	db, err := bolt.Open(c.dbfile, 0600, nil)
	if err != nil {
		log.Fatalln(err)
	}
	c.db = db
}

type GeneInfor struct {
	Feature
	Coverage, Depth, Pi float64
}

// return a map of gene to its coverage.
// coverage is calculated as the size of
// positions that are coveraged by at least min depth (default 5).
func (c *cmdReport) getGeneInfors() map[string]GeneInfor {
	geneMap := make(map[string]GeneInfor)
	c.db.View(func(tx *bolt.Tx) error {
		// assume bucket exists and has keys.
		b := tx.Bucket([]byte("gene"))
		featureBucket := tx.Bucket([]byte("feature"))
		b.ForEach(func(k, v []byte) error {
			geneID := string(k)
			snpArr := []pileup.SNP{}
			if err := msgpack.Unmarshal(v, &snpArr); err != nil {
				log.Panicln(err)
			}

			// get feature.
			feature := Feature{}
			featureBytes := featureBucket.Get(k)
			if err := msgpack.Unmarshal(featureBytes, &feature); err != nil {
				log.Panicln(err)
			}
			coverage := float64(len(snpArr)) / float64(feature.End-feature.Start+1)
			depth := 0.0
			pi := 0.0
			for _, s := range snpArr {
				depth += float64(len(s.Bases))
				pi += s.Pi()
			}
			depth /= float64(len(snpArr))
			pi /= float64(len(snpArr))
			gi := GeneInfor{Feature: feature, Coverage: coverage, Depth: depth, Pi: pi}

			geneMap[geneID] = gi
			return nil
		})

		return nil
	})

	return geneMap
}

func (c *cmdReport) getFeatures(geneIDs []string) map[string]Feature {
	m := make(map[string]Feature)
	c.db.View(func(tx *bolt.Tx) error {
		b := tx.Bucket([]byte("feature"))
		for _, id := range geneIDs {
			v := b.Get([]byte(id))
			if len(v) > 0 {
				f := Feature{}
				if err := msgpack.Unmarshal(v, &f); err != nil {
					log.Panicln(err)
				}
				m[id] = f
			}
		}

		return nil
	})
	return m
}
