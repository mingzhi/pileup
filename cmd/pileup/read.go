package main

import (
	"github.com/bmatsuo/lmdb-go/lmdb"
	"github.com/mingzhi/biogo/pileup"
	"gopkg.in/vmihailenco/msgpack.v2"
	"io"
	"log"
	"os"
	"sort"
	"strings"
)

type cmdRead struct {
	pileupFile string
	dbfile     string
	minDepth   int
	minCover   float64
	featureDB  string

	env        *lmdb.Env
	featureEnv *lmdb.Env
}

// Gene contains a group of SNP and its feature.
type Gene struct {
	*Feature
	SNPs []*pileup.SNP
}

// Genome is a group of features.
type Genome struct {
	Reference string
	Features  []*Feature
}

func (c *cmdRead) run() {
	// create environment and dbi.
	c.env = createEnv(c.dbfile)
	defer c.env.Close()

	c.featureEnv = createReadOnlyEnv(c.featureDB)
	defer c.featureEnv.Close()

	createDBI(c.env, "gene")
	// open pileup file and read SNP.
	snpChan := readPileup(c.pileupFile)
	// group SNPs into genes.
	geneChan := c.groupSNPs(snpChan)
	c.loadGenes(geneChan)

	err := c.env.View(func(tx *lmdb.Txn) error {
		dbi, err := tx.OpenDBI("gene", 0)
		if err != nil {
			return err
		}
		cur, err := tx.OpenCursor(dbi)
		if err != nil {
			return err
		}
		count := 0
		for {
			_, _, err := cur.Get(nil, nil, lmdb.Next)
			if lmdb.IsNotFound(err) {
				return nil
			}
			if err != nil {
				return err
			}
			count++
		}
		log.Printf("Total gene: %d\n", count)
		return nil
	})
	if err != nil {
		log.Panicln(err)
	}
}

func (c *cmdRead) groupSNPs(snpChan chan *pileup.SNP) chan *Gene {
	ch := make(chan *Gene, 100)
	go func() {
		defer close(ch)
		var currentGenome *Genome
		var currentGene *Gene
		var toUpdate bool

		for snp := range snpChan {
			// update genome features.
			reference := cleanAccession(snp.Reference)
			if currentGenome == nil || reference != currentGenome.Reference {
				currentGenome = c.queryGenome(reference)
			}

			if len(snp.Bases) >= c.minDepth {
				toUpdate = false
				if currentGene == nil {
					toUpdate = true
				} else {
					outBound := snp.Position < currentGene.Start || snp.Position > currentGene.End
					if outBound {
						toUpdate = true
						geneLen := currentGene.End - currentGene.Start + 1
						if float64(len(currentGene.SNPs))/float64(geneLen) >= c.minCover {
							ch <- currentGene
						}
					}
				}

				if toUpdate {
					currentGene = nil
					if currentGenome != nil {
						f := findFeature(snp.Position, currentGenome.Features)
						if f != nil {
							currentGene = &Gene{}
							currentGene.Feature = f
						}
					}
				}

				if currentGene != nil {
					currentGene.SNPs = append(currentGene.SNPs, snp)
				}
			}
		}
	}()
	return ch
}

func (c *cmdRead) queryGenome(reference string) *Genome {
	features := []*Feature{}
	fn := func(tx *lmdb.Txn) error {
		// first, we query genome bucket to
		// obtain a list of gene IDs.
		dbi, err := tx.OpenDBI("genome", 0)
		if err != nil {
			return err
		}
		v, err := tx.Get(dbi, []byte(reference))
		if err != nil {
			return err
		}
		// unpack gene IDs.
		geneIDs := []string{}
		err = msgpack.Unmarshal(v, &geneIDs)
		if err != nil {
			log.Panicln(err)
		}

		// now, we try obtain gene feature informations
		// from feature bucket.
		dbi, err = tx.OpenDBI("feature", 0)
		for _, id := range geneIDs {
			if len(id) == 0 {
				continue
			}
			v, err := tx.Get(dbi, []byte(id))
			if err != nil {
				return err
			}
			if len(v) > 0 {
				f := Feature{}
				if err := msgpack.Unmarshal(v, &f); err != nil {
					log.Panicln(err)
				}
				if f.Type == "CDS" || strings.Contains(f.Type, "RNA") {
					features = append(features, &f)
				}
			}
		}

		return nil
	}

	err := c.featureEnv.View(fn)
	if err != nil {
		if lmdb.IsNotFound(err) {
			if *debug {
				log.Printf("key %s is not found\n", reference)
			}
			return nil
		}
		log.Panicln(err)
	}

	// sort features.
	if len(features) > 0 {
		sort.Sort(ByEnd{features})
	}

	return &Genome{Reference: reference, Features: features}
}

func cleanAccession(reference string) string {
	if strings.Contains(reference, "|") {
		reference = strings.Split(reference, "|")[1]
	}
	return strings.Split(reference, ".")[0]
}

func (c *cmdRead) loadGenes(geneChan chan *Gene) {
	bufferSize := 100
	buffer := []*Gene{}
	for gene := range geneChan {
		if len(buffer) > bufferSize {
			c.loadBuffer(buffer)
			buffer = []*Gene{}
		}
		buffer = append(buffer, gene)
	}
	c.loadBuffer(buffer)
}

func (c *cmdRead) loadBuffer(buffer []*Gene) {
	fn := func(txn *lmdb.Txn) error {
		dbi, err := txn.OpenDBI("gene", 0)
		if err != nil {
			return err
		}
		for _, g := range buffer {
			key := []byte(g.PatricID)
			value, err := msgpack.Marshal(g.SNPs)
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
	if err != nil {
		log.Panicln(err)
	}
}

func findFeature(pos int, features []*Feature) *Feature {
	if len(features) == 0 {
		return nil
	}
	i := sort.Search(len(features), func(i int) bool { return features[i].End >= pos })
	if i >= len(features) {
		return nil
	}

	if features[i].Start <= pos && features[i].End >= pos {
		return features[i]
	}

	return nil
}

func readPileup(pileupFile string) chan *pileup.SNP {
	c := make(chan *pileup.SNP, 10)
	go func() {
		defer close(c)

		f, err := os.Open(pileupFile)
		if err != nil {
			log.Fatalln(err)
		}
		defer f.Close()

		pileupReader := pileup.NewReader(f)
		for {
			snp, err := pileupReader.Read()
			if err != nil {
				if err != io.EOF {
					log.Fatalln(err)
				}
				break
			}
			c <- snp
		}
	}()
	return c
}
