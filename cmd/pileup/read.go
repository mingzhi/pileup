package main

import (
	"github.com/boltdb/bolt"
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

	db *bolt.DB
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
	c.openDB()
	// open pileup file and read SNP.
	snpChan := readPileup(c.pileupFile)
	// group SNPs into genes.
	geneChan := groupSNPs(snpChan, c.db, c.minDepth, c.minCover)
	// create gene bucket.
	c.db.Update(createBucket("gene"))
	// load gene into the bucket.
	loadGenes(geneChan, c.db, "gene")

	// check number of written records.
	c.db.View(func(tx *bolt.Tx) error {
		b := tx.Bucket([]byte("gene"))
		s := b.Stats()
		log.Printf("Wrote %d records\n", s.KeyN)
		return nil
	})
}

// open a bolt db.
func (c *cmdRead) openDB() {
	if _, err := os.Stat(c.dbfile); os.IsNotExist(err) {
		log.Fatalf("can not find feature db: %s\n", c.dbfile)
	}
	db, err := bolt.Open(c.dbfile, 0600, nil)
	if err != nil {
		log.Fatal(err)
	}
	c.db = db
}

// group SNPs into genes.
func groupSNPs(snpChan chan *pileup.SNP, db *bolt.DB, minDepth int, minCover float64) chan *Gene {
	c := make(chan *Gene, 100)
	go func() {
		defer close(c)
		currentGenome := Genome{}
		var currentGene *Gene
		var toUpdate bool

		for snp := range snpChan {
			// update genome features.
			reference := cleanAccession(snp.Reference)
			if reference != currentGenome.Reference {
				currentGenome = queryGenome(db, reference)
			}

			if len(snp.Bases) >= minDepth {
				toUpdate = false
				if currentGene == nil {
					toUpdate = true
				} else {
					outBound := snp.Position < currentGene.Start || snp.Position > currentGene.End
					if outBound {
						toUpdate = true
						geneLen := currentGene.End - currentGene.Start + 1
						if float64(len(currentGene.SNPs))/float64(geneLen) >= minCover {
							c <- currentGene
						}
					}
				}

				if toUpdate {
					f := findFeature(snp.Position, currentGenome.Features)
					if f != nil {
						currentGene = &Gene{}
						currentGene.Feature = f
					} else {
						currentGene = nil
					}
				}

				if currentGene != nil {
					currentGene.SNPs = append(currentGene.SNPs, snp)
				}
			}
		}
	}()
	return c
}

// query genome from bolt db.
func queryGenome(db *bolt.DB, reference string) Genome {
	features := []*Feature{}
	fn := func(tx *bolt.Tx) error {
		// first, we query genome bucket to
		// obtain a list of gene IDs.
		b := tx.Bucket([]byte("genome"))
		v := b.Get([]byte(reference))
		// check if genome is found.
		if len(v) == 0 {
			return nil
		}
		// unpack gene IDs.
		geneIDs := []string{}
		err := msgpack.Unmarshal(v, &geneIDs)
		if err != nil {
			log.Panicln(err)
		}

		// now, we try obtain gene feature informations
		// from feature bucket.
		b = tx.Bucket([]byte("feature"))
		for _, id := range geneIDs {
			v := b.Get([]byte(id))
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

	db.View(fn)

	// sort features.
	if len(features) > 0 {
		sort.Sort(ByEnd{features})
	}

	return Genome{Reference: reference, Features: features}
}

func cleanAccession(reference string) string {
	if strings.Contains(reference, "|") {
		reference = strings.Split(reference, "|")[1]
	}
	return strings.Split(reference, ".")[0]
}

// load genes into db.
func loadGenes(geneChan chan *Gene, db *bolt.DB, bucketName string) {
	bufferSize := 1000
	buffer := []*Gene{}
	for gene := range geneChan {
		if len(buffer) > bufferSize {
			fn := func(tx *bolt.Tx) error {
				for _, g := range buffer {
					key := []byte(g.PatricID)
					value, err := msgpack.Marshal(g.SNPs)
					if err != nil {
						log.Panicln(err)
					}

					err = tx.Bucket([]byte(bucketName)).Put(key, value)
					if err != nil {
						log.Panicln(err)
					}
				}
				return nil
			}
			db.Update(fn)
			buffer = []*Gene{}
		}
		buffer = append(buffer, gene)
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
