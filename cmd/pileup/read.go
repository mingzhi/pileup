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

type Pi struct {
	PatricId string
	Position int
	Pi       float64
	Coverage int
}

type extPi struct {
	pi         Pi
	start, end int
}

func (c *cmdRead) run() {
	c.openDB()
	// open pileup file
	snpChan := readPileup(c.pileupFile)
	piChan := calcPi(snpChan, c.db, c.minDepth)
	c.db.Update(createBucket("pi"))
	gc := collectPi(piChan, c.minCover)
	loadPi(gc, c.db)

	// check number of written records.
	c.db.View(func(tx *bolt.Tx) error {
		b := tx.Bucket([]byte("pi"))
		s := b.Stats()
		log.Printf("Wrote %d records\n", s.KeyN)
		return nil
	})
}

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

type Features []*Feature

func (s Features) Len() int      { return len(s) }
func (s Features) Swap(i, j int) { s[i], s[j] = s[j], s[i] }

type ByStart struct{ Features }
type ByEnd struct{ Features }

func (s ByStart) Less(i, j int) bool { return s.Features[i].Start < s.Features[j].Start }
func (s ByEnd) Less(i, j int) bool   { return s.Features[i].End < s.Features[j].End }

func calcPi(snpChan chan *pileup.SNP, db *bolt.DB, minDepth int) chan *extPi {
	c := make(chan *extPi, 10)
	go func() {
		defer close(c)
		currentGenome := ""
		features := []*Feature{}
		for snp := range snpChan {
			genome := snp.Reference
			if genome != currentGenome {
				features = []*Feature{}
				currentGenome = genome
				acc := genome
				if strings.Contains(genome, "|") {
					acc = strings.Split(genome, "|")[1]
				}
				acc = strings.Split(acc, ".")[0]
				fn := func(tx *bolt.Tx) error {
					b := tx.Bucket([]byte("genome"))
					v := b.Get([]byte(acc))
					if len(v) == 0 {
						return nil
					}
					ids := []string{}
					err := msgpack.Unmarshal(v, &ids)
					if err != nil {
						log.Panicln(err)
					}
					b = tx.Bucket([]byte("feature"))
					for _, id := range ids {
						if id == "" {
							continue
						}
						v := b.Get([]byte(id))
						f := Feature{}
						if err := msgpack.Unmarshal(v, &f); err != nil {
							log.Panicln(err)
						}
						if f.Type == "CDS" {
							features = append(features, &f)
						}
					}

					return nil
				}

				db.View(fn)
				if len(features) > 0 {
					sort.Sort(ByEnd{features})
				}
				if *debug {
					log.Println(acc)
				}
			}

			if len(snp.Bases) >= minDepth {
				f := findFeature(snp.Position, features)
				if f != nil {
					pi := Pi{
						PatricId: f.PatricID,
						Position: snp.Position,
						Pi:       snp.Pi(),
						Coverage: len(snp.Bases),
					}
					c <- &extPi{pi: pi, start: f.Start, end: f.End}
				}
			}
		}
	}()

	return c
}

func createPiBucket(db *bolt.DB) {
	db.Update(func(tx *bolt.Tx) error {
		_, err := tx.CreateBucket([]byte("pi"))
		if err != nil {
			if err == bolt.ErrBucketExists {
				tx.DeleteBucket([]byte("pi"))
				tx.CreateBucket([]byte("pi"))
			} else {
				log.Panicln(err)
			}
		}

		return nil
	})
}

type gPiCol struct {
	gene string
	col  []Pi
}

func collectPi(piChan chan *extPi, minCover float64) chan gPiCol {
	c := make(chan gPiCol, 100)
	go func() {
		defer close(c)
		geneId := ""
		buf := []Pi{}
		geneLen := 0
		for pi := range piChan {
			if pi.pi.PatricId != geneId {
				if float64(len(buf))/float64(geneLen) >= minCover {
					c <- gPiCol{gene: geneId, col: buf}
				}
				geneId = pi.pi.PatricId
				buf = []Pi{}
			}
			buf = append(buf, pi.pi)
			geneLen = pi.end - pi.start + 1
		}
	}()
	return c
}

func loadPi(c chan gPiCol, db *bolt.DB) {
	buffer := []gPiCol{}
	for gc := range c {
		if len(buffer) > 1000 {
			db.Update(func(tx *bolt.Tx) error {
				for _, gc := range buffer {
					geneId := gc.gene
					buf := gc.col
					key := []byte(geneId)
					value, err := msgpack.Marshal(buf)
					if err != nil {
						log.Panicln(err)
					}
					if err := tx.Bucket([]byte("pi")).Put(key, value); err != nil {
						log.Panicln(err)
					}
				}

				return nil
			})
			buffer = []gPiCol{}
		}
		buffer = append(buffer, gc)
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
