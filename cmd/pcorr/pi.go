package main

import (
	"bytes"
	"encoding/json"
	"log"
	"os"
)

type cmdPi struct {
	debug                  bool
	pileupFile             string
	outFile                string
	minBQ                  int
	regionStart, regionEnd int
}

type Pi struct {
	Ref     string
	Base    string
	Pos     int
	Alleles map[string]int
}

func (p Pi) Pi() (pi float64) {
	total := 0
	nums := []int{}
	for _, n := range p.Alleles {
		total += n
		nums = append(nums, n)
	}

	cross := 0
	for i := 0; i < len(nums); i++ {
		for j := i + 1; j < len(nums); j++ {
			cross += nums[i] * nums[j]
		}
	}

	pi = float64(cross) / float64(total*(total-1)/2)

	return pi
}

func (c *cmdPi) Run() {
	f := openFile(c.pileupFile)
	defer f.Close()
	snpChan := readPileup(f, 0, 0)
	piChan := make(chan Pi)
	go func() {
		defer close(piChan)
		for s := range snpChan {
			if c.regionEnd > 0 && s.Pos > c.regionEnd {
				break
			}
			if s.Pos >= c.regionStart {
				if len(s.Alleles) > 0 {
					bases := []byte{}
					quals := []byte{}
					for _, a := range s.Alleles {
						bases = append(bases, a.Base)
						quals = append(quals, a.Qual)
					}
					bases, _ = c.filterBases(bases, quals)
					bases = bytes.ToUpper(bases)

					m := make(map[string]int)
					for _, b := range bases {
						m[string(b)]++
					}

					pi := Pi{
						Ref:     s.Ref,
						Base:    string(s.Base),
						Pos:     s.Pos,
						Alleles: m,
					}

					if c.debug {
						log.Println(string(bases))
						log.Println(pi)
					}

					piChan <- pi
				}
			}
		}
	}()

	var w *os.File
	if c.outFile != "" {
		w = createFile(c.outFile)
	} else {
		w = os.Stdout
	}
	defer w.Close()

	encoder := json.NewEncoder(w)

	for pi := range piChan {
		if err := encoder.Encode(pi); err != nil {
			log.Fatalln(err)
		}
	}
}

func (c *cmdPi) filterBases(bases []byte, quals []byte) (bases1, quals1 []byte) {
	for i := range bases {
		if int(quals[i]) > c.minBQ {
			bases1 = append(bases1, bases[i])
			quals1 = append(quals1, quals[i])
		}
	}

	return
}
