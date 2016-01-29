package main

import (
	"bytes"
	"encoding/json"
	"os"
)

type cmdPi struct {
	pileupFile string
	outFile    string
	minBQ      int
}

type Pi struct {
	Ref     string
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
	snpChan := readPileup(f)
	piChan := make(chan Pi)
	go func() {
		defer close(piChan)
		for s := range snpChan {
			if len(s.Alleles) > 0 {
				bases := []byte{}
				quals := []byte{}
				for _, a := range s.Alleles {
					bases = append(bases, a.Base)
					quals = append(quals, a.Qual)
				}
				bases1, _ := filterBases(bases, quals, c.minBQ)
				bases1 = bytes.ToUpper(bases1)
				m := make(map[string]int)
				for _, b := range bases1 {
					m[string(b)]++
				}

				pi := Pi{
					Ref:     s.Ref,
					Pos:     s.Pos,
					Alleles: m,
				}
				piChan <- pi
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
		encoder.Encode(pi)
	}
}

func filterBases(bases []byte, quals []byte, cutoff int) (bases1, quals1 []byte) {
	for i := range bases {
		if int(quals[i]) > cutoff {
			bases1 = append(bases1, bases[i])
			quals1 = append(quals1, quals[i])
		}
	}

	return
}
