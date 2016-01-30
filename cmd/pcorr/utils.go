package main

import (
	"encoding/json"
	"github.com/mingzhi/biogo/feat/gff"
	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/ncbiftp/genomes/profiling"
	"log"
	"os"
)

func openFile(filename string) *os.File {
	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}

	return f
}

func createFile(filename string) *os.File {
	w, err := os.Create(filename)
	if err != nil {
		panic(err)
	}

	return w
}

// readGenome read the genome file
// and return a genomic sequence.
func readGenome(filename string) []byte {
	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	rd := seq.NewFastaReader(f)
	ss, err := rd.ReadAll()
	if err != nil {
		panic(err)
	}

	return ss[0].Seq
}

func readGff(filename string) []*gff.Record {
	f := openFile(filename)
	defer f.Close()

	rd := gff.NewReader(f)
	ss, err := rd.ReadAll()
	if err != nil {
		panic(err)
	}

	records := []*gff.Record{}
	for _, s := range ss {
		if s.Feature == "CDS" {
			records = append(records, s)
		}
	}

	return records
}

func readPi(filename string) chan Pi {
	c := make(chan Pi)
	go func() {
		defer close(c)

		f := openFile(filename)
		defer f.Close()
		decoder := json.NewDecoder(f)

		for decoder.More() {
			var pi Pi
			err := decoder.Decode(&pi)
			if err != nil {
				log.Fatal(err)
			}
			c <- pi
		}
	}()

	return c
}

func checkPosType(posType, t1 byte) bool {
	isFirstPos := t1 == profiling.FirstPos
	isSecondPos := t1 == profiling.SecondPos
	isThirdPos := t1 == profiling.ThirdPos
	isFourFold := t1 == profiling.FourFold

	if posType == profiling.Coding {
		if isFirstPos || isSecondPos || isThirdPos || isFourFold {
			return true
		}
		return false
	}

	if posType == profiling.ThirdPos {
		if isThirdPos || isFourFold {
			return true
		}
		return false
	}

	return posType == t1
}

func convertPosType(pos int) byte {
	var p byte
	switch pos {
	case 0:
		p = profiling.NonCoding
	case 1:
		p = profiling.FirstPos
		break
	case 2:
		p = profiling.SecondPos
		break
	case 3:
		p = profiling.ThirdPos
		break
	case 4:
		p = profiling.FourFold
		break
	default:
		p = profiling.Coding
	}

	return p
}

func isATGC(b byte) bool {
	if b == 'A' {
		return true
	} else if b == 'T' {
		return true
	} else if b == 'C' {
		return true
	} else if b == 'G' {
		return true
	}

	return false
}
