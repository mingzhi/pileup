package pileup

import (
	"fmt"
	"strconv"
	"strings"
)

type SNP struct {
	Ref     string
	Base    byte
	Pos     int
	Alleles []Allele
}

func (s *SNP) String() string {
	bases := []byte{}
	quals := []byte{}
	readIDs := []int{}
	mapQs := []byte{}
	for i := range s.Alleles {
		a := s.Alleles[i]
		bases = append(bases, a.Base)
		quals = append(quals, a.Qual)
		readIDs = append(readIDs, a.ReadID)
		mapQs = append(mapQs, a.MapQ)
	}

	readIDStrs := []string{}
	for i := range readIDStrs {
		readIDStrs = append(readIDStrs, fmt.Sprintf("%d", readIDs[i]))
	}

	return fmt.Sprintf("%s\t%d\t%c\t%d\t%s\t%s\t%v\t%s", s.Ref, s.Pos+1, s.Base, len(bases), bases, quals, readIDs, mapQs)
}

type Allele struct {
	Base   byte
	Qual   byte
	ReadID int
	MapQ   byte
}

func (a Allele) String() string {
	return fmt.Sprintf("Base: %c, Qual: %v, ReadID: %v, MapQ: %v", a.Base, a.Qual, a.ReadID, a.MapQ)
}

func Parse(line string) *SNP {
	var s SNP
	terms := strings.Split(strings.TrimSpace(line), "\t")
	s.Ref = terms[0]
	s.Pos = atoi(terms[1]) - 1
	s.Base = terms[2][0]
	bases := terms[3]
	quals := terms[4]

	var readIDs []int
	var mapQs string
	if len(terms) == 7 && terms[5][0] == '[' {
		readIDs = decodeInts(terms[5])
		mapQs = terms[6]
	} else if len(terms) == 6 {
		mapQs = terms[5]
	}

	for i := range terms[3] {
		a := Allele{
			Base: bases[i],
			Qual: quals[i],
		}

		if len(readIDs) == len(bases) {
			a.ReadID = readIDs[i]
		}

		if len(mapQs) == len(bases) {
			a.MapQ = mapQs[i]
		}

		s.Alleles = append(s.Alleles, a)
	}

	return &s
}

func decodeInts(s string) []int {
	values := []int{}
	l := len(s)
	s = s[1 : l-1]
	if len(s) > 0 {
		terms := strings.Split(s, " ")
		for i := range terms {
			values = append(values, atoi(terms[i]))
		}
	}

	return values
}

func decodeBytes(s string) []byte {
	values := []byte{}
	l := len(s)
	s = s[1 : l-1]
	if len(s) > 0 {
		terms := strings.Split(s, " ")
		for i := range terms {
			values = append(values, byte(atoi(terms[i])))
		}
	}

	return values
}

func atoi(s string) int {
	v, err := strconv.Atoi(s)
	if err != nil {
		panic(err)
	}
	return v
}
