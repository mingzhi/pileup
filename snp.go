package pileup

import (
	"bytes"
	"fmt"
	"regexp"
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
	bases := decodeReadBases(terms[3], s.Base)
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

func decodeReadBases(s string, ref byte) []byte {
	r := regexp.MustCompile("\\^.")
	s = r.ReplaceAllString(s, "")
	s = strings.Replace(s, "$", "", -1)

	r2 := regexp.MustCompile("[\\+-][0-9]+")
	if r2.MatchString(s) {
		insertNumbers := r2.FindAllString(s, -1)
		insertIndex := r2.FindAllStringIndex(s, -1)
		deletedPositions := make(map[int]bool)
		for i := 0; i < len(insertNumbers); i++ {
			start := insertIndex[i][0]
			n := atoi(insertNumbers[i][1:])
			for j := start; j < start+2+n; j++ {
				deletedPositions[j] = true
			}
		}
		bs := []byte{}
		for i := 0; i < len(s); i++ {
			if !deletedPositions[i] {
				bs = append(bs, s[i])
			}
		}
		s = string(bs)
	}

	bases := []byte{}
	for i := 0; i < len(s); i++ {
		b := s[i]
		if b == '.' || b == ',' {
			bases = append(bases, ref)
		} else {
			bases = append(bases, b)
		}
	}

	bases = bytes.ToUpper(bases)

	return bases
}

func atoi(s string) int {
	v, err := strconv.Atoi(s)
	if err != nil {
		panic(err)
	}
	return v
}
