package pileup

import (
	"bytes"
	"errors"
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
	Num     int
}

func (s *SNP) String() string {
	bases := []byte{}
	quals := []byte{}
	QNames := []string{}
	for i := range s.Alleles {
		a := s.Alleles[i]
		bases = append(bases, a.Base)
		quals = append(quals, a.Qual)
		QNames = append(QNames, a.QName)
	}

	return fmt.Sprintf("%s\t%d\t%c\t%d\t%s\t%s\t%v\t", s.Ref, s.Pos+1, s.Base, len(bases), bases, quals, strings.Join(QNames, ","))
}

type Allele struct {
	Base  byte
	Qual  byte
	QName string // read name
}

func (a Allele) String() string {
	return fmt.Sprintf("Base: %c, Qual: %v, ReadID: %v", a.Base, a.Qual, a.QName)
}

func parse(line string) (*SNP, error) {
	var s SNP
	terms := strings.Split(strings.TrimSpace(line), "\t")
	s.Ref = terms[0]
	s.Pos = atoi(terms[1]) - 1
	s.Base = strings.ToUpper(terms[2])[0]
	s.Num = atoi(terms[3])
	if s.Num == 0 {
		return &s, nil
	}
	bases := decodeReadBases(terms[4], s.Base)
	quals := terms[5]

	var QNames []string
	if len(terms) == 7 {
		QNames = strings.Split(terms[6], ",")
	}

	// check bases and quals len.
	if len(bases) != s.Num {
		err := errors.New("Decode bases length did not match the indicated number.")
		return nil, err
	}

	if len(quals) != s.Num {
		err := errors.New("Decode quals length did not match the indicated number.")
		return nil, err
	}

	for i := range bases {
		a := Allele{
			Base: bases[i],
			Qual: quals[i],
		}
		s.Alleles = append(s.Alleles, a)
	}

	if len(QNames) == len(bases) {
		for i := range QNames {
			s.Alleles[i].QName = QNames[i]
		}
	}

	return &s, nil
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
			digit := insertNumbers[i][1:]
			n := atoi(digit)
			for j := start; j < start+len(digit)+n+1; j++ {
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
