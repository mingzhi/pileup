package main

import (
	"bytes"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/mingzhi/biogo/seq"
	"io"
	"log"
	"strings"
)

type cmdFilter struct {
	fnaFile     string
	bamFile     string
	outFile     string
	maxDistance float64
	mapQ        int

	filters []filter
}

type filter interface {
	Filter(r *sam.Record) bool
}

func (c *cmdFilter) run() {
	c.filters = append(c.filters, &MapQFilter{Cutoff: c.mapQ})
	c.filters = append(c.filters, &ProperPairFilter{})
	c.filters = append(c.filters, &ProperMapFilter{})
	c.filters = append(c.filters, &DiversityFilter{RefMap: readRefSeq(c.fnaFile), Cutoff: c.maxDistance})
	header := c.readBamHeader()
	ch := c.readBamFile()
	filteredChan := c.filter(ch)
	c.write(header, filteredChan)
}

func (c *cmdFilter) readBamHeader() *sam.Header {
	f := openFile(c.bamFile)
	defer f.Close()
	rd, err := bam.NewReader(f, 0)
	raiseError(err)
	defer rd.Close()
	header := rd.Header()
	return header
}

func (c *cmdFilter) readBamFile() (ch chan *sam.Record) {
	ch = make(chan *sam.Record)

	go func() {
		defer close(ch)
		f := openFile(c.bamFile)
		defer f.Close()
		rd, err := bam.NewReader(f, 0)
		raiseError(err)
		defer rd.Close()
		for {
			record, err := rd.Read()
			if err != nil {
				if err != io.EOF {
					raiseError(err)
				}
				break
			}
			ch <- record
		}
	}()

	return
}

func (c *cmdFilter) filter(ch chan *sam.Record) chan *sam.Record {
	filteredChan := make(chan *sam.Record)
	go func() {
		defer close(filteredChan)
		for r := range ch {
			good := true
			for _, f := range c.filters {
				good = f.Filter(r)
				if !good {
					break
				}
			}

			if good {
				filteredChan <- r
			}
		}
	}()
	return filteredChan
}

func (c *cmdFilter) write(header *sam.Header, ch chan *sam.Record) {
	w := createFile(c.outFile)
	defer w.Close()

	writer, err := bam.NewWriter(w, header, 0)
	raiseError(err)
	defer writer.Close()

	for r := range ch {
		writer.Write(r)
	}
}

type ProperPairFilter struct {
}

func (p *ProperPairFilter) Filter(r *sam.Record) bool {
	properPair := r.Flags&sam.ProperPair == sam.ProperPair
	unmapped := r.Flags&sam.Unmapped == sam.Unmapped
	mateUnmapped := r.Flags&sam.MateUnmapped == sam.MateUnmapped
	secondary := r.Flags&sam.Secondary == sam.Secondary
	supplementary := r.Flags&sam.Supplementary == sam.Supplementary
	return properPair && (!unmapped) && (!secondary) && (!supplementary) && (!mateUnmapped)
}

type ProperMapFilter struct {
}

func (p *ProperMapFilter) Filter(r *sam.Record) bool {
	var readLen int
	var matchLen int
	var containInsertation bool = false
	for _, co := range r.Cigar {
		switch co.Type() {
		case sam.CigarMatch, sam.CigarMismatch, sam.CigarEqual:
			readLen += co.Len()
			matchLen += co.Len()
		case sam.CigarInsertion:
			containInsertation = true
			readLen += co.Len()
		case sam.CigarSoftClipped:
			readLen += co.Len()
		}
	}

	if containInsertation {
		return false
	}

	if readLen-matchLen > readLen/10 {
		return false
	}

	return true
}

type DiversityFilter struct {
	RefMap map[string][]byte
	Cutoff float64
}

func (s *DiversityFilter) Filter(r *sam.Record) bool {
	acc := r.Ref.Name()
	genome, found := s.RefMap[acc]
	if !found {
		if *debug {
			log.Printf("%s not found\n", acc)
		}
		return false
	}

	start := r.Start()
	end := r.End()
	if start < 0 || end > len(genome) {
		if *debug {
			text, err := r.MarshalSAM(sam.FlagDecimal)
			raiseError(err)
			log.Printf("acc: %s, genome length %d, read starts at %d and ends at %d: %s\n", acc, len(genome), start, end, text)
		}
		return false
	}
	refSeq := genome[start:end]
	diff := 0
	read := map2Ref(r)
	for i := 0; i < len(read); i++ {
		if read[i] != refSeq[i] {
			diff++
		}
	}

	if float64(diff)/float64(len(read)) > s.Cutoff {
		return false
	}

	return true
}

type MapQFilter struct {
	Cutoff int
}

func (m *MapQFilter) Filter(r *sam.Record) bool {
	if int(r.MapQ) >= m.Cutoff {
		return true
	}

	return false
}

func readRefSeq(fnaFile string) map[string][]byte {
	f := openFile(fnaFile)
	fastaReader := seq.NewFastaReader(f)
	fastaReader.DeflineParser = func(s string) string { return strings.Split(strings.TrimSpace(s), " ")[0] }
	m := make(map[string][]byte)
	sequences, err := fastaReader.ReadAll()
	raiseError(err)
	for _, s := range sequences {
		m[s.Id] = bytes.ToUpper(s.Seq)
	}
	return m
}

// Obtain the sequence of a read mapping to the reference genome.
// Return the mapped sequence.
func map2Ref(r *sam.Record) []byte {
	s := []byte{}
	p := 0                                // position in the read sequence.
	read := bytes.ToUpper(r.Seq.Expand()) // read sequence.
	for _, c := range r.Cigar {
		switch c.Type() {
		case sam.CigarMatch, sam.CigarMismatch, sam.CigarEqual:
			s = append(s, read[p:p+c.Len()]...)
			p += c.Len()
		case sam.CigarInsertion, sam.CigarSoftClipped, sam.CigarHardClipped:
			p += c.Len()
		case sam.CigarDeletion, sam.CigarSkipped:
			s = append(s, bytes.Repeat([]byte{'*'}, c.Len())...)
		}
	}

	return s
}
