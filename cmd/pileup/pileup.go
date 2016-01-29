package main

import (
	"bufio"
	"bytes"
	"github.com/alecthomas/kingpin"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/mingzhi/biogo/seq"
	. "github.com/mingzhi/pileup"
	"io"
	"log"
	"os"
)

var (
	debug     = kingpin.Flag("debug", "Enable debug mode.").Bool()
	minMQ     = kingpin.Flag("min-mq", "minimum mapping quality.").Default("0").Short('q').Int()
	minBQ     = kingpin.Flag("min-bq", "minimum base quality").Default("13").Short('Q').Int()
	fastaFile = kingpin.Flag("fasta-file", "reference genome in FASTA format").Short('f').String()
	outFile   = kingpin.Flag("output", "output file").Short('o').String()
	bamFile   = kingpin.Arg("bam-file", "bam or sam file").Required().String()
)

func init() {
	kingpin.Version("0.3")
	kingpin.Parse()
}

func main() {
	_, readChan := readBamFile(*bamFile)
	filteredReadChan := filterReads(readChan)
	mappedReadChan := mapReads(filteredReadChan)
	genome := []byte{}
	if *fastaFile != "" {
		genome = readGenome(*fastaFile)
	}
	snpChan := pileupReads(mappedReadChan, genome)
	var f *os.File
	if *outFile != "" {
		f = createFile(*outFile)
	} else {
		f = os.Stdout
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	writeSNP(snpChan, w)
	w.Flush()
}

type MappedRead struct {
	Ref  string
	ID   string // ID
	Pos  int    // 0-based leftmost mapping POSition of the first matching base.
	Seq  []byte // base sequence.
	Qual []byte // qualities for each base.
	MapQ byte
}

func pileupReads(mappedReadChan chan *MappedRead, genome []byte) chan *SNP {
	c := make(chan *SNP)
	go func() {
		defer close(c)
		var currentPos int
		var currentRef string
		var currentReadID int
		snpBuf := []*SNP{}
		for mr := range mappedReadChan {
			if mr.Pos < currentPos {
				log.Fatalf("Current Position: %d, Read Position: %d\n", currentPos, mr.Pos)
			}

			if currentRef != mr.Ref {
				for _, s := range snpBuf {
					c <- s
				}
				snpBuf = []*SNP{}
				currentRef = mr.Ref
				currentPos = mr.Pos
			}

			for currentPos < mr.Pos {
				if len(snpBuf) > 0 {
					s := snpBuf[0]
					if s != nil {
						c <- s
					}
					snpBuf = snpBuf[1:]
				}
				currentPos++
			}

			currentReadID++

			for i := range mr.Seq {
				pos := mr.Pos + i
				a := Allele{
					Base:   mr.Seq[i],
					Qual:   mr.Qual[i],
					ReadID: currentReadID,
					MapQ:   mr.MapQ,
				}

				index := pos - currentPos
				if len(snpBuf) > index {
					if int(a.Qual) > *minBQ {
						snpBuf[index].Alleles = append(snpBuf[index].Alleles, a)
					}
				} else {
					for len(snpBuf) <= index {
						s := SNP{
							Ref: mr.Ref,
							Pos: pos,
						}
						if len(genome) > s.Pos {
							s.Base = genome[s.Pos]
						} else {
							s.Base = 'N'
						}
						snpBuf = append(snpBuf, &s)
					}
				}
			}
		}

		for i := range snpBuf {
			c <- snpBuf[i]
		}
	}()

	return c
}

// mapReads maps read to the reference genome and obtain the only mapped part.
func mapReads(readChan chan *sam.Record) chan *MappedRead {
	c := make(chan *MappedRead)
	go func() {
		defer close(c)
		for r := range readChan {
			s, q := mapRead2Ref(r)
			mr := MappedRead{
				Ref:  r.Ref.Name(),
				ID:   r.Name,
				Pos:  r.Pos,
				Seq:  s,
				Qual: q,
				MapQ: r.MapQ,
			}
			c <- &mr
		}
	}()

	return c
}

// mapRead2Ref Obtains a read mapping to the reference genome.
func mapRead2Ref(r *sam.Record) (s []byte, q []byte) {
	p := 0                 // position in the read sequence.
	read := r.Seq.Expand() // read sequence.
	qual := r.Qual
	for _, c := range r.Cigar {
		switch c.Type() {
		case sam.CigarMatch, sam.CigarMismatch, sam.CigarEqual:
			s = append(s, read[p:p+c.Len()]...)
			q = append(q, qual[p:p+c.Len()]...)
			p += c.Len()
		case sam.CigarInsertion, sam.CigarSoftClipped, sam.CigarHardClipped:
			p += c.Len()
		case sam.CigarDeletion, sam.CigarSkipped:
			for i := 0; i < c.Len(); i++ {
				s = append(s, '*')
				q = append(q, 0)
			}
		}
	}

	s = bytes.ToUpper(s)

	return
}

// filterReads filter low quality reads.
func filterReads(readChan chan *sam.Record) (c chan *sam.Record) {
	c = make(chan *sam.Record)
	go func() {
		defer close(c)
		var goodReads, badReads int
		for r := range readChan {
			// checking mapping quality
			mapQ := int(r.MapQ)
			if mapQ <= *minMQ || mapQ == 255 {
				badReads++
			} else {
				goodReads++
				c <- r
			}

		}

		if *debug {
			log.Printf("Total used reads: %d\n", goodReads)
			log.Printf("Total discard reads: %d\n", badReads)
		}

	}()

	return c
}

// ReadBamFile reads bam file, and return the header and a channel of sam records.
func readBamFile(fileName string) (h *sam.Header, c chan *sam.Record) {
	// Initialize the channel of sam records.
	c = make(chan *sam.Record)

	// Create a new go routine to read the records.
	go func() {
		// Close the record channel when finished.
		defer close(c)

		// Open file stream, and close it when finished.
		f, err := os.Open(fileName)
		if err != nil {
			panic(err)
		}
		defer f.Close()

		type SamReader interface {
			Header() *sam.Header
			Read() (*sam.Record, error)
		}

		var reader SamReader
		if fileName[len(fileName)-3:] == "bam" {
			bamReader, err := bam.NewReader(f, 0)
			if err != nil {
				panic(err)
			}
			defer bamReader.Close()
			reader = bamReader
		} else {
			reader, err = sam.NewReader(f)
			if err != nil {
				panic(err)
			}
		}

		// Read and assign header.
		h = reader.Header()

		// Read sam records and send them to the channel,
		// until it hit an error, which raises a panic
		// if it is not a IO EOF.
		for {
			rec, err := reader.Read()
			if err != nil {
				if err != io.EOF {
					panic(err)
				}
				break
			}
			c <- rec
		}
		if *debug {
			log.Println("Finished reading bam file!")
		}
	}()

	return
}

func writeSNP(snpChan chan *SNP, w io.Writer) {
	buf := bufio.NewWriter(w)
	for s := range snpChan {
		buf.WriteString(s.String() + "\n")
	}
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

func createFile(filename string) *os.File {
	w, err := os.Create(filename)
	if err != nil {
		panic(err)
	}

	return w
}
