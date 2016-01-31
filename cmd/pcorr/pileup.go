package main

import (
	"bufio"
	"encoding/json"
	"github.com/biogo/hts/sam"
	. "github.com/mingzhi/pileup"
	"io"
	"log"
	"os"
	"sort"
)

type cmdPileup struct {
	debug                       bool
	minBQ, minMQ                int
	bamFile, fastaFile, outFile string
}

func (cmd *cmdPileup) Run() {
	_, readChan := readBamFile(cmd.bamFile)
	filteredReadChan := cmd.filterReads(readChan)
	mappedReadChan := cmd.mapReads(filteredReadChan)
	genome := []byte{}
	if cmd.fastaFile != "" {
		genome = readGenome(cmd.fastaFile)
	}
	snpChan := cmd.pileupReads(mappedReadChan, genome)
	var f *os.File
	if cmd.outFile != "" {
		f = createFile(cmd.outFile)
	} else {
		f = os.Stdout
	}
	defer f.Close()

	cmd.writeSNP(snpChan, f)
}

type MappedRead struct {
	Ref  string
	ID   string // ID
	Pos  int    // 0-based leftmost mapping POSition of the first matching base.
	Seq  []byte // base sequence.
	Qual []byte // qualities for each base.
	MapQ byte
}

func (cmd *cmdPileup) pileupReads(mappedReadChan chan MappedRead, genome []byte) chan SNP {
	c := make(chan SNP)
	go func() {
		defer close(c)
		var currentReadID int

		buffer := make(map[int]*SNP)
		for mr := range mappedReadChan {
			currentReadID++
			for i := 0; i < len(mr.Seq); i++ {
				pos := mr.Pos + i
				if int(mr.Qual[i]) > cmd.minBQ {
					a := Allele{
						Base:   mr.Seq[i],
						Qual:   mr.Qual[i],
						ReadID: currentReadID,
						MapQ:   mr.MapQ,
					}

					if _, found := buffer[pos]; !found {
						s := SNP{
							Ref: mr.Ref,
							Pos: pos,
						}
						if len(genome) > s.Pos {
							s.Base = genome[s.Pos]
						} else {
							s.Base = 'N'
						}
						buffer[pos] = &s
					}
					buffer[pos].Alleles = append(buffer[pos].Alleles, a)
				}
			}

			if len(buffer) > 200 {
				positions := []int{}
				for pos := range buffer {
					positions = append(positions, pos)
				}
				sort.Ints(positions)
				for _, pos := range positions {
					snp := buffer[pos]
					if pos < mr.Pos {
						if snp.Base != 'N' {
							bases := []byte{}
							for _, a := range snp.Alleles {
								bases = append(bases, a.Base)
							}
							c <- *snp
						}
						delete(buffer, pos)
					}
				}
			}
		}

		for _, snp := range buffer {
			c <- *snp
		}
	}()

	return c
}

// mapReads maps read to the reference genome and obtain the only mapped part.
func (cmd *cmdPileup) mapReads(readChan chan sam.Record) chan MappedRead {
	c := make(chan MappedRead)
	go func() {
		defer close(c)
		for r := range readChan {
			s, q := cmd.mapRead2Ref(r)
			if len(s) > 0 {
				mr := MappedRead{
					Ref:  r.Ref.Name(),
					ID:   r.Name,
					Pos:  r.Pos,
					Seq:  s,
					Qual: q,
					MapQ: r.MapQ,
				}
				c <- mr
			}
		}
	}()

	return c
}

// mapRead2Ref Obtains a read mapping to the reference genome.
func (cmd *cmdPileup) mapRead2Ref(r sam.Record) (s []byte, q []byte) {
	read := r.Seq.Expand() // read sequence.
	qual := r.Qual
	for _, c := range r.Cigar {
		if c.Type() != sam.CigarMatch && c.Type() != sam.CigarMismatch && c.Type() != sam.CigarEqual {
			if cmd.debug {
				log.Println(r.Cigar)
			}

			return
		}
	}
	s = read
	q = qual
	return
}

// filterReads filter low quality reads.
func (cmd *cmdPileup) filterReads(readChan chan sam.Record) (c chan sam.Record) {
	c = make(chan sam.Record)
	go func() {
		defer close(c)
		var goodReads, badReads int
		for r := range readChan {
			// checking mapping quality
			mapQ := int(r.MapQ)
			if mapQ <= cmd.minMQ || mapQ == 255 {
				badReads++
			} else {
				goodReads++
				c <- r
			}

		}

		if cmd.debug {
			log.Printf("Total used reads: %d\n", goodReads)
			log.Printf("Total discard reads: %d\n", badReads)
		}

	}()

	return c
}

func (cmd *cmdPileup) writeSNP(snpChan chan SNP, w io.Writer) {
	bw := bufio.NewWriter(w)
	defer bw.Flush()
	encoder := json.NewEncoder(bw)
	for s := range snpChan {
		s.Num = len(s.Alleles)
		if s.Num > 0 {
			if err := encoder.Encode(s); err != nil {
				log.Fatalln(err)
			}
		}
	}
}
