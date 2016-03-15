package main

import (
	"bytes"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/bmatsuo/lmdb-go/lmdb"
	"github.com/mingzhi/biogo/seq"
	"io"
	"io/ioutil"
	"log"
	"os"
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

	header := c.readBamHeader()

	ch := c.readBamFile()
	filteredChan1 := c.filter(ch)

	tempDir, err := ioutil.TempDir("temp", "filter")
	raiseError(err)
	defer os.RemoveAll(tempDir)

	env, err := lmdb.NewEnv()
	raiseError(err)
	err = env.SetMaxDBs(10)
	raiseError(err)
	err = env.SetMapSize(1 << 31)
	raiseError(err)
	env.Open(tempDir, 0, 0644)
	defer env.Close()

	diverFiler := DiversityFilter{
		RefFile: c.fnaFile,
		Cutoff:  c.maxDistance,
		MDB:     env,
	}
	diverFiler.Init()
	filteredChan2 := diverFiler.FilterAll(filteredChan1)

	c.write(header, filteredChan2)
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
	RefFile string
	Cutoff  float64
	MDB     *lmdb.Env
}

func (d *DiversityFilter) Init() {
	createDBI(d.MDB, "read")
	createDBI(d.MDB, "genome")
	d.loadGenome()
}

func (d *DiversityFilter) loadGenome() {
	f := openFile(d.RefFile)
	defer f.Close()

	fastaReader := seq.NewFastaReader(f)
	fastaReader.DeflineParser = func(s string) string { return strings.Split(strings.TrimSpace(s), " ")[0] }
	buf := make([]*seq.Sequence, 100)
	k := 0
	for {
		s, err := fastaReader.Read()
		if err != nil {
			if err != io.EOF {
				raiseError(err)
			}
			break
		}

		if k >= len(buf) {
			d.loadGenomeBuf(buf)
			k = 0
		}
		buf[k] = s
		k++
	}

	d.loadGenomeBuf(buf)
	log.Println("finish loading genome")
}

func (d *DiversityFilter) loadGenomeBuf(genomes []*seq.Sequence) {
	fn := func(txn *lmdb.Txn) error {
		dbi, err := txn.OpenDBI("genome", 0)
		if err != nil {
			return err
		}
		for _, s := range genomes {
			key := []byte(s.Id)
			val := bytes.ToUpper(s.Seq)
			err = txn.Put(dbi, key, val, 0)
			if err != nil {
				return err
			}
		}
		return nil
	}

	err := d.MDB.Update(fn)
	raiseError(err)
}

func (d *DiversityFilter) filter(buf []*sam.Record, acc string, genome []byte) (out []*sam.Record, acc1 string, genome1 []byte) {
	fn := func(txn *lmdb.Txn) error {
		dbi, err := txn.OpenDBI("read", 0)
		if err != nil {
			return err
		}
		genomeDbi, err := txn.OpenDBI("genome", 0)
		if err != nil {
			return err
		}
		for _, r := range buf {
			key := []byte(r.Name)
			val, err := txn.Get(dbi, key)
			if err != nil {
				if lmdb.IsNotFound(err) {
					val, err = r.MarshalText()
					if err != nil {
						return err
					}
					err = txn.Put(dbi, key, val, 0)
					if err != nil {
						return err
					}
				} else {
					return err
				}
			} else {
				var mate *sam.Record = &sam.Record{}
				err := mate.UnmarshalText(val)
				raiseError(err)
				if r.Ref.Name() == mate.Ref.Name() {
					if acc != r.Ref.Name() {
						genome, err = d.findGenome(r, txn, genomeDbi)
						raiseError(err)
						acc = r.Ref.Name()
						if *debug {
							log.Println(acc)
						}
					}

					diff1, len1 := d.Diff(r, genome)
					diff2, len2 := d.Diff(mate, genome)
					if len1 > 0 && len2 > 0 && float64(diff1+diff2)/float64(len1+len2) <= d.Cutoff {
						out = append(out, r)
						out = append(out, mate)
					} else {
						if *debug {
							log.Printf("%d, %d, %d, %d\n", diff1, diff2, len1, len2)
						}
					}
				}

				txn.Del(dbi, key, val)
			}

		}
		return nil
	}
	err := d.MDB.Update(fn)
	raiseError(err)

	genome1 = genome
	acc1 = acc
	return
}

func (d *DiversityFilter) findGenome(r *sam.Record, txn *lmdb.Txn, dbi lmdb.DBI) (val []byte, err error) {
	key := []byte(r.Ref.Name())
	val, err = txn.Get(dbi, key)
	return
}

func (d *DiversityFilter) FilterAll(ch chan *sam.Record) chan *sam.Record {
	filteredChan := make(chan *sam.Record)

	go func() {
		defer close(filteredChan)
		buf := make([]*sam.Record, 10000)
		var out []*sam.Record
		var acc string
		var genome []byte
		k := 0
		for r := range ch {
			if k >= len(buf) {
				out, acc, genome = d.filter(buf, acc, genome)
				for _, r1 := range out {
					filteredChan <- r1
				}
				k = 0
			}
			buf[k] = r
			k++
		}

		out, _, _ = d.filter(buf, acc, genome)
		for _, r1 := range out {
			filteredChan <- r1
		}
	}()

	return filteredChan
}

func (d *DiversityFilter) Diff(r *sam.Record, genome []byte) (diff, length int) {
	start := r.Start()
	end := r.End()
	if start < 0 || end > len(genome) {
		if *debug {
			text, err := r.MarshalSAM(sam.FlagDecimal)
			raiseError(err)
			log.Printf("acc: %s, genome length %d, read starts at %d and ends at %d: %s\n", r.Ref.Name(), len(genome), start, end, text)
		}
		length = 0
		return
	}
	refSeq := genome[start:end]
	diff = 0
	read := map2Ref(r)
	length = len(read)
	for i := 0; i < length; i++ {
		if read[i] != refSeq[i] {
			diff++
		}
	}
	return
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
