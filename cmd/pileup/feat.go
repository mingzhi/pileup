package main

import (
	"bufio"
	"github.com/bmatsuo/lmdb-go/lmdb"
	"github.com/mingzhi/biogo/seq"
	"gopkg.in/vmihailenco/msgpack.v2"
	"io"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

type cmdFeat struct {
	dir string
	out string

	env *lmdb.Env
}

func (c *cmdFeat) run() {
	// create an environment and make sure it is eventually closed.
	c.env = createEnv(c.out)
	defer c.env.Close()

	createDBI(c.env, "feature")
	createDBI(c.env, "genome")

	featureFileList := c.walk(".PATRIC.features.tab")
	for _, featureFile := range featureFileList {
		features := readFeatures(featureFile)
		c.loadFeatures(features)
	}
}

func (c *cmdFeat) loadFeatures(features []Feature) {
	err := c.env.Update(func(txn *lmdb.Txn) error {
		genomeGeneMap := make(map[string][]string)
		for _, f := range features {
			genomeGeneMap[f.Genome] = append(genomeGeneMap[f.Genome], f.PatricID)
		}

		var dbi lmdb.DBI
		var err error
		dbi, err = txn.OpenDBI("genome", 0)
		if err != nil {
			return err
		}
		for genome, ids := range genomeGeneMap {
			key := []byte(genome)
			value, err := msgpack.Marshal(ids)
			if err != nil {
				return err
			}
			if err := txn.Put(dbi, key, value, 0); err != nil {
				return err
			}
		}

		dbi, err = txn.OpenDBI("feature", 0)
		if err != nil {
			return err
		}
		for _, f := range features {
			if f.PatricID == "" {
				continue
			}
			key := []byte(f.PatricID)
			value, err := msgpack.Marshal(f)
			if err != nil {
				return err
			}
			if err := txn.Put(dbi, key, value, 0); err != nil {
				return err
			}
		}

		return nil
	})

	if err != nil {
		log.Fatalln(err)
	}
}

// walk returns a list of feature files.
func (f *cmdFeat) walk(pattern string) []string {
	featureFileList := []string{}
	fileInfoList := readDir(f.dir)
	for _, fi := range fileInfoList {
		if fi.IsDir() {
			dirPath := filepath.Join(f.dir, fi.Name())
			filelist := readDir(dirPath)
			for _, f := range filelist {
				if strings.Contains(f.Name(), pattern) {
					filePath := filepath.Join(dirPath, f.Name())
					featureFileList = append(featureFileList, filePath)
				}
			}
		}
	}

	return featureFileList
}

func readFna(filename string, getAcc func(string) string) (accessions []string, sequences [][]byte) {
	f, err := os.Open(filename)
	if err != nil {
		log.Panicln(err)
	}
	defer f.Close()
	fastaReader := seq.NewFastaReader(f)
	ss, err := fastaReader.ReadAll()
	if err != nil {
		log.Panicln(err)
	}
	for _, a := range ss {
		id := getAcc(a.Id)
		accessions = append(accessions, id)
		sequences = append(sequences, a.Seq)
	}
	return
}

func readFeatures(filename string) []Feature {
	f, err := os.Open(filename)
	if err != nil {
		log.Panicln(err)
	}
	defer f.Close()

	features := []Feature{}
	rd := bufio.NewReader(f)
	rd.ReadString('\n')
	for {
		line, err := rd.ReadString('\n')
		if err != nil {
			if err != io.EOF {
				log.Panicln(err)
			}
			break
		}
		features = append(features, parseFeatureLine(line))
	}
	return features
}

func parseFeatureLine(line string) Feature {
	l := line[0 : len(line)-1]
	terms := strings.Split(l, "\t")
	f := Feature{}
	f.TaxID = strings.Split(terms[0], ".")[0]
	f.Species = terms[1]
	f.Genome = terms[2]
	f.Annotation = terms[3]
	f.Type = terms[4]
	f.PatricID = terms[5]
	f.LocusTag = terms[6]
	f.AltLocusTag = terms[7]
	f.Uniprotkb = terms[8]
	f.Start = atoi(terms[9])
	f.End = atoi(terms[10])
	f.Strand = terms[11]
	f.Length = atoi(terms[12])
	if len(terms) > 13 {
		f.Gene = terms[13]
		f.Product = terms[14]
		f.FigfamID = terms[15]
		f.PlfamID = terms[16]
		f.PgfamID = terms[17]
		f.Go = terms[18]
		f.Ec = terms[19]
		f.Pathway = terms[20]
	}

	return f
}

func isDirExist(dir string) bool {
	if _, err := os.Stat(dir); err != nil {
		return false
	} else {
		return true
	}
}

func readDir(dir string) []os.FileInfo {
	if !isDirExist(dir) {
		log.Panicf("%s does not exist!\n", dir)
	}

	fileInfoList, err := ioutil.ReadDir(dir)
	if err != nil {
		log.Panicln(err)
	}
	return fileInfoList
}

func atoi(s string) int {
	i, err := strconv.Atoi(s)
	if err != nil {
		log.Println(s)
		log.Panicln(err)
	}
	return i
}

func getID(id string) string {
	terms := strings.Split(id, "|")
	if strings.Contains(id, "fig|") {
		return strings.Join(terms[:2], "|")
	} else if strings.Contains(id, "accn|") {
		return terms[1]
	} else {
		if *debug {
			log.Println(id)
		}
		return id
	}

}
