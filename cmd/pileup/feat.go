package main

import (
	"bufio"
	"github.com/boltdb/bolt"
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
}

type dbfunc func(tx *bolt.Tx) error

func (f *cmdFeat) run() {
	db, err := bolt.Open(f.out, 0600, nil)
	if err != nil {
		log.Fatal(err)
	}
	defer db.Close()

	if err := db.Update(createBucket("feature")); err != nil {
		log.Panicln(err)
	}

	if err := db.Update(createBucket("genome")); err != nil {
		log.Panicln(err)
	}

	featureFileList := f.walk()
	for _, featureFile := range featureFileList {
		features := readFeatures(featureFile)
		loadFeatures(db, features)
	}
}

func createBucket(bucketName string) dbfunc {
	f := func(tx *bolt.Tx) error {
		_, err := tx.CreateBucket([]byte(bucketName))
		if err != nil {
			if err == bolt.ErrBucketExists {
				tx.DeleteBucket([]byte(bucketName))
				tx.CreateBucket([]byte(bucketName))
			} else {
				log.Panicln(err)
			}
		}
		return nil
	}
	return f
}

func loadFeatures(db *bolt.DB, features []Feature) {
	err := db.Update(func(tx *bolt.Tx) error {
		genomeGeneMap := make(map[string][]string)
		for _, f := range features {
			genomeGeneMap[f.Genome] = append(genomeGeneMap[f.Genome], f.PatricID)
		}

		b := tx.Bucket([]byte("genome"))
		for genome, ids := range genomeGeneMap {
			key := []byte(genome)
			value, err := msgpack.Marshal(ids)
			if err != nil {
				return err
			}
			if err := b.Put(key, value); err != nil {
				return err
			}
		}

		b = tx.Bucket([]byte("feature"))
		for _, f := range features {
			if f.PatricID == "" {
				continue
			}
			key := []byte(f.PatricID)
			value, err := msgpack.Marshal(f)
			if err != nil {
				return err
			}
			if err := b.Put(key, value); err != nil {
				return err
			}
		}

		return nil
	})

	if err != nil {
		log.Panicln(err)
	}
}

// walk returns a list of feature files.
func (f *cmdFeat) walk() []string {
	featureFileList := []string{}
	fileInfoList := readDir(f.dir)
	for _, fi := range fileInfoList {
		if fi.IsDir() {
			dirPath := filepath.Join(f.dir, fi.Name())
			filelist := readDir(dirPath)
			for _, f := range filelist {
				if strings.Contains(f.Name(), ".PATRIC.features.tab") {
					filePath := filepath.Join(dirPath, f.Name())
					featureFileList = append(featureFileList, filePath)
				}
			}
		}
	}

	return featureFileList
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
