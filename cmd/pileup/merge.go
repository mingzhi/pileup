package main

import (
	"bufio"
	"github.com/bmatsuo/lmdb-go/lmdb"
	"io"
	"log"
	"os"
	"strings"
)

type cmdMerge struct {
	sampleFile string
	dbiName    string
	dbOut      string
}

func (c *cmdMerge) run() {
	env := createEnv(c.dbOut)
	defer env.Close()

	fn := func(txn *lmdb.Txn) error {
		dbi, err := txn.OpenDBI(c.dbiName, 0)
		if err != nil {
			return err
		}
		kvChan := c.readAllCr()
		for kv := range kvChan {
			key, val := kv.Key, kv.Value
			err := txn.Put(dbi, key, val, 0)
			if err != nil {
				return err
			}
		}

		return nil
	}

	err := env.Update(fn)
	if err != nil {
		if *debug {
			log.Panicln(err)
		} else {
			log.Fatalln(err)
		}
	}
}

func (c *cmdMerge) readSamples() []string {
	f, err := os.Open(c.sampleFile)
	if err != nil {
		log.Fatalln(err)
	}
	defer f.Close()

	list := []string{}
	r := bufio.NewReader(f)
	for {
		line, err := r.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			}
			if *debug {
				log.Println(err)
			} else {
				log.Fatalln(err)
			}
		}
		l := strings.TrimSpace(line)
		list = append(list, l)
	}

	return list
}

func (c *cmdMerge) readAllCr() chan KeyValue {
	kvChan := make(chan KeyValue)
	go func() {
		defer close(kvChan)
		samples := c.readSamples()
		for _, s := range samples {
			path := s + "_mdb"
			env := createReadOnlyEnv(path)
			defer env.Close()
			ch := getAllCr(env, c.dbiName)
			for kv := range ch {
				k := []byte(s + "_" + string(kv.Key))
				kv.Key = k
				kvChan <- kv
			}

			if *debug {
				log.Println(s)
			}
		}
	}()

	return kvChan
}

func getAllCr(env *lmdb.Env, dbname string) chan KeyValue {
	ch := make(chan KeyValue)
	go func() {
		defer close(ch)
		fn := func(txn *lmdb.Txn) error {
			dbi, err := txn.OpenDBI(dbname, 0)
			if err != nil {
				return err
			}

			cur, err := txn.OpenCursor(dbi)
			if err != nil {
				return err
			}

			for {
				k, v, err := cur.Get(nil, nil, lmdb.Next)
				if lmdb.IsNotFound(err) {
					break
				}

				if err != nil {
					return err
				}

				ch <- KeyValue{Key: k, Value: v}
			}
			return nil
		}

		err := env.View(fn)
		if err != nil {
			if *debug {
				log.Panicln(err)
			} else {
				log.Fatalln(err)
			}
		}
	}()
	return ch
}
