package main

import (
	"github.com/bmatsuo/lmdb-go/lmdb"
	"github.com/mingzhi/biogo/pileup"
	"gopkg.in/vmihailenco/msgpack.v2"
)

type SNPArr struct {
	Key []byte
	Arr []pileup.SNP
}

func getAllSNPs(env *lmdb.Env) chan SNPArr {
	ch := make(chan SNPArr)
	fn := func(tx *lmdb.Txn) error {
		dbi, err := tx.OpenDBI("gene", 0)
		if err != nil {
			return err
		}

		cur, err := tx.OpenCursor(dbi)
		if err != nil {
			return err
		}
		defer cur.Close()

		for {
			k, v, err := cur.Get(nil, nil, lmdb.Next)
			if lmdb.IsNotFound(err) {
				return nil
			} else if err != nil {
				return err
			}

			arr := []pileup.SNP{}
			if err := msgpack.Unmarshal(v, &arr); err != nil {
				return err
			}

			ch <- SNPArr{Key: k, Arr: arr}
		}
	}

	go func() {
		defer close(ch)
		env.View(fn)
	}()

	return ch
}
