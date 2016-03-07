package main

import (
	"github.com/bmatsuo/lmdb-go/lmdb"
	"github.com/mingzhi/gomath/stat/desc/meanvar"
	"log"
)

func newEnv() *lmdb.Env {
	env, err := lmdb.NewEnv()
	if err != nil {
		log.Fatalln(err)
	}

	err = env.SetMaxDBs(10)
	if err != nil {
		log.Fatalln(err)
	}

	err = env.SetMapSize(1 << 32)
	if err != nil {
		log.Fatalln(err)
	}

	return env
}

func createReadOnlyEnv(path string) *lmdb.Env {
	env := newEnv()
	err := env.Open(path, lmdb.Readonly, 0644)
	if err != nil {
		log.Fatalln(err)
	}
	return env
}

func createEnv(path string) *lmdb.Env {
	env := newEnv()
	err := env.Open(path, 0, 0644)
	if err != nil {
		log.Panicln(err)
	}

	return env
}

func createDBI(env *lmdb.Env, name string) {
	fn := func(txn *lmdb.Txn) error {
		var dbi lmdb.DBI
		var err error
		var del bool = false

		if dbi, err = txn.CreateDBI(name); err != nil {
			return err
		}

		if err = txn.Drop(dbi, del); err != nil {
			return err
		}

		return nil
	}
	err := env.Update(fn)
	if err != nil {
		log.Fatalln(err)
	}
}

func newMeanVars(size int) []*meanvar.MeanVar {
	mvs := make([]*meanvar.MeanVar, size)
	for i := range mvs {
		mvs[i] = meanvar.New()
	}
	return mvs
}
