package main

import (
	"github.com/bmatsuo/lmdb-go/lmdb"
	"github.com/mingzhi/gomath/stat/desc/meanvar"
	"log"
	"os"
)

type KeyValue struct {
	Key, Value []byte
}

func newEnv() *lmdb.Env {
	env, err := lmdb.NewEnv()
	raiseError(err)
	err = env.SetMaxDBs(10)
	raiseError(err)

	err = env.SetMapSize(1 << 42)
	raiseError(err)

	return env
}

func createReadOnlyEnv(path string) *lmdb.Env {
	env := newEnv()
	err := env.Open(path, lmdb.Readonly, 0644)
	raiseError(err)
	return env
}

func createNoLockEnv(path string) *lmdb.Env {
	env := newEnv()
	err := env.Open(path, lmdb.NoLock, 0644)
	raiseError(err)
	return env
}

func createEnv(path string) *lmdb.Env {
	env := newEnv()
	err := env.Open(path, 0, 0644)
	raiseError(err)

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
	raiseError(err)
}

func newMeanVars(size int) []*meanvar.MeanVar {
	mvs := make([]*meanvar.MeanVar, size)
	for i := range mvs {
		mvs[i] = meanvar.New()
	}
	return mvs
}

func openFile(filename string) *os.File {
	f, err := os.Open(filename)
	raiseError(err)

	return f
}

func createFile(filename string) *os.File {
	f, err := os.Create(filename)
	raiseError(err)

	return f
}

func raiseError(err error) {
	if err != nil {
		if *debug {
			log.Panic(err)
		} else {
			log.Fatalln(err)
		}
	}
}
