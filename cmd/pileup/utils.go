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

func newEnv(numDB int, sizeDB int64) *lmdb.Env {
	env, err := lmdb.NewEnv()
	raiseError(err)
	err = env.SetMaxDBs(numDB)
	raiseError(err)

	err = env.SetMapSize(sizeDB)
	raiseError(err)

	return env
}

func createReadOnlyEnv(path string, numDB int, sizeDB int64) *lmdb.Env {
	env := newEnv(numDB, sizeDB)
	err := env.Open(path, lmdb.Readonly, 0644)
	raiseError(err)
	return env
}

func createNoLockEnv(path string, numDB int, sizeDB int64) (env *lmdb.Env, err error) {
	env = newEnv(numDB, sizeDB)
	err = env.Open(path, lmdb.NoLock, 0644)
	return
}

func createEnv(path string, numDB int, sizeDB int64) (env *lmdb.Env, err error) {
	env = newEnv(numDB, sizeDB)
	if _, err := os.Stat(path); err != nil {
		if os.IsNotExist(err) {
			os.Mkdir(path, 0644)
		} else {
			raiseError(err)
		}
	}
	err = env.Open(path, 0, 0644)
	return
}

func createDBI(env *lmdb.Env, name string) error {
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
	return err
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
