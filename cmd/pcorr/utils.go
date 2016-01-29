package main

import (
	"os"
)

func openFile(filename string) *os.File {
	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}

	return f
}

func createFile(filename string) *os.File {
	w, err := os.Create(filename)
	if err != nil {
		panic(err)
	}

	return w
}
