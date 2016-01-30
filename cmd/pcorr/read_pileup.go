package main

import (
	"encoding/json"
	"github.com/mingzhi/pileup"
	"log"
	"os"
)

func readPileup(f *os.File) chan pileup.SNP {
	c := make(chan pileup.SNP)
	go func() {
		defer close(c)
		defer f.Close()
		decoder := json.NewDecoder(f)
		for decoder.More() {
			var s pileup.SNP
			err := decoder.Decode(&s)
			if err != nil {
				log.Fatalln(err)
			}
			c <- s
		}
	}()

	return c
}
