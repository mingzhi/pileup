package main

import (
	"encoding/json"
	"github.com/mingzhi/pileup"
	"io"
	"log"
	"os"
)

func readPileup(f *os.File) chan *pileup.SNP {
	c := make(chan *pileup.SNP)
	go func() {
		defer close(c)
		decoder := json.NewDecoder(f)
		for decoder.More() {
			var s *pileup.SNP
			err := decoder.Decode(&s)
			if err != nil {
				log.Fatalln(err)
			}
			c <- s
		}
	}()

	return c
}

func readMPileup(f *os.File) chan *pileup.SNP {
	c := make(chan *pileup.SNP)
	go func() {
		defer close(c)
		reader := pileup.NewReader(f)
		for {
			s, err := reader.Read()
			if err != nil {
				if err != io.EOF {
					log.Panic(err)
				}
			}
			c <- s
		}
	}()
	return c
}
