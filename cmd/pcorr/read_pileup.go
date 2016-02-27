package main

import (
	"encoding/json"
	"io"
	"log"
	"os"

	"github.com/mingzhi/pileup"
)

func readPileup(f *os.File, regionStart, regionEnd int, fileFormate string) <-chan *pileup.SNP {
	c := make(chan *pileup.SNP)
	go func() {
		defer close(c)
		var snpChan <-chan *pileup.SNP
		done := make(chan struct{})
		defer close(done)
		switch fileFormate {
		case "json":
			snpChan = readPileupJSON(f, done)
			break
		case "tab":
			snpChan = readPileupTab(f, done)
			break
		default:
			log.Fatalf("Can not recognize the pileup format: %s\n", fileFormate)
		}

		for s := range snpChan {
			if regionEnd > 0 {
				if s.Pos >= regionEnd {
					done <- struct{}{}
				} else if s.Pos >= regionStart {
					c <- s
				}
			} else {
				if s.Pos >= regionStart {
					c <- s
				}
			}
		}
	}()
	return c
}

func readPileupJSON(f *os.File, done <-chan struct{}) <-chan *pileup.SNP {
	c := make(chan *pileup.SNP)
	go func() {
		defer close(c)
		decoder := json.NewDecoder(f)
		for decoder.More() {
			var s *pileup.SNP
			err := decoder.Decode(&s)
			if err != nil {
				if *debug {
					log.Panic(err)
				} else {
					log.Fatalln(err)
				}
			}

			select {
			case c <- s:
			case <-done:
				return
			}
		}
	}()

	return c
}

func readPileupTab(f *os.File, done <-chan struct{}) <-chan *pileup.SNP {
	c := make(chan *pileup.SNP)
	go func() {
		defer close(c)
		reader := pileup.NewReader(f)
		for {
			s, err := reader.Read()
			if err != nil {
				if err != io.EOF {
					if *debug {
						log.Panic(err)
					} else {
						log.Fatalln(err)
					}

				}
				break
			}

			select {
			case c <- s:
			case <-done:
				return
			}
		}
	}()
	return c
}
