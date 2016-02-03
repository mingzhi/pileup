package main

import (
	"encoding/json"
	"io"
	"log"
	"os"
	"strings"

	"github.com/mingzhi/pileup"
)

func readPileup(f *os.File, regionStart, regionEnd int) <-chan *pileup.SNP {
	c := make(chan *pileup.SNP)
	go func() {
		defer close(c)
		// check format.
		fields := strings.Split(f.Name(), ".")
		format := fields[len(fields)-1]
		var snpChan <-chan *pileup.SNP
		done := make(chan struct{})
		defer close(done)
		switch format {
		case "pileup":
			snpChan = readPileup1(f, done)
			break
		case "mpileup":
			snpChan = readPileup2(f, done)
			break
		default:
			log.Fatalf("Can not recognize the pileup format: %s.\nFile should be ended with .pileup or .mpileup\n", format)
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

func readPileup1(f *os.File, done <-chan struct{}) <-chan *pileup.SNP {
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

func readPileup2(f *os.File, done <-chan struct{}) <-chan *pileup.SNP {
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
