package main

import (
	"bufio"
	"github.com/mingzhi/pileup"
	"io"
	"os"
)

func readPileup(f *os.File) chan *pileup.SNP {
	c := make(chan *pileup.SNP)
	go func() {
		defer close(c)
		defer f.Close()
		rd := bufio.NewReader(f)
		for {
			line, err := rd.ReadString('\n')
			if err != nil {
				if err != io.EOF {
					panic(err)
				}
				break
			}
			s := pileup.Parse(line)
			c <- s
		}
	}()

	return c
}
