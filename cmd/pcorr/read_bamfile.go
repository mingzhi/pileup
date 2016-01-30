package main

import (
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"io"
	"log"
	"os"
)

// ReadBamFile reads bam file, and return the header and a channel of sam records.
func readBamFile(fileName string) (h *sam.Header, c chan sam.Record) {
	// Initialize the channel of sam records.
	c = make(chan sam.Record)

	// Create a new go routine to read the records.
	go func() {
		// Close the record channel when finished.
		defer close(c)

		// Open file stream, and close it when finished.
		f, err := os.Open(fileName)
		if err != nil {
			log.Fatalln(err)
		}
		defer f.Close()

		type SamReader interface {
			Header() *sam.Header
			Read() (*sam.Record, error)
		}

		var reader SamReader
		if fileName[len(fileName)-3:] == "bam" {
			bamReader, err := bam.NewReader(f, 0)
			if err != nil {
				log.Fatalln(err)
			}
			defer bamReader.Close()
			reader = bamReader
		} else {
			reader, err = sam.NewReader(f)
			if err != nil {
				log.Fatalln(err)
			}
		}

		// Read and assign header.
		h = reader.Header()

		// Read sam records and send them to the channel,
		// until it hit an error, which raises a panic
		// if it is not a IO EOF.
		for {
			rec, err := reader.Read()
			if err != nil {
				if err != io.EOF {
					log.Fatalln(err)
				}
				break
			}
			c <- *rec
		}
	}()

	return
}
