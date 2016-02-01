package pileup

import (
	"bufio"
	"io"
)

type Reader struct {
	r *bufio.Reader
}

func NewReader(r io.Reader) *Reader {
	d := Reader{}
	d.r = bufio.NewReader(r)
	return &d
}

func (d *Reader) Read() (s *SNP, err error) {
	var line string
	line, err = d.r.ReadString('\n')
	if err != nil {
		return nil, err
	}

	s, err = parse(line)
	return s, err
}
