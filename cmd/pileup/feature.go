package main

type Feature struct {
	TaxID       string
	Species     string
	Genome      string
	Annotation  string
	Type        string
	PatricID    string
	LocusTag    string
	AltLocusTag string
	Uniprotkb   string
	Start       int
	End         int
	Strand      string
	Length      int
	Gene        string
	Product     string
	FigfamID    string
	PlfamID     string
	PgfamID     string
	Go          string
	Ec          string
	Pathway     string
}

type Features []*Feature

func (s Features) Len() int      { return len(s) }
func (s Features) Swap(i, j int) { s[i], s[j] = s[j], s[i] }

type ByStart struct{ Features }
type ByEnd struct{ Features }

func (s ByStart) Less(i, j int) bool { return s.Features[i].Start < s.Features[j].Start }
func (s ByEnd) Less(i, j int) bool   { return s.Features[i].End < s.Features[j].End }
