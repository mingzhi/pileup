package main

import (
	"log"
	"os"
	"runtime"
	"runtime/pprof"

	"github.com/alecthomas/kingpin"
)

var (
	app     = kingpin.New("pcorr", "A command-line application for correlation calculation.")
	debug   = app.Flag("debug", "Enable debug mode.").Bool()
	ncpu    = app.Flag("ncpu", "number of CPUs for using").Default("1").Int()
	profile = app.Flag("profile", "cpu and heap profile file").Default("").String()

	pileupApp       = app.Command("pileup", "pileup reads")
	pileupMinBQ     = pileupApp.Flag("min-BQ", "minimum base quality").Short('Q').Default("13").Int()
	pileupMinMQ     = pileupApp.Flag("min-MQ", "minimum mapping quality").Short('q').Default("0").Int()
	pileupOutFile   = pileupApp.Flag("outfile", "output file").Short('o').Default("").String()
	pileupFastaFile = pileupApp.Flag("fastafile", "genome fasta file").Short('f').Default("").String()
	pileupBamFile   = pileupApp.Arg("bamfile", "bam file of reads").Required().String()

	piApp         = app.Command("pi", "calculate pi")
	piMinBQ       = piApp.Flag("min-BQ", "minimum base quality").Short('Q').Default("13").Int()
	piOutFile     = piApp.Flag("output", "output file").Short('o').Default("").String()
	piRegionStart = piApp.Flag("region-start", "region start").Default("0").Int()
	piRegionEnd   = piApp.Flag("region-end", "region end").Default("0").Int()
	piPileupFile  = piApp.Arg("pileupfile", "pileup file").Required().String()

	ctApp           = app.Command("ct", "calculate total correlation")
	ctCondonTableID = ctApp.Flag("codon", "condon table ID").Default("11").String()
	ctMaxL          = ctApp.Flag("maxl", "max length of correlation").Default("100").Int()
	ctPos           = ctApp.Flag("pos", "position").Default("4").Int()
	ctMinCoverage   = ctApp.Flag("min-coverage", "minimum read coverage").Default("10").Int()
	ctRegionStart   = ctApp.Flag("region-start", "region start").Default("0").Int()
	ctRegionEnd     = ctApp.Flag("region-end", "region end").Default("0").Int()
	ctChunckSize    = ctApp.Flag("chunck-size", "chunck size").Default("10000").Int()
	ctPileupFile    = ctApp.Arg("pileup", "pileup file").Required().String()
	ctFastaFile     = ctApp.Arg("fasta", "genome fasta file").Required().String()
	ctGffFile       = ctApp.Arg("gff", "GFF file").Required().String()
	ctOutFile       = ctApp.Arg("out", "output file").Required().String()

	crApp           = app.Command("cr", "calculate total correlation")
	crCondonTableID = crApp.Flag("codon", "condon table ID").Default("11").String()
	crMaxL          = crApp.Flag("maxl", "max length of correlation").Default("100").Int()
	crPos           = crApp.Flag("pos", "position").Default("4").Int()
	crMinCoverage   = crApp.Flag("min-coverage", "minimum read coverage").Default("10").Int()
	crRegionStart   = crApp.Flag("region-start", "region start").Default("0").Int()
	crRegionEnd     = crApp.Flag("region-end", "region end").Default("0").Int()
	crChunckSize    = crApp.Flag("chunck-size", "chunck size").Default("10000").Int()
	crPrefix        = crApp.Arg("prefix", "prefix").Required().String()
	crGenomeDir     = crApp.Arg("genome-dir", "genome directory").Required().String()
)

func main() {
	// Parse flags and arguments.
	command := kingpin.MustParse(app.Parse(os.Args[1:]))

	// Profile CPU usage.
	if *profile != "" {
		cpuprofile := *profile + ".cpu"
		heapprofile := *profile + ".heap"
		f1 := createFile(cpuprofile)
		defer f1.Close()
		f2 := createFile(heapprofile)
		defer f2.Close()
		if *debug {
			log.Printf("Created CPU profile file: %s\nand Heap profile: %s\n", cpuprofile, heapprofile)
		}
		pprof.StartCPUProfile(f1)
		pprof.WriteHeapProfile(f2)
		defer pprof.StopCPUProfile()
	}

	switch command {
	case pileupApp.FullCommand():
		pileupCmd := cmdPileup{
			minBQ:     *pileupMinBQ,
			minMQ:     *pileupMinMQ,
			outFile:   *pileupOutFile,
			fastaFile: *pileupFastaFile,
			bamFile:   *pileupBamFile,
		}
		pileupCmd.Run()
		break
	case piApp.FullCommand():
		piCmd := cmdPi{
			outFile:    *piOutFile,
			minBQ:      *piMinBQ,
			pileupFile: *piPileupFile,
			debug:      *debug,
		}
		piCmd.regionStart = *piRegionStart
		piCmd.regionEnd = *piRegionEnd
		piCmd.Run()
		break
	case ctApp.FullCommand():
		runtime.GOMAXPROCS(*ncpu)
		ctCmd := cmdCt{
			codonTableID: *ctCondonTableID,
			maxl:         *ctMaxL,
			pos:          *ctPos,
			minCoverage:  *ctMinCoverage,
			regionStart:  *ctRegionStart,
			regionEnd:    *ctRegionEnd,
			chunckSize:   *ctChunckSize,
			pileupFile:   *ctPileupFile,
			fastaFile:    *ctFastaFile,
			gffFile:      *ctGffFile,
			outFile:      *ctOutFile,
			debug:        *debug,
		}
		ctCmd.Run()
		break
	case crApp.FullCommand():
		crCmd := cmdCr{
			codonTableID: *crCondonTableID,
			maxl:         *crMaxL,
			pos:          *crPos,
			minCoverage:  *crMinCoverage,
			regionStart:  *crRegionStart,
			regionEnd:    *crRegionEnd,
			chunckSize:   *crChunckSize,
			genomeDir:    *crGenomeDir,
			prefix:       *crPrefix,
		}
		crCmd.Run()
		break
	}
}
