package main

import (
	"github.com/alecthomas/kingpin"
	"os"
)

var (
	app   = kingpin.New("pcorr", "A command-line application for correlation calculation.")
	debug = app.Flag("debug", "Enable debug mode.").Bool()

	piApp        = app.Command("pi", "calculate pi")
	piMinBQ      = piApp.Flag("min-BQ", "minimum base quality").Short('Q').Int()
	piOutFile    = piApp.Flag("output", "output file").Short('o').String()
	piPileupFile = piApp.Arg("pileupfile", "pileup file").Required().String()
)

func main() {
	switch kingpin.MustParse(app.Parse(os.Args[1:])) {
	case piApp.FullCommand():
		piCmd := cmdPi{
			outFile:    *piOutFile,
			minBQ:      *piMinBQ,
			pileupFile: *piPileupFile,
		}
		piCmd.Run()
	}
}
