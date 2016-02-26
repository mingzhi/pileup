package main

import (
	"github.com/alecthomas/kingpin"
	"os"
)

var (
	app   = kingpin.New("pileup", "A command-line application for reading pileup file.")
	debug = app.Flag("debug", "Enable debug mode.").Bool()

	featApp = app.Command("feat", "read genome features.")
	featDir = featApp.Arg("dir", "features diretory.").Required().String()
	featOut = featApp.Arg("out", "features db file.").Required().String()

	readApp    = app.Command("read", "read samtools mpileup results.")
	pileupFile = readApp.Arg("pileup_file", "pileup file.").Required().String()
)

func main() {
	command := kingpin.MustParse(app.Parse(os.Args[1:]))

	switch command {
	case featApp.FullCommand():
		featcmd := featCmd{
			out: *featOut,
			dir: *featDir,
		}
		featcmd.run()
		break
	}
}
