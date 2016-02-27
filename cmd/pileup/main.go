package main

import (
	"github.com/alecthomas/kingpin"
	"os"
	"runtime"
)

var (
	app   = kingpin.New("pileup", "A command-line application for reading pileup file.")
	debug = app.Flag("debug", "Enable debug mode.").Bool()
	ncpu  = app.Flag("ncpu", "number of CPUs for using.").Default("0").Int()

	featApp = app.Command("feat", "read genome features.")
	featDir = featApp.Arg("dir", "features diretory.").Required().String()
	featOut = featApp.Arg("db", "db file.").Required().String()

	readApp      = app.Command("read", "read samtools mpileup results.")
	pileupFile   = readApp.Arg("pileup_file", "pileup file.").Required().String()
	readOut      = readApp.Arg("db", "db file.").Required().String()
	readMinDepth = readApp.Flag("min_depth", "min depth").Default("10").Int()
	readMinCover = readApp.Flag("min_coverage", "min coverage").Default("0.8").Float64()
)

func main() {
	command := kingpin.MustParse(app.Parse(os.Args[1:]))
	runtime.GOMAXPROCS(*ncpu)

	switch command {
	case featApp.FullCommand():
		featcmd := featCmd{
			out: *featOut,
			dir: *featDir,
		}
		featcmd.run()
		break
	case readApp.FullCommand():
		readcmd := cmdRead{
			pileupFile: *pileupFile,
			dbfile:     *readOut,
		}
		readcmd.run()
		break
	}
}
