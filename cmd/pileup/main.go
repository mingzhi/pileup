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
	readFeature  = readApp.Arg("feature_db_path", "feature db path").Required().String()
	readMinDepth = readApp.Flag("min_depth", "min depth").Default("5").Int()
	readMinCover = readApp.Flag("min_coverage", "min coverage").Default("0.8").Float64()

	reportApp    = app.Command("report", "report db statistics.")
	reportDB     = reportApp.Arg("db", "db file").Required().String()
	reportPrefix = reportApp.Flag("prefix", "prefix").Required().String()

	covApp       = app.Command("cov", "calculate rate covariance.")
	covDB        = covApp.Arg("db", "db file").Required().String()
	covFeatureDb = covApp.Arg("feature_db_path", "feature db path").Required().String()
	covGC        = covApp.Flag("codon", "codon table id").Default("11").String()
)

func main() {
	command := kingpin.MustParse(app.Parse(os.Args[1:]))
	runtime.GOMAXPROCS(*ncpu)

	switch command {
	case featApp.FullCommand():
		featcmd := cmdFeat{
			out: *featOut,
			dir: *featDir,
		}
		featcmd.run()
		break
	case readApp.FullCommand():
		readcmd := cmdRead{
			pileupFile: *pileupFile,
			dbfile:     *readOut,
			minCover:   *readMinCover,
			minDepth:   *readMinDepth,
			featureDB:  *readFeature,
		}
		readcmd.run()
		break
	case reportApp.FullCommand():
		reportcmd := cmdReport{
			dbfile: *reportDB,
			prefix: *reportPrefix,
		}
		reportcmd.run()
		break
	case covApp.FullCommand():
		crcmd := cmdCr{
			dbfile:        *covDB,
			codonID:       *covGC,
			featureDbPath: *covFeatureDb,
		}
		crcmd.run()
		break
	}
}
