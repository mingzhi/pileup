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
	featDir = featApp.Arg("genome_dir", "genome directory").Required().String()
	featOut = featApp.Arg("feature_db_path", "feature db path").Required().String()

	readApp      = app.Command("read", "read samtools mpileup results.")
	pileupFile   = readApp.Arg("pileup_file", "pileup file.").Required().String()
	readFeature  = readApp.Arg("feature_db_path", "feature db path").Required().String()
	readOut      = readApp.Arg("results_db_path", "results db path").Required().String()
	readMinDepth = readApp.Flag("min_depth", "min depth").Default("5").Int()
	readMinCover = readApp.Flag("min_coverage", "min coverage").Default("0.8").Float64()

	reportApp       = app.Command("report", "report db statistics.")
	reportFeatureDB = reportApp.Arg("feature_db_path", "feature db path").Required().String()
	reportResultsDB = reportApp.Arg("results_db_path", "results db path").Required().String()
	reportPrefix    = reportApp.Flag("prefix", "prefix").Required().String()
	reportMaxl      = reportApp.Flag("maxl", "max length of correlations").Default("300").Int()

	covApp       = app.Command("cov", "calculate rate covariance.")
	covFeatureDb = covApp.Arg("feature_db_path", "feature db path").Required().String()
	covResultsDb = covApp.Arg("results_db_path", "results db file").Required().String()
	covGC        = covApp.Flag("codon", "codon table id").Default("11").String()
	covMinDepth  = covApp.Flag("min_depth", "min depth").Default("5").Int()
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
		featureDB := createLMDBEnv(*reportFeatureDB)
		defer featureDB.Close()
		resultsDB := createLMDBEnv(*reportResultsDB)
		defer resultsDB.Close()
		reportcmd := cmdReport2{
			featureDB: featureDB,
			resultsDB: resultsDB,
			prefix:    *reportPrefix,
			maxl:      *reportMaxl,
		}
		reportcmd.run()
		break
	case covApp.FullCommand():
		crcmd := cmdCr{
			dbfile:        *covResultsDb,
			codonID:       *covGC,
			featureDbPath: *covFeatureDb,
			minDepth:      *covMinDepth,
		}
		crcmd.run()
		break
	}
}
