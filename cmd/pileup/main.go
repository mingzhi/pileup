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

	filterApp     = app.Command("filter", "filter reads.")
	filterFna     = filterApp.Arg("fna_file", "reference genome file.").Required().String()
	filterBam     = filterApp.Arg("bam_file", "bam file.").Required().String()
	filterOut     = filterApp.Arg("out_file", "out file.").Required().String()
	filterMaxDist = filterApp.Flag("max_dist", "max distance.").Default("0.05").Float64()
	filterMapQ    = filterApp.Flag("mapQ", "min mapQ").Default("30").Int()

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

	report2App       = app.Command("report2", "report db statistics.")
	report2FeatureDB = report2App.Arg("feature_db_path", "feature db path").Required().String()
	report2ResultsDB = report2App.Arg("results_db_path", "results db path").Required().String()
	report2Prefix    = report2App.Flag("prefix", "prefix").Required().String()

	covApp       = app.Command("cov", "calculate rate covariance.")
	covFeatureDb = covApp.Arg("feature_db_path", "feature db path").Required().String()
	covResultsDb = covApp.Arg("results_db_path", "results db file").Required().String()
	covGC        = covApp.Flag("codon", "codon table id").Default("11").String()
	covMinDepth  = covApp.Flag("min_depth", "min depth").Default("5").Int()

	mergeApp        = app.Command("merge", "merge mutliple mdb.")
	mergeSampleFile = mergeApp.Arg("sample", "sample list file").Required().String()
	mergeDbiName    = mergeApp.Flag("dbi", "dbi name").Default("cr").String()
	mergeOutDb      = mergeApp.Arg("out", "out db path").Required().String()
)

func main() {
	command := kingpin.MustParse(app.Parse(os.Args[1:]))
	runtime.GOMAXPROCS(*ncpu)

	switch command {
	case filterApp.FullCommand():
		filtercmd := cmdFilter{
			fnaFile:     *filterFna,
			bamFile:     *filterBam,
			outFile:     *filterOut,
			maxDistance: *filterMaxDist,
			mapQ:        *filterMapQ,
		}
		filtercmd.run()
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
		featureDB := createNoLockEnv(*reportFeatureDB)
		defer featureDB.Close()
		resultsDB := createReadOnlyEnv(*reportResultsDB)
		defer resultsDB.Close()
		reportcmd := cmdReport{
			featureDB: featureDB,
			resultsDB: resultsDB,
			prefix:    *reportPrefix,
			maxl:      *reportMaxl,
		}
		reportcmd.run()
		break
	case report2App.FullCommand():
		featureDB := createNoLockEnv(*report2FeatureDB)
		defer featureDB.Close()
		resultsDB := createReadOnlyEnv(*report2ResultsDB)
		defer resultsDB.Close()
		reportcmd2 := cmdReport2{
			featureDB: featureDB,
			resultsDB: resultsDB,
			prefix:    *report2Prefix,
		}
		reportcmd2.run()
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
	case mergeApp.FullCommand():
		mergecmd := cmdMerge{
			sampleFile: *mergeSampleFile,
			dbiName:    *mergeDbiName,
			dbOut:      *mergeOutDb,
		}
		mergecmd.run()
		break
	}
}
