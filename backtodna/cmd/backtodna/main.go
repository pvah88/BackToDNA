//export GOPATH=`pwd`
package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"time"
)

var (
	noLogFlag     = flag.Bool("nolog", false, "use this flag to disable log print outs")
	queryFlag     = flag.String("query", "", "your amino acid alignment as query file inluding its directory")
	outputDirFlag = flag.String("outdir", ".", "directory to amino acid alignment, DNA database and blast output ")
	geneticCode   = map[string]string{
		"TCA": "S",
		"TCC": "S",
		"TCG": "S",
		"TCT": "S",
		"TTC": "F",
		"TTT": "F",
		"TTA": "L",
		"TTG": "L",
		"TAC": "Y",
		"TAT": "Y",
		"TAA": "_",
		"TAG": "_",
		"TGC": "C",
		"TGT": "C",
		"TGA": "_",
		"TGG": "W",
		"CTA": "L",
		"CTC": "L",
		"CTG": "L",
		"CTT": "L",
		"CCA": "P",
		"CAT": "H",
		"CAA": "Q",
		"CAG": "Q",
		"CGA": "R",
		"CGC": "R",
		"CGG": "R",
		"CGT": "R",
		"ATA": "I",
		"ATC": "I",
		"ATT": "I",
		"ATG": "M",
		"ACA": "T",
		"ACC": "T",
		"ACG": "T",
		"ACT": "T",
		"AAC": "N",
		"AAT": "N",
		"AAA": "K",
		"AAG": "K",
		"AGC": "S",
		"AGT": "S",
		"AGA": "R",
		"AGG": "R",
		"CCC": "P",
		"CCG": "P",
		"CCT": "P",
		"CAC": "H",
		"GTA": "V",
		"GTC": "V",
		"GTG": "V",
		"GTT": "V",
		"GCA": "A",
		"GCC": "A",
		"GCG": "A",
		"GCT": "A",
		"GAC": "D",
		"GAT": "D",
		"GAA": "E",
		"GAG": "E",
		"GGA": "G",
		"GGC": "G",
		"GGG": "G",
		"GGT": "G",
	}
)

type blastResult struct {
	SubjectID         string
	QueryStart        int
	QueryEnd          int
	SubjectStart      int
	SubjectEnd        int
	Frame             int
	AllignedQuery     string
	AllignedSubjTrans string
	QueryID           string
	Bitscore          float64
	Pident            float64
}

type blastResults []blastResult

func (a blastResults) Len() int {
	return len(a)
}

func (a blastResults) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

func (a blastResults) Less(i, j int) bool {
	if a[i].QueryID != a[j].QueryID {
		return i < j
	}
	return a[i].Bitscore > a[j].Bitscore
}

func readNamedFiles(suffix string) map[string][]string {
	files, err := ioutil.ReadDir(*outputDirFlag)
	if err != nil {
		panic(err)
	}

	fastas := make(map[string][]string)

	for _, f := range files {
		fName := f.Name()
		if strings.HasSuffix(fName, suffix) {
			fileBytes, err := ioutil.ReadFile(filepath.Join(*outputDirFlag, fName))
			if err != nil {
				panic(err)
			}
			lines := strings.Split(string(fileBytes), "\n")
			fastas[strings.TrimSuffix(fName, suffix)] = lines
		}
	}

	return fastas
}

func readOneNamedFile(fName string) []string {
	fileBytes, err := ioutil.ReadFile(filepath.Join(*outputDirFlag, fName+".fa"))
	if err != nil {
		panic(err)
	}
	lines := strings.Split(string(fileBytes), "\n")
	return lines
}

func parseInt(s string) int {
	num, err := strconv.Atoi(strings.Trim(s, " "))
	if err != nil {
		panic(err)
	}
	return num
}

func parseFloat(s string) float64 {
	num, err := strconv.ParseFloat(strings.Trim(s, " "), 64)
	if err != nil {
		numInt, err := strconv.Atoi(strings.Trim(s, " "))
		if err != nil {
			panic(err)
		}
		return float64(numInt)
	}
	return num
}

func parseAllBlastResults(namedFileLines map[string][]string) map[string][]blastResult {

	res := make(map[string][]blastResult)

	for fName, v := range namedFileLines {
		res[fName] = make([]blastResult, 0)
		for _, resultLine := range v {
			if len(resultLine) == 0 {
				continue
			}
			outSplit := strings.Split(resultLine, ",")
			curRes := blastResult{}
			curRes.SubjectID = outSplit[0]
			curRes.QueryStart = parseInt(outSplit[1])

			curRes.QueryEnd = parseInt(outSplit[2])
			curRes.SubjectStart = parseInt(outSplit[3])
			curRes.SubjectEnd = parseInt(outSplit[4])

			curRes.Frame = parseInt(strings.Split(outSplit[5], "/")[1])
			curRes.AllignedQuery = outSplit[6]
			curRes.AllignedSubjTrans = outSplit[7]
			curRes.QueryID = outSplit[8]

			curRes.Bitscore = parseFloat(outSplit[9])
			curRes.Pident = parseFloat(outSplit[10])
			if curRes.Pident < 100.0 {
				continue
			}

			res[fName] = append(res[fName], curRes)
		}

		sort.Sort(blastResults(res[fName]))
		for i := len(res[fName]) - 1; i >= 1; i-- {
			if res[fName][i].QueryID == res[fName][i-1].QueryID {
				res[fName] = append(res[fName][:i], res[fName][i+1:]...)
			}
		}
	}

	return res
}

func queryMatch(Query []string) string {
	queryContent := strings.Join(Query, "")
	return queryContent
}

func subjectMatch(blastOut []blastResult, fname string) map[string]string {
	file, err := os.Open(filepath.Join(*outputDirFlag, fname+".fa"))
	for err != nil {
		log.Println("error opening .fa file, retrying.", fname, err.Error())
		time.Sleep(3 * time.Second)
		file, err = os.Open(filepath.Join(*outputDirFlag, fname+".fa"))
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	subjectIDToDNALine := make(map[string]string)
	blastOutMap := make(map[string]bool)

	for _, v := range blastOut {
		blastOutMap[v.SubjectID] = true
	}

	var curDna bytes.Buffer
	curSubject := ""
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if len(line) == 0 {
			continue
		}
		if (line[0] == '>') && (curSubject == "") {
			curSubject = line[1:]
		} else if (line[0] == '>') && (curSubject != "") {
			if _, ok := blastOutMap[curSubject]; ok {
				subjectIDToDNALine[curSubject] = curDna.String()
			}
			curDna.Reset()
			curSubject = line[1:]
			if len(blastOutMap) == len(subjectIDToDNALine) {
				break
			}

		} else {
			curDna.WriteString(strings.TrimSpace(line))
		}
	}
	return subjectIDToDNALine
}

func subjectSpecific(subjectMatch string, subjectStart int, subjectEnd int, frame int) []string {

	if frame == 1 {
		return codonGroup(subjectMatch[subjectStart-1 : subjectEnd])
	} else if frame == 2 {
		return codonGroup(subjectMatch[subjectStart-1 : subjectEnd])
	} else if frame == 3 {
		return codonGroup(subjectMatch[subjectStart-1 : subjectEnd])
	} else if -3 <= frame && frame <= -1 {
		complementStrand := strings.Replace(
			strings.Replace(
				strings.Replace(
					strings.Replace(
						strings.Replace(
							strings.Replace(subjectMatch, "A", "W", -1),
							"T", "A", -1),
						"W", "T", -1),
					"C", "Z", -1),
				"G", "C", -1),
			"Z", "G", -1)

		switch {
		case frame == -1:
			return codonGroup(reverseComplement(complementStrand[subjectEnd-1 : subjectStart]))
		case frame == -2:
			return codonGroup(reverseComplement(complementStrand[subjectEnd-1 : subjectStart]))
		case frame == -3:
			return codonGroup(reverseComplement(complementStrand[subjectEnd-1 : subjectStart]))
		default:
			return codonGroup("Nothing")
		}
	} else {
		return codonGroup("Ignore")
	}
}

func reverseComplement(complementStrand string) string {
	r := []byte(complementStrand)
	for i, j := 0, len(r)-1; i < len(r)/2; i, j = i+1, j-1 {
		r[i], r[j] = r[j], r[i]
	}
	return string(r)
}

func codonGroup(s string) []string {
	triples := make([]string, len(s)/3)
	for i := 0; i < len(s)/3; i++ {
		triples[i] = s[i*3 : (i+1)*3]
	}
	return triples
}

func processResults(queryFileLines []string, namedBlastResults map[string][]blastResult) {
	queryIDToQueryContent := make(map[string]string)
	for _, data := range queryFileLines[1:] {
		queryLines := strings.Split(data, "\n")
		queryIDToQueryContent[queryLines[0]] = queryMatch(queryLines[1:])
	}

	nonDashRegexp := regexp.MustCompile("[^-]+")
	alreadyPrintedQueryIds := make(map[string]bool)

	for fName, blastResults := range namedBlastResults {
		subjectToDNAStringMap := subjectMatch(blastResults, fName)
		for _, oneBlastResult := range blastResults {

			query := queryIDToQueryContent[oneBlastResult.QueryID]

			subjectMatch, ok := subjectToDNAStringMap[oneBlastResult.SubjectID]
			if !ok {
				panic(fmt.Sprintf("Missing subjectId %s\n%+ v\n", oneBlastResult.SubjectID, subjectToDNAStringMap))
			}

			codon := subjectSpecific(subjectMatch, oneBlastResult.SubjectStart,
				oneBlastResult.SubjectEnd, oneBlastResult.Frame)

			nonDashIndexes := nonDashRegexp.FindStringIndex(query)
			var buf bytes.Buffer
			adjustedQueryStart := nonDashIndexes[0] + oneBlastResult.QueryStart - 1

			if ok, _ := alreadyPrintedQueryIds[oneBlastResult.QueryID]; ok {
				continue
			} else {
				alreadyPrintedQueryIds[oneBlastResult.QueryID] = true
			}

			buf.WriteString(fmt.Sprintf(">%s|%s\n", oneBlastResult.QueryID, oneBlastResult.SubjectID))
			buf.WriteString(strings.Repeat("---", adjustedQueryStart))

			queryIndex := adjustedQueryStart
			for _, oneCodon := range codon {
				for query[queryIndex] == '-' {
					buf.WriteString("---")
					queryIndex++
				}

				buf.WriteString(oneCodon)
				queryIndex++
			}

			buf.WriteString(strings.Repeat("---", len(query)-queryIndex))
			fmt.Println(buf.String())
		}
	}
}

func main() {
	flag.Parse()
	if *noLogFlag {
		log.SetOutput(ioutil.Discard)
	}
	if *queryFlag == "" {
		flag.Usage()
		return
	}

	outs := readNamedFiles(".out")

	queryFileContent, err := ioutil.ReadFile(*queryFlag)
	if err != nil {
		panic(err)
	}

	queryFileLines := strings.Split(string(queryFileContent), ">")

	namedBlastResults := parseAllBlastResults(outs)
	processResults(queryFileLines, namedBlastResults)

}
