/*
This file is part of HiCdat.

HiCdat is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

HiCdat is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

See <http://www.gnu.org/licenses/> for a a copy of the GNU General Public License.
*/

#ifndef HELPERCLASSES_H
#define HELPERCLASSES_H

#include <QtCore>
#include <iostream>
#include <fstream>

class readPair
{
public:
    QString chromA;
    QChar strandA;
    int posA;
    QString chromB;
    QChar strandB;
    int posB;

    inline readPair() : chromA(""), strandA(QChar('.')), posA(0), chromB(""), strandB(QChar('.')), posB(0) {}

    inline void initFirst(QString& cA, QChar& sA, int& pA) {
        chromA = cA;
        strandA = sA;
        posA = pA;
        chromB = "";
        strandB = '.';
        posB = 0;
    }

    inline void addSecond(QString& cB, QChar& sB, int& pB) {
        chromB = cB;
        strandB = sB;
        posB = pB;
    }

    inline bool hasB() const{
        if (chromB == "") {return(false);}
        return(true);
    }
};

inline bool myChromosomeSorterLessThan(const QString &s1, const QString &s2) {
    // if no letter in name, this wins (is smaller) over a letter in the name
    // same for letter + number versus only letter
    // same for only number versus letter + number
    // if only letter, then normal string comparison
    // if only number, then normal number comparison
    // if letter and number, first letter comparison, then number comparison.
    QRegExp rxNumberA("(\\d+)");
    QRegExp rxLetterA("(\\D+)");
    QRegExp rxNumberB("(\\d+)");
    QRegExp rxLetterB("(\\D+)");
    int rxNumberApos = rxNumberA.indexIn(s1);
    int rxLetterApos = rxLetterA.indexIn(s1);
    int rxNumberBpos = rxNumberB.indexIn(s2);
    int rxLetterBpos = rxLetterB.indexIn(s2);
    int NOMATCH = -1;
    if ((rxNumberApos == NOMATCH) && (rxNumberBpos == NOMATCH)) { return(s1.toLower() < s2.toLower()); }
    if ((rxLetterApos == NOMATCH) && (rxLetterBpos == NOMATCH)) { return(s1.toUInt() < s2.toUInt()); }
    if ((rxLetterApos == NOMATCH) && (rxLetterBpos != NOMATCH)) { return(true); }
    if ((rxLetterApos != NOMATCH) && (rxLetterBpos == NOMATCH)) { return(false); }
    if ((rxNumberApos == NOMATCH) && (rxNumberBpos != NOMATCH)) { return(false); }
    if ((rxNumberApos != NOMATCH) && (rxNumberBpos == NOMATCH)) { return(true); }
    QString lettersA = "";
    QString lettersB = "";
    uint numberA = 0;
    uint numberB = 0;
    if ((rxLetterApos != NOMATCH) && (rxLetterBpos != NOMATCH) && (rxNumberApos == NOMATCH) && (rxNumberBpos == NOMATCH)) { return(s1.toLower() < s2.toLower()); } // only letters
    if ((rxLetterApos == NOMATCH) && (rxLetterBpos == NOMATCH) && (rxNumberApos != NOMATCH) && (rxNumberBpos != NOMATCH)) { return(s1.toUInt() < s2.toUInt()); } // only numbers
    if ((rxLetterApos != NOMATCH) && (rxLetterBpos != NOMATCH) && (rxNumberApos != NOMATCH) && (rxNumberBpos == NOMATCH)) { return(true); } // number and letter in A, only letter in B
    if ((rxLetterApos != NOMATCH) && (rxLetterBpos == NOMATCH) && (rxNumberApos != NOMATCH) && (rxNumberBpos != NOMATCH)) { return(false); } // number and letter in A, only number in B
    if ((rxLetterApos != NOMATCH) && (rxLetterBpos != NOMATCH) && (rxNumberApos == NOMATCH) && (rxNumberBpos != NOMATCH)) { return(false); } // number and letter in B, only letter in A
    if ((rxLetterApos == NOMATCH) && (rxLetterBpos != NOMATCH) && (rxNumberApos != NOMATCH) && (rxNumberBpos != NOMATCH)) { return(true); } // number and letter in B, only number in A
    if ((rxLetterApos != NOMATCH) && (rxLetterBpos != NOMATCH) && (rxNumberApos != NOMATCH) && (rxNumberBpos != NOMATCH)) { // both in both
        lettersA = rxLetterA.capturedTexts().first().toLower();
        lettersB = rxLetterB.capturedTexts().first().toLower();
        numberA = rxNumberA.capturedTexts().first().toUInt();
        numberB = rxNumberB.capturedTexts().first().toUInt();
        if (lettersA != lettersB) { return(lettersA < lettersB); };
        if (numberA != numberB) { return(numberA < numberB); };
    }
    std::cerr << "check the chromosome names - they should only consist of simple letters, numbers or a one-to-one combination of this - eg chr1, but not c1x4" << std::endl << std::flush;
    std::cerr << s1.toStdString() << "\t" << s2.toStdString() << std::endl << std::flush;
    std::cerr << lettersA.toStdString() << "\t" << lettersB.toStdString() << std::endl << std::flush;
    std::cerr << numberA << "\t" << numberB << std::endl << std::flush;
    return(false); // this are all cases which should never happen
}

/*!
 *
 *
 *
 */

class matrixLine
{
public:
    QString _chromA;
    QString _strandA;
    uint _posA;
    uint _fragA;
    QString _chromB;
    QString _strandB;
    uint _posB;
    uint _fragB;

    inline matrixLine() : _chromA(""), _strandA(""), _posA(0), _fragA(0), _chromB(""), _strandB(""), _posB(0), _fragB(0) {}
    inline void init(const QString& chromA, const QString& strandA, const QString posA, const QString& chromB, const QString& strandB, const QString posB) {
        this->_chromA = chromA;
        this->_strandA = strandA;
        this->_posA = posA.toUInt();
        this->_fragA = 0;
        this->_chromB = chromB;
        this->_strandB = strandB;
        this->_posB = posB.toUInt();
        this->_fragB = 0;
    }
};

/*!
 *
 *
 *
 */


class mini_fasta_reader
{
private:
    std::ifstream _infile;
    std::string _cur_line;
public:
    std::string _name;
    std::string _sequence;
public:
    // constructor should open the file and tell if something is wrong - read until the first line with a > at beginning
    mini_fasta_reader(std::string fastafile) {
        _infile.open(fastafile.c_str());
        if (!_infile) { std::cerr << "file not found" << std::endl << std::flush; exit(8); }
        else {
            bool searching = true;
            while (searching) {
                if (!_infile.eof()){
                    std::getline(_infile, _cur_line);
                    if (_cur_line[0] != '>') { std::cerr << "no fasta header - searching more" << _cur_line << std::endl << std::flush; }
                    else { searching = false; }
                }
                else {
                    std::cerr << "did not find any valid fasta header" << std::endl << std::flush;
                    exit(8);
                }
            }
        }
    }

    // destructor should close the file
    ~mini_fasta_reader() { _infile.close(); }

    // returns true if entry was read successfully
    inline bool read_entry(){
        bool rval = false;
        bool searching = true;
        std::string cur_seq = "";

        // read the name
        if (!_infile.eof()){
            if (_cur_line[0] != '>') { std::cerr << "no fasta entry\t" << _cur_line << std::endl << std::flush; exit(8); }
            else { _name = _cur_line.substr(1, _cur_line.find(' ')-1); rval = true; searching = true; } //! TODO skip white spaces between the > and the name
        }
        else { rval = false; searching = false; }
        // get the sequence
        while (searching) {
            if (!_infile.eof()) {
                std::getline(_infile, _cur_line);
                if (_cur_line.length() == 0) { continue; }
                if ( _cur_line[0] == '>' ) { _sequence = cur_seq; searching = false; }
                else { cur_seq.append(_cur_line); }
            }
            else { _sequence = cur_seq; searching = false;}
        }
        return(rval);
    }
};


/*!
 *
 *
 *
 */


class SIMPLEgenomeFragment
{
public:
    QString _chrom;
    uint _start;
    uint _end;
    uint _fragNum;
    inline SIMPLEgenomeFragment(QString& chrom, uint& start, uint& end, uint& fragNum) : _chrom(chrom), _start(start), _end(end), _fragNum(fragNum) {}
};


/*!
 *
 *
 *
 */

inline bool SIMPLEmyGenomeFragmentSorterLessThan(const SIMPLEgenomeFragment *A, const SIMPLEgenomeFragment *B) {
    if (A->_chrom == B->_chrom) { return(false); }
    QStringList irrelevantChromosomesIdentifiers;
    irrelevantChromosomesIdentifiers << "Mt" << "Pt" << "MtDNA" << "PtDNA" << "Un" << "ChrM" << "ChrC" << "M" << "C" << "MT" << "CP" << "PT";
    if ((!irrelevantChromosomesIdentifiers.contains(A->_chrom)) && (irrelevantChromosomesIdentifiers.contains(B->_chrom))) { return(true); }
    if ((irrelevantChromosomesIdentifiers.contains(A->_chrom)) && (!irrelevantChromosomesIdentifiers.contains(B->_chrom))) { return(false); }
    return(myChromosomeSorterLessThan(A->_chrom, B->_chrom));
}

class SIMPLEgenomeFragmentCollection
{
private:
    //! the rootItem of the database and the index
    QList<SIMPLEgenomeFragment*> genomeFragments;
    QHash< QPair<QString,uint>, QVector<SIMPLEgenomeFragment*> > CPindex;
    QHash< QString, uint> chromToSize; // for the organism specific R-file
    QHash< QString, uint> chromToFragNum; // for the organism specific R-file
    QStringList relevantChromosomes; // for the organism specific R-file
    QStringList irrelevantChromosomes; // for the organism specific R-file
    QStringList irrelevantChromosomesIdentifiers; // for the organism specific R-file
    QList<uint> chromosomesWithNumbers; // for the organism specific R-file
    QStringList chromosomesWithLetters; // for the organism specific R-file
    QString fastaFile;
    QString fragInfile; // instead of the fasta file
    QString restrictionPatternOrWinSize;
    uint stepSize;

private:
    inline uint getOffset(const uint& position) {
        uint x = (position / this->stepSize); // uses integer division
        return(x);
    }

    void createIndex() {
        QString chrom;
        uint start;
        uint end;
        uint startStep;
        uint endStep;
        uint stepCounter;
        QPair<QString,uint> addKey;
        QVector<SIMPLEgenomeFragment*> emptyVec;
        emptyVec.reserve(100);
        foreach (SIMPLEgenomeFragment* frag, this->genomeFragments) {
            chrom = frag->_chrom;
            start = frag->_start;
            end = frag->_end;
            startStep = this->getOffset(start);
            endStep = this->getOffset(end);
            for (stepCounter = startStep; stepCounter <= endStep; ++stepCounter) {
                addKey = qMakePair(chrom, stepCounter);
                if (!this->CPindex.contains(addKey)) { this->CPindex.insert(addKey, emptyVec); }
                this->CPindex[addKey].append(frag);
            }
        }
    }

    bool setupUsingWindowSize() {
        bool rval = false;
        QVector<qint64> matches;
        QString chrom = "";
        uint start = 0;
        uint end = 0;
        uint fragNumber = 0;
        uint windowSize = this->restrictionPatternOrWinSize.toUInt();
        // open a fasta reader to get the sizes of the chromosomes
        mini_fasta_reader fasta(this->fastaFile.toStdString());
        while (fasta.read_entry()) { //! NOTE THAT THE HEADERS SHOULD BE >CHR1
            chrom = QString::fromStdString(fasta._name);
            // well - not an optimal solution, but ok
            matches.clear();
            matches.reserve(50000);
            for (int i = 0; i*windowSize < fasta._sequence.length(); ++i) {
                matches.push_back(i*windowSize);
            }
            matches.push_back(fasta._sequence.length()-1);
            // add restriction fragments
            for (int i = 0; i < matches.size()-1; ++i) {
                start = matches.at(i);
                end = matches.at(i+1);
                SIMPLEgenomeFragment *frag = new SIMPLEgenomeFragment(chrom, start, end, fragNumber); //! NOTE THAT THE NUMBER NEEDS TO BE ADDED LATER ON
                this->genomeFragments += frag;
            }
            std::cerr << "added " << (matches.size()-1) << " windows of chromosome " << chrom.toStdString() << " with size " << fasta._sequence.length() << std::endl << std::flush; //! TEST
            rval = true;
        }
        //! very important - sort the fragments according to the chromosome - use qStableSort to keep the "within-chromosome" order correct
        qStableSort(this->genomeFragments.begin(), this->genomeFragments.end(), SIMPLEmyGenomeFragmentSorterLessThan);
        foreach (SIMPLEgenomeFragment * frag, this->genomeFragments) {
            ++fragNumber;
            frag->_fragNum = fragNumber;
        }
        return(rval);
    }

    bool setupUsingOneRestrictionSite() {
        bool rval = false;
        QVector<qint64> matches;
        QString chrom = "";
        uint start = 0;
        uint end = 0;
        qint64 currentMatch = -1;
        uint fragNumber = 0;
        std::string resPat = this->restrictionPatternOrWinSize.toStdString();
        // open a fasta reader to get chromosome sequences
        mini_fasta_reader fasta(this->fastaFile.toStdString());
        while (fasta.read_entry()) { //! NOTE THAT THE HEADERS SHOULD BE >CHR1
            chrom = QString::fromStdString(fasta._name);
            // get first all the restriction pattern matches
            matches.clear();
            matches.reserve(50000);
            matches.push_back(0);
            currentMatch = fasta._sequence.find(resPat);
            while (currentMatch != -1) {
                matches.push_back(currentMatch);
                currentMatch = fasta._sequence.find(resPat, currentMatch+1);
            }
            matches.push_back(fasta._sequence.length()-1); //! NOTE THAT THE LEFT RESTRICTION SITE IS INSIDE, THE RIGHT OUTSIDE >> except here for the last fragment
            // add restriction fragments
            for (int i = 0; i < matches.size()-1; ++i) {
                start = matches.at(i);
                end = matches.at(i+1);
                SIMPLEgenomeFragment *frag = new SIMPLEgenomeFragment(chrom, start, end, fragNumber);//! NOTE THAT THE NUMBER NEEDS TO BE ADDED LATER ON
                this->genomeFragments += frag;
            }
            std::cerr << "added " << (matches.size()-1) << " restriction fragments of chromosome " << chrom.toStdString() << " with size " << fasta._sequence.length() << std::endl << std::flush; //! TEST
            rval = true;
            // for the organism specific R-file
            this->chromToSize[chrom] = fasta._sequence.length();
            this->chromToFragNum[chrom] = matches.size()-1;
            if (this->irrelevantChromosomesIdentifiers.contains(chrom)) { this->irrelevantChromosomes << chrom; }
            else {this->relevantChromosomes << chrom;}
        }
        //! very important - sort the fragments according to the chromosome - use qStableSort to keep the "within-chromosome" order correct
        qStableSort(this->genomeFragments.begin(), this->genomeFragments.end(), SIMPLEmyGenomeFragmentSorterLessThan);
        foreach (SIMPLEgenomeFragment * frag, this->genomeFragments) {
            ++fragNumber;
            frag->_fragNum = fragNumber;
        }
        return(rval);
    }

    bool setupUsingTwoRestrictionSites() {
        bool rval = false;
        QVector<qint64> matches;
        QString chrom = "";
        uint start = 0;
        uint end = 0;
        qint64 currentMatch = -1;
        qint64 currentMatchA = -1;
        qint64 currentMatchB = -1;
        uint fragNumber = 0;
        QStringList resPatterns = this->restrictionPatternOrWinSize.split(",");
        std::string resPatA = resPatterns.at(0).toStdString();
        std::string resPatB = resPatterns.at(1).toStdString();
        // open a fasta reader to get chromosome sequences
        mini_fasta_reader fasta(this->fastaFile.toStdString());
        while (fasta.read_entry()) { //! NOTE THAT THE HEADERS SHOULD BE >CHR1
            chrom = QString::fromStdString(fasta._name);
            // get first all the restriction pattern matches
            matches.clear();
            matches.reserve(50000);
            matches.push_back(0);
            currentMatchA = fasta._sequence.find(resPatA);
            currentMatchB = fasta._sequence.find(resPatB);
            if ((currentMatchA == -1) && (currentMatchB == -1)) { currentMatch = -1; }
            else if ((currentMatchA != -1) && (currentMatchB == -1)) { currentMatch = currentMatchA; }
            else if ((currentMatchA == -1) && (currentMatchB != -1)) { currentMatch = currentMatchB; }
            else {
                if (currentMatchA < currentMatchB) { currentMatch = currentMatchA; }
                else {currentMatch = currentMatchB; }
            }
            while (currentMatch != -1) {
                matches.push_back(currentMatch);
                currentMatchA = fasta._sequence.find(resPatA, currentMatch+1);
                currentMatchB = fasta._sequence.find(resPatB, currentMatch+1);
                if ((currentMatchA == -1) && (currentMatchB == -1)) { currentMatch = -1; }
                else if ((currentMatchA != -1) && (currentMatchB == -1)) { currentMatch = currentMatchA; }
                else if ((currentMatchA == -1) && (currentMatchB != -1)) { currentMatch = currentMatchB; }
                else {
                    if (currentMatchA < currentMatchB) { currentMatch = currentMatchA; }
                    else {currentMatch = currentMatchB; }
                }
            }
            matches.push_back(fasta._sequence.length()-1); //! NOTE THAT THE LEFT RESTRICTION SITE IS INSIDE, THE RIGHT OUTSIDE >> except here for the last fragment
            // add restriction fragments
            for (int i = 0; i < matches.size()-1; ++i) {
                start = matches.at(i);
                end = matches.at(i+1);
                SIMPLEgenomeFragment *frag = new SIMPLEgenomeFragment(chrom, start, end, fragNumber);//! NOTE THAT THE NUMBER NEEDS TO BE ADDED LATER ON
                this->genomeFragments += frag;
            }
            std::cerr << "added " << (matches.size()-1) << " restriction fragments of chromosome " << chrom.toStdString() << " with size " << fasta._sequence.length() << std::endl << std::flush; //! TEST
            rval = true;
            // for the organism specific R-file
            this->chromToSize[chrom] = fasta._sequence.length();
            this->chromToFragNum[chrom] = matches.size()-1;
            if (this->irrelevantChromosomesIdentifiers.contains(chrom)) { this->irrelevantChromosomes << chrom; }
            else {this->relevantChromosomes << chrom;}
        }
        //! very important - sort the fragments according to the chromosome - use qStableSort to keep the "within-chromosome" order correct
        qStableSort(this->genomeFragments.begin(), this->genomeFragments.end(), SIMPLEmyGenomeFragmentSorterLessThan);
        foreach (SIMPLEgenomeFragment * frag, this->genomeFragments) {
            ++fragNumber;
            frag->_fragNum = fragNumber;
        }
        return(rval);
    }

    bool setupUsingFragFile() {
        bool rval = true;

        QFile infile(fragInfile);
        if (!infile.open(QIODevice::ReadOnly)) { return(false); }
        QTextStream in(&infile);
        in.setCodec("UTF-8");

        QString line;
        QStringList fields;
        QString chrom = "";
        uint start = 0;
        uint end = 0;
        uint fragNum = 0;

        line = in.readLine(); //! there is a header which needs to be ignored

        while ( !in.atEnd() ) {
            line = in.readLine();
            fields = line.split("\t");
            chrom = fields.at(1);
            start = fields.at(2).toUInt();
            end = fields.at(3).toUInt();
            fragNum = fields.at(0).toUInt();
            SIMPLEgenomeFragment *frag = new SIMPLEgenomeFragment(chrom, start, end, fragNum);
            this->genomeFragments += frag;
        }

        infile.close();
        if (infile.error() != QFile::NoError) {
            std::cerr << "Error: Cannot read file " << qPrintable(fragInfile)
                      << ": " << qPrintable(infile.errorString())
                      << std::endl;
            rval = false;
        }

        return(rval);
    }


public:
    explicit SIMPLEgenomeFragmentCollection(QString fastaFile, QString restrictionPatternOrWinSize, uint stepSize) :
        fastaFile(fastaFile),
        fragInfile(""),
        restrictionPatternOrWinSize(restrictionPatternOrWinSize),
        stepSize(stepSize) { irrelevantChromosomesIdentifiers << "Mt" << "Pt" << "MtDNA" << "PtDNA" << "Un" << "ChrM" << "ChrC" << "M" << "C" << "MT" << "CP" << "PT"; }

    explicit SIMPLEgenomeFragmentCollection(QString fragFile, uint stepSize) :
        fastaFile(""),
        fragInfile(fragFile),
        restrictionPatternOrWinSize(""),
        stepSize(stepSize) { irrelevantChromosomesIdentifiers << "Mt" << "Pt" << "MtDNA" << "PtDNA" << "Un" << "ChrM" << "ChrC" << "M" << "C" << "MT" << "CP" << "PT"; }

    ~SIMPLEgenomeFragmentCollection() { qDeleteAll(genomeFragments); }

    bool setup() {
        QRegExp rx("(\\d+)");
        bool rval = true;
        if (fragInfile != "") {
            rval = this->setupUsingFragFile();
        } else {
            if (restrictionPatternOrWinSize.contains(rx)) {
                rval = this->setupUsingWindowSize();
            } else {
                if (restrictionPatternOrWinSize.contains(",")) {
                    rval = this->setupUsingTwoRestrictionSites();
                } else {
                    rval = this->setupUsingOneRestrictionSite();
                }
            }
        }
        if (rval) { this->createIndex(); }
        return(rval);
    }

    bool writeFragments(QString fragListFile) {
        QFile outfile(fragListFile);
        if (!outfile.open(QFile::WriteOnly | QFile::Text)) { return(false); }
        QTextStream out(&outfile);
        out.setCodec("UTF-8");

        out << "fragmentNumber\tchrom\tstart\tend\n";
        foreach (SIMPLEgenomeFragment* frag, this->genomeFragments) {
            out << frag->_fragNum << "\t" <<
                   frag->_chrom << "\t" <<
                   frag->_start << "\t" <<
                   frag->_end << "\n";
        }

        out.flush();
        outfile.close();
        if (outfile.error() != QFile::NoError) {
            std::cerr << "Error: Cannot write file " << qPrintable(fragListFile)
                      << ": " << qPrintable(outfile.errorString())
                      << std::endl;
            return(false);
        }


        return(true);
    }

    bool writeOrganismFileForR(QString outfileName) {
        // check if it is a restriction site - if not, interrupt it!
        QRegExp rx("(\\d+)");
        if (restrictionPatternOrWinSize.contains(rx)) { std::cerr << "do not try to write an R-file without providing a restriction site" << std::endl << std::flush; return(false); }

        QFile outfile(outfileName);
        if (!outfile.open(QFile::WriteOnly | QFile::Text)) { return(false); }
        QTextStream out(&outfile);
        out.setCodec("UTF-8");

        //! sort the chromosomes in a stable and senseful way
        qSort(this->relevantChromosomes.begin(), this->relevantChromosomes.end(), myChromosomeSorterLessThan);
        qSort(this->irrelevantChromosomes.begin(), this->irrelevantChromosomes.end(), myChromosomeSorterLessThan);
        QStringList chromOrder;
        chromOrder << this->relevantChromosomes << this->irrelevantChromosomes;

        QString relevantString = "";
        QTextStream(&relevantString) << "'" << this->relevantChromosomes.join("', '") << "'";
        QString chromSizeString = "";
        foreach (const QString &chrom, chromOrder) {
            if (chromSizeString.isEmpty()) {
                QTextStream(&chromSizeString) << "\n\t\t'" << chrom << "' = " << this->chromToSize.value(chrom);
            } else {
                QTextStream(&chromSizeString) << ",\n\t\t'" << chrom << "' = " << this->chromToSize.value(chrom);
            }
        }

        QString fragListString = "";
        uint prevEnd = 0;
        uint curS = 0;
        uint curE = 0;
        foreach (const QString &chrom, chromOrder) {
            curS = prevEnd + 1;
            curE = curS + this->chromToFragNum.value(chrom) - 1;
            if (fragListString.isEmpty()) {
                QTextStream(&fragListString) << "\n\t\t'" << chrom << "' = c(" << curS << ", " << curE << ")";
            } else {
                QTextStream(&fragListString) << ",\n\t\t'" << chrom << "' = c(" << curS << ", " << curE << ")";
            }
            prevEnd = curE;
        }
        //QTextStream(&fragListString) << ",\n\t\t'ALL' = c(1," << prevEnd << ")";

        out << "### organism specific functions for the HiCdat\n### automatically generated by HiCdat\n\n";
        out << "f.get.chrom.sizes <- function() {\n\tout <- list(";
        out << chromSizeString;
        out << ")\n\treturn(out)\n}\n\n";
        out << "f.get.frag.list <- function() {\n\tout <- list(";
        out << fragListString;
        out << ")\n\treturn(out)\n}\n\n";
        out << "f.get.relevant.chromosomes <- function() {\n\tout <- c(";
        out << relevantString;
        out << ")\n\treturn(out)\n}\n\n";

        out.flush();
        outfile.close();
        if (outfile.error() != QFile::NoError) {
            std::cerr << "Error: Cannot write file " << qPrintable(outfileName)
                      << ": " << qPrintable(outfile.errorString())
                      << std::endl;
            return(false);
        }

        return(true);

    }

    //! NOTE THAT FRAGMENTS ARE PER SE NON-OVERLAPPING - THUS THERE IS ALWAYS ONLY ONE TO REPORT - otherwise see database for RNA-Seq

    void mapPosition(const QString &chrom, const uint &position, uint &result) {
        QPair<QString,uint> offset = qMakePair(chrom, this->getOffset(position));
        foreach (SIMPLEgenomeFragment* frag, this->CPindex.value(offset)) {
            if ((frag->_start <= position) && (frag->_end >= position)) {
                result = frag->_fragNum;
                break;
            }
        }
    }

    void mapRange(const QString &chrom, const uint &start, const uint &end, uint &result) {
        QPair<QString,uint> offset = qMakePair(chrom, this->getOffset(start));
        foreach (SIMPLEgenomeFragment* frag, this->CPindex.value(offset)) {
            if ((frag->_start <= start) && (frag->_end >= end)) {
                result = frag->_fragNum;
                break;
            }
        }
    }

    /*
    uint mapPosition(const QString &chrom, const uint &position) {
        uint out = 0;
        QPair<QString,uint> offset = qMakePair(chrom, this->getOffset(position));
        foreach (SIMPLEgenomeFragment* frag, this->CPindex.value(offset)) {
            if ((frag->_start <= position) && (frag->_end >= position)) {
                out = frag->_fragNum;
                break;
            }
        }
        return(out);
    }

    uint mapRange(const QString &chrom, const uint &start, const uint &end) {
        uint out = 0;
        QPair<QString,uint> offset = qMakePair(chrom, this->getOffset(start));
        foreach (SIMPLEgenomeFragment* frag, this->CPindex.value(offset)) {
            if ((frag->_start <= start) && (frag->_end >= end)) {
                out = frag->_fragNum;
                break;
            }
        }
        return(out);
    }
    */
};

#endif // HELPERCLASSES_H
