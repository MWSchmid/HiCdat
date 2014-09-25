/*
This file is part of HiCat.

HiCat is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

HiCat is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

See <http://www.gnu.org/licenses/> for a a copy of the GNU General Public License.
*/

#include "workprocessor.h"
#include <iostream>
//! general functions

#include <QtGui>
#ifdef Q_WS_MAC
#include <QtWidgets> // Qt5 instead of 4
#endif

#define SEQAN_HAS_ZLIB 1
#include <zlib.h>

#include "../seqan/stream.h"
#include "../seqan/bam_io.h"

#include "annotation_reader.h"


workProcessor::workProcessor()
{
    this->_allFine = true;
    this->_isCancelled = false;
    this->_progress = 0;
    this->start();
}

workProcessor::~workProcessor()
{
    // empty the queue of projects and add STOPTHREAD - this tells the tread to exit the forever loop
    {
        QMutexLocker locker(&this->_mutex);
        while (!this->_work.isEmpty()) {
            this->_work.dequeue();
        }
        this->_work.enqueue("STOPTHREAD");
        this->_workAdded.wakeOne();
    }
    // wait before calling the base class destructor
    this->wait();
}

void workProcessor::addWork(QString workString) // will add some work to the queue
{
    QMutexLocker locker(&this->_mutex);
    this->_work.enqueue(workString);
    this->_workAdded.wakeOne();
}

void workProcessor::run()
{
    QString workString;

    forever {
        {// check if there is something to process
            QMutexLocker locker(&this->_mutex);
            if (this->_work.isEmpty()) {
                this->_workAdded.wait(&this->_mutex);
            }
            // take a project from the queue
            workString = this->_work.dequeue();
            // here is a keyword that cancels the run and stops everything
            if (workString == "STOPTHREAD") {
                break;
            }
        }

        emit this->processStatus("starting");
        emit this->processProgress(static_cast<int>(this->_progress));
        this->_allFine = this->doProcess(workString);
        if (!this->_allFine) {
            if (this->_isCancelled) {
                this->_isCancelled = false;
                this->_allFine = true;
            } else {
                this->_errorMessage = "somewhere an error";
                emit this->errorMessage(this->_errorMessage);
            }
        }
        // for the queue display one could send now the workString back
        emit this->workFinished(workString);
        emit this->processStatus("idle");
        emit this->processProgress(0);
        // check if all are processed
        {
            QMutexLocker locker(&this->_mutex);
            if (this->_work.isEmpty()) {
                emit this->idleAgain();
            }
        }
    }
}

void workProcessor::cancelProcessing() {
    this->_isCancelled = true;
}

bool workProcessor::doProcess(QString workString) {
    bool rval = false;
    QStringList fields = workString.split("|");
    QString subPart = fields.takeFirst();
    QString newWorkString = fields.join("|");
    if (subPart == "RUNmergeReadPairs") {rval = this->MERGEPAIRSdoProcess(newWorkString);}
    else if (subPart == "RUNmapReadPairs") {rval = this->MAPREADPAIRSdoProcess(newWorkString);}
    else if (subPart == "RUNcreateFragments") {rval = this->CREATEFRAGMENTSdoProcess(newWorkString);}
    else if (subPart == "RUNaddTracks") {rval = this->ADDTRACKSdoProcess(newWorkString);}
    else if (subPart == "RUNcreateRcode") {rval = this->CREATERCODEdoProcess(newWorkString);}
    else {std::cerr << "unknown work " << subPart.toStdString() << std::endl << std::flush;}
    return(rval);
}

quint64 workProcessor::getNumberLinesFromTxt(QString fileName)
{
    QFile file;
    QTextStream in;
    file.setFileName(fileName);
    quint64 numLines = 0;

    // check if file available (interrupt if not)
    if (!file.open(QIODevice::ReadOnly)) {return(0);}

    // set some stuff
    in.setDevice(&file);
    in.setCodec("UTF-8");

    // read the file
    QString line;
    while ( !in.atEnd() ) {
        line = in.readLine();
        ++numLines;
        if ((numLines % 1000000) == 0) {
            if (this->_isCancelled) { break; }
        }
    }
    // close and check for errors
    file.close();
    if (file.error() != QFile::NoError) {return(0);}
    return(numLines);
}

quint64 workProcessor::getNumberLinesFromBam(QString fileName)
{
    seqan::BamStream reader(fileName.toStdString().c_str());
    seqan::BamAlignmentRecord record;
    quint64 numLines = 0;
    while (!atEnd(reader)) {
        seqan::readRecord(record, reader);
        if (seqan::hasFlagUnmapped(record)) { continue; }
        ++numLines;
        if (this->_isCancelled) { break; }
    }
    seqan::close(reader);
    return(numLines);
}

QMap<uint, QString> workProcessor::getIDtoNAMEfromBam(QString fileName)
{
    QMap<uint, QString> out;
    std::string chromName = "*";
    out[-1] = QString::fromStdString(chromName);
    seqan::BamStream reader(fileName.toStdString().c_str());
    for (uint chromID = 0; chromID < length(reader.header.sequenceInfos); ++chromID) {
        chromName = seqan::toCString(reader.header.sequenceInfos[chromID].i1);
        out[chromID] = QString::fromStdString(chromName);
    }
    seqan::close(reader);
    return(out);
}


//! functions previously implemented in the readprocessor class


void workProcessor::MERGEPAIRSaddPartner(QString fileName, QMap<std::string, readPair> &readCollection, QMap<uint, QString> &idToName, quint64 numLines)
{
    seqan::BamStream reader(fileName.toStdString().c_str());
    seqan::BamAlignmentRecord record;
    QChar strand;
    std::string readName;
    QString chrom;
    quint64 numAdded = 0;
    while (!atEnd(reader)) {
        seqan::readRecord(record, reader);
        if (seqan::hasFlagUnmapped(record)) { continue; }
        if (seqan::hasFlagRC(record)) {strand = '-';}
        else {strand = '+';}
        chrom = idToName[record.rID];
        readName = seqan::toCString(record.qName);
        if (readCollection.count(readName) > 0) {
            ++numAdded;
            readCollection[readName].addSecond(chrom, strand, record.beginPos);
            if ((numAdded % 1000000) == 0) {
                this->_progress += (100000000.0) / (numLines*2);
                emit processProgress(static_cast<int>(this->_progress));
            }
        }
    }
    seqan::close(reader);
}

void workProcessor::MERGEPAIRSwriteCollection(QTextStream& out, QMap<std::string, readPair> &readCollection)
{
    QMapIterator<std::string, readPair> iter(readCollection);
    while(iter.hasNext()) {
        iter.next();
        if (iter.value().hasB()) {
            out << iter.value().chromA << "\t" <<
                   iter.value().strandA << "\t" <<
                   iter.value().posA << "\t" <<
                   iter.value().chromB << "\t" <<
                   iter.value().strandB << "\t" <<
                   iter.value().posB << "\n";
        }
    }
    readCollection.clear();
}

bool workProcessor::MERGEPAIRSdoProcess(QString workString)
{
    QStringList temp = workString.split("|");
    QString infileNameR1 = temp.at(0);
    QString infileNameR2 = temp.at(1);
    QString outfileName = temp.at(2);
    quint64 numMaxReads = temp.at(3).toULongLong();
    QMap<uint, QString> idToName;
    quint64 numLines;
    QChar strand;
    std::string readName;
    QString chrom;
    QMap<std::string, readPair> readCollection;
    readPair rp;
    bool rval = true;
    this->_progress = 0;
    emit processProgress(static_cast<int>(this->_progress));

    //! get line number for the process status
    emit this->processStatus("reading line numbers");
    numLines = this->getNumberLinesFromBam(infileNameR1);
    //std::cerr << numLines << std::endl << std::flush; //! TEST
    if (this->_isCancelled) {
        rval = false;
    } else {
        //! get the id to name map
        emit this->processStatus("obtaining map");
        idToName = this->getIDtoNAMEfromBam(infileNameR1);

        //! load a certain number of reads and then get their partner
        QFile file(outfileName);
        if (!file.open(QFile::WriteOnly | QFile::Text)) {
            std::cerr << "Error: Cannot read file " << qPrintable(outfileName)
                      << ": " << qPrintable(file.errorString())
                      << std::endl;
            rval = false;
        }
        else {
            QTextStream out(&file);
            out.setCodec("UTF-8");

            seqan::BamStream reader(infileNameR1.toStdString().c_str());
            seqan::BamAlignmentRecord record;
            quint64 numRead = 0;
            emit processStatus("reading forward reads");
            while (!atEnd(reader)) {
                seqan::readRecord(record, reader);
                if (seqan::hasFlagUnmapped(record)) { continue; }
                ++numRead;
                // ADD FIRST READ
                if (seqan::hasFlagRC(record)) {strand = '-';}
                else {strand = '+';}
                chrom = idToName[record.rID];
                readName = seqan::toCString(record.qName);
                rp.initFirst(chrom, strand, record.beginPos);
                readCollection[readName] = rp;
                // if many reads read, get their partner and write them out
                if ((numRead % numMaxReads) == 0) {
                    // switch the status
                    emit processStatus("reading reverse reads");
                    // add the second pair
                    this->MERGEPAIRSaddPartner(infileNameR2, readCollection, idToName, numLines);
                    // iterate over map and write out - then remove all the reads
                    this->MERGEPAIRSwriteCollection(out, readCollection);
                    // signal for progress
                    emit processStatus("reading forward reads");
                }
                if ((numRead % 1000000) == 0) {
                    // signal for progress
                    this->_progress += (100000000.0) / (numLines*2);
                    emit processProgress(static_cast<int>(this->_progress));
                    //std::cerr << numRead << std::endl << std::flush; //! TEST
                    //std::cerr << this->_progress << std::endl << std::flush; //! TEST
                    //std::cerr << static_cast<int>(this->_progress) << std::endl << std::flush; //! TEST
                    if (this->_isCancelled) {
                        readCollection.clear();
                        rval = false;
                        break;
                    }
                }
            }
            seqan::close(reader);
            if (!readCollection.isEmpty()) {
                emit processStatus("reading reverse reads");
                this->MERGEPAIRSaddPartner(infileNameR2, readCollection, idToName, numLines);
                this->MERGEPAIRSwriteCollection(out, readCollection);
            }

            // close the file and check if it worked
            out.flush();
            file.close();
            if (file.error() != QFile::NoError) {
                std::cerr << "Error: Cannot write file " << qPrintable(outfileName)
                          << ": " << qPrintable(file.errorString())
                          << std::endl;
                rval = false;
            }
        }
    }
    return(rval);
}

//! functions previously implemented in the mapprocessor class

bool workProcessor::MAPREADPAIRSdoProcess(QString workString)
{
    QStringList temp = workString.split("|");
    QString fragListFile = temp.at(0);
    QString preMatrixFile = temp.at(1);
    QString postMatrixFile = temp.at(2);
    QString doReducedMatrixStr = temp.at(3);
    QString reducedMatrixFile = temp.at(4);
    QString useFilterStr = temp.at(5);
    QString filterSelfCircleStr = temp.at(6);
    uint inwardClose = temp.at(7).toUInt();
    uint outwardClose = temp.at(8).toUInt();
    bool doReducedMatrix = false;
    bool useFilter = false;
    bool filterSelfCircle = false;
    if (doReducedMatrixStr == "y") { doReducedMatrix = true; }
    if (useFilterStr == "y") { useFilter = true; }
    if (filterSelfCircleStr == "y") { filterSelfCircle = true; }
    bool rval = false;
    quint64 numLines = 0;
    quint64 numRead = 0;
    uint readDist = 0; // this is the distance between two read ends - for the filter. It is safer to use C = A > B ? A-B : B-A; and if (C < xyz) than abs(uint-uint). abs() works somehow on linux - on windows not
    this->_progress = 0;
    emit processProgress(static_cast<int>(this->_progress));

    emit processStatus("setting up fragments");
    std::cerr << fragListFile.toStdString() << std::endl << std::flush;
    SIMPLEgenomeFragmentCollection fragmentCollection(fragListFile, 10000);
    rval = fragmentCollection.setup();

    // interrupt if something is not ok
    if (!rval) { return(false); }

    // get the number of lines
    emit processStatus("reading number of lines");
    numLines = this->getNumberLinesFromTxt(preMatrixFile);

    // interrupt if cancelled
    if (this->_isCancelled) { return(false); }

    // read the merged matrix, map the positions, and write the new matrix file
    emit processStatus("processing read pairs");

    // open files
    // infile related
    QFile infile(preMatrixFile);
    if (!infile.open(QIODevice::ReadOnly)) { return(false); }
    QTextStream in(&infile);
    in.setCodec("UTF-8");
    // outfile related
    QFile outfile(postMatrixFile);
    if (!outfile.open(QFile::WriteOnly | QFile::Text)) { return(false); }
    QTextStream out(&outfile);
    out.setCodec("UTF-8");

    // process
    QHash<QPair<uint,uint>, uint> redMat; // for reduced matrix
    QPair<uint,uint> curPair; // for reduced matrix
    QString line;
    QStringList fields;
    matrixLine curLine;
    if (!useFilter) {
        while ( !in.atEnd() ) {
            line = in.readLine();
            ++numRead;
            if ((numRead % 1000000) == 0) {
                this->_progress += (100000000.0) / numLines;
                emit processProgress(static_cast<int>(this->_progress));
                if (this->_isCancelled) {
                    rval = false;
                    break;
                }
            }
            fields = line.split("\t");
            curLine.init(fields.at(0), fields.at(1), fields.at(2), fields.at(3), fields.at(4), fields.at(5));
            fragmentCollection.mapPosition(curLine._chromA, curLine._posA, curLine._fragA);
            fragmentCollection.mapPosition(curLine._chromB, curLine._posB, curLine._fragB);
            out << curLine._chromA << "\t" <<
                   curLine._strandA << "\t" <<
                   curLine._posA << "\t" <<
                   curLine._fragA << "\t" <<
                   curLine._chromB << "\t" <<
                   curLine._strandB << "\t" <<
                   curLine._posB << "\t" <<
                   curLine._fragB << "\n";
            if (doReducedMatrix) {
                if (curLine._fragA < curLine._fragB) { curPair = qMakePair(curLine._fragA, curLine._fragB); }
                else { curPair = qMakePair(curLine._fragB, curLine._fragA); }
                if (redMat.count(curPair) == 0) { redMat[curPair] = 1;}
                else { redMat[curPair] += 1;}
            }
        }
    } else {
        while ( !in.atEnd() ) {
            line = in.readLine();
            ++numRead;
            if ((numRead % 1000000) == 0) {
                this->_progress += (100000000.0) / numLines;
                emit processProgress(static_cast<int>(this->_progress));
                if (this->_isCancelled) {
                    rval = false;
                    break;
                }
            }
            fields = line.split("\t");
            curLine.init(fields.at(0), fields.at(1), fields.at(2), fields.at(3), fields.at(4), fields.at(5));
            fragmentCollection.mapPosition(curLine._chromA, curLine._posA, curLine._fragA);
            fragmentCollection.mapPosition(curLine._chromB, curLine._posB, curLine._fragB);
            if (curLine._chromA == curLine._chromB) {
                if (curLine._strandA != curLine._strandB) {
                    readDist = curLine._posA > curLine._posB ? curLine._posA-curLine._posB : curLine._posB-curLine._posA;
                    // check inward
                    if (((curLine._posA < curLine._posB) && (curLine._strandA == "+")) || ((curLine._posB < curLine._posA) && (curLine._strandB == "+"))) {
                        if (readDist < inwardClose) { continue; }
                    }
                    // check outward
                    if (readDist < outwardClose) { continue; }
                    else { if (filterSelfCircle && (curLine._fragA == curLine._fragB)) { continue; } } //! NOTE THAT SELF CIRCLES SHOULD BE OUTWARD - IT IS DISABLED NOW
                }
            }
            // write the rest
            out << curLine._chromA << "\t" <<
                   curLine._strandA << "\t" <<
                   curLine._posA << "\t" <<
                   curLine._fragA << "\t" <<
                   curLine._chromB << "\t" <<
                   curLine._strandB << "\t" <<
                   curLine._posB << "\t" <<
                   curLine._fragB << "\n";
            if (doReducedMatrix) {
                if (curLine._fragA < curLine._fragB) { curPair = qMakePair(curLine._fragA, curLine._fragB); }
                else { curPair = qMakePair(curLine._fragB, curLine._fragA); }
                if (redMat.count(curPair) == 0) { redMat[curPair] = 1;}
                else { redMat[curPair] += 1;}
            }
        }
    }
    // close files
    // infile related
    infile.close();
    if (infile.error() != QFile::NoError) {
        std::cerr << "Error: Cannot read file " << qPrintable(preMatrixFile)
                  << ": " << qPrintable(infile.errorString())
                  << std::endl;
        rval = false;
    }
    // outfile related
    out.flush();
    outfile.close();
    if (outfile.error() != QFile::NoError) {
        std::cerr << "Error: Cannot write file " << qPrintable(postMatrixFile)
                  << ": " << qPrintable(outfile.errorString())
                  << std::endl;
        rval = false;
    }

    // interrupt if cancelled
    if (this->_isCancelled) { return(false); }

    // make the reduced matrix if necessary
    if (doReducedMatrix) {
        emit processStatus("writing reduced matrix");
        QFile reducedOutfile(reducedMatrixFile);
        if (!reducedOutfile.open(QFile::WriteOnly | QFile::Text)) { return(false); }
        out.setDevice(&reducedOutfile);
        out.setCodec("UTF-8");
        for (QHash<QPair<uint,uint>, uint>::Iterator iter = redMat.begin(); iter != redMat.end(); ++iter) {
            out << iter.key().first << "\t" <<
                   iter.key().second << "\t" <<
                   iter.value() << "\n";
        }
        out.flush();
        reducedOutfile.close();
        if (reducedOutfile.error() != QFile::NoError) {
            std::cerr << "Error: Cannot write file " << qPrintable(reducedMatrixFile)
                      << ": " << qPrintable(reducedOutfile.errorString())
                      << std::endl;
            rval = false;
        }
    }

    this->_progress = 0;
    emit processProgress(static_cast<int>(this->_progress));
    emit processStatus("idle");

    // check if it was cancelled
    if (this->_isCancelled) { rval = false; }

    return(rval);
}

bool workProcessor::CREATEFRAGMENTSdoProcess(QString workString) {
    QStringList temp = workString.split("|");
    QString restrictionPatternOrWinSize = temp.at(0);
    QString fastaFile = temp.at(1);
    QString fragListFile = temp.at(2);
    bool rval = false;
    this->_progress = 0;
    emit processProgress(static_cast<int>(this->_progress));

    emit processStatus("setting up fragments");
    SIMPLEgenomeFragmentCollection fragmentCollection(fastaFile, restrictionPatternOrWinSize, 10000);
    rval = fragmentCollection.setup();

    emit processStatus("writing fragments code");
    this->_progress = 50;
    emit processProgress(static_cast<int>(this->_progress));
    rval = fragmentCollection.writeFragments(fragListFile);

    this->_progress = 0;
    emit processProgress(static_cast<int>(this->_progress));
    emit processStatus("idle");

    // check if it was cancelled
    if (this->_isCancelled) { rval = false; }

    return(rval);
}


bool workProcessor::CREATERCODEdoProcess(QString workString) {
    QStringList temp = workString.split("|");
    QString restrictionPatternOrWinSize = temp.at(0);
    QString fastaFile = temp.at(1);
    QString rFile = temp.at(2);
    bool rval = false;
    this->_progress = 0;
    emit processProgress(static_cast<int>(this->_progress));

    emit processStatus("setting up fragments");
    SIMPLEgenomeFragmentCollection fragmentCollection(fastaFile, restrictionPatternOrWinSize, 10000);
    rval = fragmentCollection.setup();

    emit processStatus("writing R code");
    this->_progress = 50;
    emit processProgress(static_cast<int>(this->_progress));
    rval = fragmentCollection.writeOrganismFileForR(rFile);

    this->_progress = 0;
    emit processProgress(static_cast<int>(this->_progress));
    emit processStatus("idle");

    // check if it was cancelled
    if (this->_isCancelled) { rval = false; }

    return(rval);
}

//! functions previously implemented in the addtracksprocessor class

bool workProcessor::ADDTRACKSdoProcess(QString workString)
{
    QStringList temp = workString.split("|");
    bool rval = false;
    QString command = temp.at(0);

    if (command == "LOADFRAGMENTS") {
        QString fileName = temp.at(1);
        emit processStatus("loading fragments");
        rval = this->ADDTRACKSloadFragments(fileName);
        emit ADDTRACKSupdateFragmentsView(ADDTRACKSgenomeFragments.getHeader());
    }
    if (command == "SAVEFRAGMENTS") {
        QString fileName = temp.at(1);
        emit processStatus("saving fragments");
        rval = this->ADDTRACKSsaveFragments(fileName);
    }
    if (command == "CLOSEFRAGMENTS") {
        rval = this->ADDTRACKScloseFragments();
        emit ADDTRACKSupdateFragmentsView("");
    }
    if (command == "PROCESSTRACKS") {
        //! std::cerr << QTime::currentTime().toString().toStdString() << std::endl << std::flush;
        QStringList::Iterator iter;
        QString curFile;
        QString messageToSend;
        QStringList genomeAnnoFiles = temp.at(1).split(";", QString::SkipEmptyParts);
        QStringList countFeatureFiles = temp.at(2).split(";", QString::SkipEmptyParts);
        QStringList densityFeatureFiles = temp.at(3).split(";", QString::SkipEmptyParts);
        QStringList DNAmethylationFiles = temp.at(4).split(";", QString::SkipEmptyParts);
        float numTracks = genomeAnnoFiles.count() + countFeatureFiles.count() + densityFeatureFiles.count() + DNAmethylationFiles.count();
        uint numProcessed = 0;
        this->_progress = 0;
        emit processProgress(static_cast<int>(this->_progress));
        if (genomeAnnoFiles.count() > 0) {
            emit processStatus("adding genome annotation tracks");
            for (iter = genomeAnnoFiles.begin(); iter != genomeAnnoFiles.end(); ++iter ) {
                curFile = (*iter).split(QDir::separator()).last();
                messageToSend = "processing " + curFile;
                emit processStatus(messageToSend);
                rval = this->ADDTRACKSaddGenomeAnnotationTrack((*iter));
                //emit ADDTRACKSgenomeAnnotationTrackAdded((*iter));
                ++numProcessed;
                this->_progress = (numProcessed*100 / numTracks);
                emit processProgress(static_cast<int>(this->_progress));
            }
        }
        if (countFeatureFiles.count() > 0) {
            emit processStatus("adding count feature tracks");
            for (iter = countFeatureFiles.begin(); iter != countFeatureFiles.end(); ++iter ) {
                curFile = (*iter).split(QDir::separator()).last();
                messageToSend = "processing " + curFile;
                emit processStatus(messageToSend);
                rval = this->ADDTRACKSaddCountFeatureTrack((*iter));
                //emit ADDTRACKScountFeatureTrackAdded((*iter));
                ++numProcessed;
                this->_progress = (numProcessed*100 / numTracks);
                emit processProgress(static_cast<int>(this->_progress));
            }
        }
        if (densityFeatureFiles.count() > 0) {
            emit processStatus("adding density feature tracks");
            for (iter = densityFeatureFiles.begin(); iter != densityFeatureFiles.end(); ++iter ) {
                curFile = (*iter).split(QDir::separator()).last();
                messageToSend = "processing " + curFile;
                emit processStatus(messageToSend);
                rval = this->ADDTRACKSaddDensityFeatureTrack((*iter));
                //emit ADDTRACKSdensityFeatureTrackAdded((*iter));
                ++numProcessed;
                this->_progress = (numProcessed*100 / numTracks);
                emit processProgress(static_cast<int>(this->_progress));
            }
        }
        if (DNAmethylationFiles.count() > 0) {
            emit processStatus("adding DNA methylation tracks");
            for (iter = DNAmethylationFiles.begin(); iter != DNAmethylationFiles.end(); ++iter ) {
                curFile = (*iter).split(QDir::separator()).last();
                messageToSend = "processing " + curFile;
                emit processStatus(messageToSend);
                rval = this->ADDTRACKSaddDNAmethylationTrack((*iter));
                //emit ADDTRACKSDNAmethylationTrackAdded((*iter));
                ++numProcessed;
                this->_progress = (numProcessed*100 / numTracks);
                emit processProgress(static_cast<int>(this->_progress));
            }
        }
        //! std::cerr << QTime::currentTime().toString().toStdString() << std::endl << std::flush;
    }

    /* some templates
      if (this->_isCancelled) { return(false); }
        if ((numRead % 1000000) == 0) {
            this->_progress += (100000000.0) / numLines;
            emit processProgress(static_cast<int>(this->_progress));
            if (this->_isCancelled) {
                rval = false;
                break;
            }
        }
    emit processStatus("finished processing something");
    */
    return(rval);
}

bool workProcessor::ADDTRACKSloadFragments(QString fileName)
{
    bool rval = false;
    rval = this->ADDTRACKSgenomeFragments.setup(fileName);
    return(rval);
}

bool workProcessor::ADDTRACKSsaveFragments(QString fileName)
{
    bool rval = false;
    rval = this->ADDTRACKSgenomeFragments.writeFragments(fileName);
    return(rval);
}

bool workProcessor::ADDTRACKScloseFragments()
{
    this->ADDTRACKSgenomeFragments.deleteAll();
    return(true);
}

bool workProcessor::ADDTRACKSaddGenomeAnnotationTrack(QString fileName)
{
    bool rval = false;
    if (this->_isCancelled) { return(false); }
    QString fileType;
    QMap<QString, QList<genomeAnnotationElement> > featureMap;
    QList<genomeAnnotationElement> emptyList;
    genomeAnnotationElement curElement;
    if ( fileName.contains("gtf", Qt::CaseInsensitive) ) { fileType = "GTF"; }
    else if ( fileName.contains("gff", Qt::CaseInsensitive) ) { fileType = "GFF3"; }
    else { return(false); }
    annotationReader reader(fileName, fileType);
    while ( reader.readLine() ) {
        if (featureMap.count(reader.feature) == 0) {
            featureMap[reader.feature] = emptyList;
            featureMap[reader.feature].reserve(100000);
        }
        curElement.init(reader.chrom, reader.start, reader.end);
        featureMap[reader.feature] << curElement;
    }
    if (this->_isCancelled) { return(false); }
    QString curTrackName;
    QMap<QString, QList<genomeAnnotationElement> >::Iterator iter;
    for (iter = featureMap.begin(); iter != featureMap.end(); ++iter) {
        curTrackName = "ann_" + iter.key();
        if (this->_isCancelled) { return(false); }
        this->ADDTRACKSgenomeFragments.addGenomeAnnotationTrack(curTrackName, iter.value());
        emit ADDTRACKSupdateFragmentsView(ADDTRACKSgenomeFragments.getHeader());
    }
    emit ADDTRACKSgenomeAnnotationTrackAdded(fileName);
    return(rval);
}

bool workProcessor::ADDTRACKSaddCountFeatureTrack(QString fileName)
{
    bool rval = false;
    if (this->_isCancelled) { return(false); }
    QString curTrackName = "sum_" + fileName.split(QDir::separator()).last().split(".").first();
    seqan::BamStream reader(fileName.toStdString().c_str());
    std::string curChromStd;
    uint curLen;
    QVector<float> curVec;
    QMap<QString, QVector<float> > staPoCov;
    QMap<int,QString> IDtoNAME;
    QString curChrom;
    for (uint i = 0; i < length(reader.header.sequenceInfos); ++i) {
        curChromStd = seqan::toCString(reader.header.sequenceInfos[i].i1);
        curLen = reader.header.sequenceInfos[i].i2;
        curVec.clear();
        curVec.fill(0, curLen);
        curChrom = QString::fromStdString(curChromStd);
        staPoCov[curChrom] = curVec;
        IDtoNAME[i] = curChrom;
    }
    seqan::BamAlignmentRecord record;
    uint numLines = 0;
    float w;
    uint tagID;
    while (!atEnd(reader)) {
        ++numLines;
        if ((numLines % 1000000) == 0) { if (this->_isCancelled) { break; }}
        seqan::readRecord(record, reader);
        if (seqan::hasFlagUnmapped(record)) { continue; }
        seqan::BamTagsDict tags(record.tags);
        if (seqan::findTagKey(tagID, tags, "XW")) {seqan::extractTagValue(w, tags, tagID);}
        else { w = 1; }
        curChrom = IDtoNAME.value(static_cast<int>(record.rID));
        staPoCov[curChrom][static_cast<int>(record.beginPos)] += w;
    }
    seqan::close(reader);
    if (this->_isCancelled) { return(false); }
    this->ADDTRACKSgenomeFragments.addCountFeatureTrack(curTrackName, staPoCov);
    emit ADDTRACKSupdateFragmentsView(ADDTRACKSgenomeFragments.getHeader());
    emit ADDTRACKScountFeatureTrackAdded(fileName);
    return(rval);
}

bool workProcessor::ADDTRACKSaddDensityFeatureTrack(QString fileName)
{
    bool rval = false;
    if (this->_isCancelled) { return(false); }
    QString curTrackName = "den_" + fileName.split(QDir::separator()).last().split(".").first();
    seqan::BamStream reader(fileName.toStdString().c_str());
    std::string curChromStd;
    uint curLen;
    QVector<float> curVec;
    QMap<QString, QVector<float> > cov;
    QMap<int,QString> IDtoNAME;
    QString curChrom;
    uint i;
    uint k;
    for (i = 0; i < length(reader.header.sequenceInfos); ++i) {
        curChromStd = seqan::toCString(reader.header.sequenceInfos[i].i1);
        curLen = reader.header.sequenceInfos[i].i2;
        curVec.clear();
        curVec.fill(0, curLen);
        curChrom = QString::fromStdString(curChromStd);
        cov[curChrom] = curVec;
        IDtoNAME[i] = curChrom;
    }
    seqan::BamAlignmentRecord record;
    uint numLines = 0;
    float w;
    uint tagID;
    uint start;
    uint end;
    while (!atEnd(reader)) {
        ++numLines;
        if ((numLines % 1000000) == 0) { if (this->_isCancelled) { break; }}
        seqan::readRecord(record, reader);
        if (seqan::hasFlagUnmapped(record)) { continue; }
        seqan::BamTagsDict tags(record.tags);
        if (seqan::findTagKey(tagID, tags, "XW")) {seqan::extractTagValue(w, tags, tagID);}
        else { w = 1; }
        curChrom = IDtoNAME.value(static_cast<int>(record.rID));
        start = static_cast<uint>(record.beginPos);
        for (i = 0; i < seqan::length(record.cigar); ++i) {
            end = start + record.cigar[i].count - 1; // minus one because be need closed interval (as in the annotation)
            if (record.cigar[i].operation == 'M') {
                for (k = start; k <= end; ++k) { cov[curChrom][k] += w; }
            }
            start = end + 1; // we need the plus one again, as the start needs to be included
        }
    }
    seqan::close(reader);
    if (this->_isCancelled) { return(false); }
    uint minCov = 1;
    this->ADDTRACKSgenomeFragments.addDensityFeatureTrack(curTrackName, cov, minCov);
    emit ADDTRACKSupdateFragmentsView(ADDTRACKSgenomeFragments.getHeader());
    emit ADDTRACKSdensityFeatureTrackAdded(fileName);
    return(rval);
}

bool workProcessor::ADDTRACKSaddDNAmethylationTrack(QString fileName)
{
    if (this->_isCancelled) { return(false); }
    QString curTrackName = "den_" + fileName.split(QDir::separator()).last().split(".").first();
    QMap<QString, uint> chromSizes = this->ADDTRACKSgenomeFragments.getChromosomeSizes();
    if (chromSizes.count() == 0) {return(false);}
    QMap<QString, QVector<QString> > methState;
    QVector<QString> curVec;
    for (QMap<QString, uint>::Iterator iter = chromSizes.begin(); iter != chromSizes.end(); ++iter) {
        curVec.clear();
        curVec.fill("n", iter.value());
        methState[iter.key()] = curVec;
    }
    if (this->_isCancelled) { return(false); }
    // read the file
    QFile file(fileName);
    QTextStream in;
    // check if file available (interrupt if not)
    if (!file.open(QIODevice::ReadOnly)) {return(false);}
    // set some stuff
    in.setDevice(&file);
    in.setCodec("UTF-8");
    // read the file
    QString line;
    QStringList fields;
    QString chrom;
    uint pos;
    QString state;
    // the entries
    while ( !in.atEnd() ) {
        line = in.readLine();
        fields = line.split("\t");
        chrom = fields.at(0);
        pos = fields.at(1).toUInt();
        state = fields.at(2);
        if (methState.count(chrom) == 0) {
            this->_errorMessage = "DNA methylation track had a chromosome ID that was not in the fragments";
            return(false);
        }
        methState[chrom][pos] = state;
    }
    // close and check for errors
    file.close();
    if (file.error() != QFile::NoError) {return(false);}
    if (this->_isCancelled) { return(false); }
    this->ADDTRACKSgenomeFragments.addDNAmethylationTrack(curTrackName, methState);
    emit ADDTRACKSupdateFragmentsView(ADDTRACKSgenomeFragments.getHeader());
    emit ADDTRACKSDNAmethylationTrackAdded(fileName);
    return(true);
}




