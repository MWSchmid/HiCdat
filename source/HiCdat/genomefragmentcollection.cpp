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

#include <iostream>
#include "genomefragmentcollection.h"

genomeFragmentCollection::genomeFragmentCollection() : stepSize(10000) {}

genomeFragmentCollection::~genomeFragmentCollection() {
    qDeleteAll(genomeFragments);
}

void genomeFragmentCollection::deleteAll() {
    this->header.clear();
    CPindex.clear();
    qDeleteAll(genomeFragments);
    genomeFragments.clear();
}

QString genomeFragmentCollection::getHeader() {
    return (this->header.join("|"));
}

void genomeFragmentCollection::prepareForNewTrack (const QString& trackName) {
    if (!this->header.contains(trackName)) {
        this->header << trackName;
        foreach (genomeFragment* item, genomeFragments) { item->_data << 0; }
    }
}

bool genomeFragmentCollection::setup(QString fileName) {
    bool rval = true;
    rval = this->loadFragments(fileName);
    if (rval) { this->createIndex(); }
    return(rval);
}

uint genomeFragmentCollection::getOffset(const uint& position) {
    uint x = (position / this->stepSize); // uses integer division
    return(x);
}

void genomeFragmentCollection::createIndex() {
    QString chrom;
    uint start;
    uint end;
    uint startStep;
    uint endStep;
    uint stepCounter;
    QPair<QString,uint> addKey;
    QVector<genomeFragment*> emptyVec;
    emptyVec.reserve(100);
    foreach (genomeFragment* frag, this->genomeFragments) {
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

QMap<QString, uint> genomeFragmentCollection::getChromosomeSizes() {
    QMap<QString, uint> out;
    foreach (genomeFragment* item, this->genomeFragments) {
        if (out.count(item->_chrom) == 0) { out[item->_chrom] = 0; }
        if (out[item->_chrom] < item->_end) { out[item->_chrom] = item->_end; }
    }
    return(out);
}

bool genomeFragmentCollection::loadFragments(QString fileName) {
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
    uint fragNum;
    uint start;
    uint end;
    // header
    if (!in.atEnd()) { line = in.readLine(); }
    else { return(false); }
    this->header = line.split("\t");
    if (header.count() > 4) {
        for (int i = 4; i < header.count(); ++i) {
            onlyTracksHeader << header.at(i);
        }
    }
    // the entries
    while ( !in.atEnd() ) {
        line = in.readLine();
        fields = line.split("\t");
        fragNum = fields.takeFirst().toUInt();
        chrom = fields.takeFirst();
        start = fields.takeFirst().toUInt();
        end = fields.takeFirst().toUInt();
        genomeFragment *frag = new genomeFragment(chrom, start, end, fragNum);
        foreach (QString cur, fields) {
            frag->_data += cur;
        }
        this->genomeFragments += frag;
    }
    // close and check for errors
    file.close();
    if (file.error() != QFile::NoError) {return(false);}
    return(true);
}


bool genomeFragmentCollection::writeFragments(QString fragListFile) {
    QFile outfile(fragListFile);
    if (!outfile.open(QFile::WriteOnly | QFile::Text)) { return(false); }
    QTextStream out(&outfile);
    out.setCodec("UTF-8");

    out << this->header.join("\t") << "\n";
    foreach (genomeFragment* frag, this->genomeFragments) {
        out << frag->_fragNum << "\t" <<
               frag->_chrom << "\t" <<
               frag->_start << "\t" <<
               frag->_end;
        foreach (QVariant dat, frag->_data) {
            out << "\t" << dat.toString();
        }
        out << "\n";
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

void genomeFragmentCollection::mapLongRange(const QString &chrom, const uint &start, const uint &end, QList<genomeFragment*> &results) {
    uint startOffset = this->getOffset(start);
    uint endOffset = this->getOffset(end);
    QPair<QString,uint> offset;
    for (uint i = startOffset; i <= endOffset; ++i) {
        offset = qMakePair(chrom, i);
        foreach (genomeFragment* frag, this->CPindex.value(offset)) {
            if (((frag->_start <= start) && (frag->_end >= start)) || ((frag->_start <= end) && (frag->_end >= end))) {
                results << frag;
            }
        }
    }
}

bool genomeFragmentCollection::mapRange(const QString &chrom, const uint &start, const uint &end, genomeFragment* result) {
    QPair<QString,uint> offset = qMakePair(chrom, this->getOffset(start));
    foreach (genomeFragment* frag, this->CPindex.value(offset)) {
        if ((frag->_start <= start) && (frag->_end >= end)) {
            result = frag;
            return(true);
        }
    }
    return(false);
}

bool genomeFragmentCollection::mapPosition(const QString &chrom, const uint &pos, genomeFragment* result) {
    QPair<QString,uint> offset = qMakePair(chrom, this->getOffset(pos));
    foreach (genomeFragment* frag, this->CPindex.value(offset)) {
        if ((frag->_start <= pos) && (frag->_end >= pos)) {
            result = frag;
            return(true);
        }
    }
    return(false);
}


void genomeFragmentCollection::addGenomeAnnotationTrack(const QString& trackName, const QList<genomeAnnotationElement> &genomeAnnotationElements)
{
    this->prepareForNewTrack(trackName);
    float curValue;
    QList<genomeFragment*> curFragments;
    QList<genomeAnnotationElement>::const_iterator iter;
    for (iter = genomeAnnotationElements.begin(); iter != genomeAnnotationElements.end(); ++iter) {
        curFragments.clear();
        this->mapLongRange((*iter)._chrom, (*iter)._start, (*iter)._end, curFragments);
        foreach (genomeFragment* item, curFragments) {
            curValue = item->_data.last().toFloat();
            if (((*iter)._start >= item->_start) && ((*iter)._end <= item->_end)) {curValue += 1;}
            else {curValue += 0.5;}
            item->_data.last() = curValue;
        }
    }
}

void genomeFragmentCollection::addCountFeatureTrack(const QString& trackName, const QMap<QString, QVector<float> >& startPosCoverage)
{
    this->prepareForNewTrack(trackName);
    float startPosSum;
    foreach (genomeFragment* item, genomeFragments) {
        startPosSum = 0;
        if (startPosCoverage.count(item->_chrom) == 0) { std::cerr << "chromosome missing" << std::endl << std::flush; break;}
        for (uint i = item->_start; i <= item->_end; ++i) {
            startPosSum += startPosCoverage[item->_chrom].at(i);
        }
        item->_data.last() = startPosSum;
    }
}

void genomeFragmentCollection::addDensityFeatureTrack(const QString& trackName, const QMap<QString, QVector<float> >& coverage, const uint &minCov)
{
    this->prepareForNewTrack(trackName);
    uint numCov;
    uint fragSize;
    float density;
    foreach (genomeFragment* item, genomeFragments) {
        numCov = 0;
        if (coverage.count(item->_chrom) == 0) { std::cerr << "chromosome missing" << std::endl << std::flush; break;}
        fragSize = item->_end-item->_start;
        for (uint i = item->_start; i <= item->_end; ++i) {
            if (coverage[item->_chrom].at(i) >= minCov) { ++numCov;}
        }
        density = (numCov/static_cast<float>(fragSize)) * 100;
        item->_data.last() = density;
    }
}

void genomeFragmentCollection::addDNAmethylationTrack(const QString& trackName, const QMap<QString, QVector<QString> >& methState)
{
    this->prepareForNewTrack(trackName);
    uint numMeth;
    uint numTot;
    float density;
    foreach (genomeFragment* item, genomeFragments) {
        numMeth = 0;
        numTot = 0;
        if (methState.count(item->_chrom) == 0) { std::cerr << "chromosome missing" << std::endl << std::flush; break;}
        for (uint i = item->_start; i <= item->_end; ++i) {
            if (i >= methState[item->_chrom].count()) { continue; }
            if (methState[item->_chrom].at(i) == "m") { ++numMeth; ++numTot; }
            else if (methState[item->_chrom].at(i) == "u") { ++numTot; }
        }
        density = (numMeth/static_cast<float>(numTot)) * 100;
        item->_data.last() = density;
    }
}








