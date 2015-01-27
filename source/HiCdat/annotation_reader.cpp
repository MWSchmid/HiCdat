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

#include <QtGui>

#ifndef ANNOTATION_READER_H
#include "annotation_reader.h" //! do not include while compiling (as it is inline)!
#endif



inline annotationReader::annotationReader(QString& annofile, QString& annotype) {
    filetype = annotype;
    file.setFileName(annofile);
    if(!file.open(QIODevice::ReadOnly)) {
        std::cerr << "Can not open file for reading: "
                  << qPrintable(file.errorString()) << std::endl;
        return;
    }
    in.setDevice(&file);
    in.setCodec("UTF-8");
    linenumber = 0;
    //initialize also others variables?
}

inline annotationReader::~annotationReader() {
    file.close();
}

inline bool annotationReader::readLine() {
    bool rval = false;
    if (filetype == "BED") { rval = readBedLine(); }
    else if (filetype == "GTF") { rval = readGtfLine(); }
    else if (filetype == "GFF3") { rval = readGffLine(); }
    else { std::cerr << "unknown file type" << std::endl; }
    ++linenumber;
    return(rval);
}

//! for BED reading
inline bool annotationReader::readBedLine() {
    bool rval = false;
    if(!in.atEnd()) {
        QString line = in.readLine();
        QStringList fields = line.split('\t');
        if (fields.size() >= 6) {
            chrom = fields.takeFirst();
            start = fields.takeFirst().toUInt();
            end = fields.takeFirst().toUInt();
            name = fields.takeFirst();
            QString tempscore = fields.takeFirst();
            if (tempscore == ".") { score = 0; }
            else { score = tempscore.toFloat(); }
            strand = fields.takeFirst();
            source = "BED";
            feature = "none";
            phase = ".";
            parent = "none";
            rval = true;
        } else {
            if (line.length() == 0) {
                std::cerr << "skipping empty line" << std::endl << std::flush;
                rval = this->readBedLine();
            }
            else { if (line.at(0) == '#') {
                    std::cerr << "skipping comment line:" << std::endl << line.toStdString() << std::endl << std::flush;
                    rval = this->readBedLine(); }
            }
        }
    }
    return(rval);
}

//! for GTF reading
inline bool annotationReader::getGtfAttributes(QString& str)
{
    bool rval = false;
    Attr.clear(); //very important
    QRegExp rx("(\\w+|\\w+\\.\\w+|\\w+\\.\\w+\\.\\w+|\\w+\\-\\w+|\\w+\\(\\w+\\)\\w+)");
    QString key = "none";
    int counter = 1;
    int pos = 0;
    while ((pos = rx.indexIn(str, pos)) != -1) {
        if ((counter % 2) != 0) { key = rx.cap(1); }
        else { Attr.insert(key, rx.cap(1)); }
        pos += rx.matchedLength();
        ++counter;
    }
    rval = (counter != 1);
    if (Attr.contains("gene_id") && Attr.contains("transcript_id")) {
        if (Attr["transcript_id"] != "none") {
            if ((feature == "exon") || (feature == "CDS")) {
                name = "none";
                parent = Attr["transcript_id"];
            } else {
                name = Attr["transcript_id"];
                parent = Attr["gene_id"];
            }
        } else {
            name = Attr["gene_id"];
            parent = "none";
        }
    } else {
        name = "none";
        parent = "none";
    }
    return(rval);
}

inline bool annotationReader::readGtfLine() {
    bool rval = false;
    if(!in.atEnd()) {
        QString line = in.readLine();
        QStringList fields = line.split('\t');
        if (fields.length() == 9) {
            chrom = fields.takeFirst();
            source = fields.takeFirst();
            feature = fields.takeFirst();
            start = fields.takeFirst().toUInt();
            end = fields.takeFirst().toUInt();
            QString tempscore = fields.takeFirst();
            if (tempscore == ".") { score = 0; }
            else { score = tempscore.toFloat(); }
            strand = fields.takeFirst();
            phase = fields.takeFirst();
            QString attributes = fields.takeFirst();
            if (getGtfAttributes(attributes)) { rval = true; }
        } else {
            if (line.length() == 0) {
                std::cerr << "skipping empty line" << std::endl << std::flush;
                rval = this->readGtfLine();
            }
            else { if (line.at(0) == '#') {
                    std::cerr << "skipping comment line:" << std::endl << line.toStdString() << std::endl << std::flush;
                    rval = this->readGtfLine(); }
            }
        }
    }
    return(rval);
}

//! for GFF3 reading
inline bool annotationReader::getGffAttributes(QString& str) {
    bool rval = false;
    Attr.clear(); //very important
    QStringList pairs = str.split(';', QString::SkipEmptyParts);
    foreach (QString pair, pairs) {
        QStringList keyval = pair.split("=", QString::SkipEmptyParts);
        Attr.insert(keyval.at(0), keyval.at(1));
    }
    rval = (pairs.length() > 0);
    if (rval) {
        bool has_ID = Attr.contains("ID");
        bool has_Name = Attr.contains("Name");
        bool has_Parent = Attr.contains("Parent");
        if (has_ID && has_Name) {
            name = Attr["ID"];
            if (Attr["ID"] == Attr["Name"]) {
                if (has_Parent) { parent = Attr["Parent"]; }
                else { parent = "none"; }
            } else {
                parent = Attr["Name"];
            }
        } else {
            name = "none";
            if (has_Parent) {
                int pos = Attr["Parent"].indexOf(',');
                if (pos != -1) {
                    parent = Attr["Parent"].split(',').at(1);
                } else {
                    parent = Attr["Parent"];
                }
            } else {
                parent = "none";
            }
        }
    }
    return(rval);
}

inline bool annotationReader::readGffLine() {
    bool rval = false;
    if(!in.atEnd()) {
        QString line = in.readLine();
        QStringList fields = line.split('\t');
        if (fields.length() == 9) {
            chrom = fields.takeFirst();
            source = fields.takeFirst();
            feature = fields.takeFirst();
            start = fields.takeFirst().toUInt();
            end = fields.takeFirst().toUInt();
            QString tempscore = fields.takeFirst();
            if (tempscore == ".") { score = 0; }
            else { score = tempscore.toFloat(); }
            strand = fields.takeFirst();
            phase = fields.takeFirst();
            QString attributes = fields.takeFirst();
            if (getGffAttributes(attributes)) { rval = true; }
        } else {
            if (line.length() == 0) {
                std::cerr << "skipping empty line" << std::endl << std::flush;
                rval = this->readGffLine();
            }
            else { if (line.at(0) == '#') {
                    std::cerr << "skipping comment line:" << std::endl << line.toStdString() << std::endl << std::flush;
                    rval = this->readGffLine(); }
            }
        }
    }
    return(rval);
}

