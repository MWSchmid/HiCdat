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

#ifndef GENOMEFRAGMENTCOLLECTION_H
#define GENOMEFRAGMENTCOLLECTION_H

#include <QtCore>
#include "genomefragment.h"
#include "genomeannotationelement.h"

class genomeFragmentCollection
{
private:
    //! the rootItem of the database and the index
    QList<genomeFragment*> genomeFragments;
    QStringList header;
    QStringList onlyTracksHeader;
    QHash< QPair<QString,uint>, QVector<genomeFragment*> > CPindex;
    uint stepSize;

private:
    uint getOffset(const uint& position);
    void createIndex();
    bool loadFragments(QString fileName);
    void prepareForNewTrack (const QString& trackName);

    //! FOR GENOMEANNOFEATURES - map longer regions - ergo also regions spanning a whole fragment - here there might also be multiple mappings by the way
    void mapLongRange(const QString &chrom, const uint &start, const uint &end, QList<genomeFragment*> &results);
    //! FOR COUNTFEATURES - map short regions and add a simple count - careful - there might be no mapping
    //! well - as we map only positions, for the short count features, one can also do it with a coverage vector
    bool mapRange(const QString &chrom, const uint &start, const uint &end, genomeFragment* result);
    bool mapPosition(const QString &chrom, const uint &pos, genomeFragment* result);
    //! FOR DENSITYFEATURES - supply vector with coverage (no mapping)
    //! FOR DNAMETHYLATION - supply vector with 'm', 'n', 'u' for methylated, noCinContext, unmethylated (no mapping)

public:
    genomeFragmentCollection();
    ~genomeFragmentCollection();
    void deleteAll();
    QString getHeader();

    bool setup(QString fileName);
    bool writeFragments(QString fragListFile);

    QMap<QString, uint> getChromosomeSizes(); //! this is approximate - it takes the highes end value of a fragment

    //! add tracks
    //! genome annotations - give a vector with all elements of a certain feature (GTF/GFF allows multiple features)
    void addGenomeAnnotationTrack(const QString& trackName, const QList<genomeAnnotationElement> &genomeAnnotationElements);
    //! short count features - as we map only positions, for the short count features, one can also do it with a coverage map (chrom to vector) (mark only the start position!)
    void addCountFeatureTrack(const QString& trackName, const QMap<QString, QVector<float> >& startPosCoverage);
    //! density features - supply a coverage map (chrom to vector)
    void addDensityFeatureTrack(const QString& trackName, const QMap<QString, QVector<float> >& coverage, const uint &minCov);
    //! DNA methylation - supply a state map (chrom to vector that marks all Cs as meth or unmeth)
    void addDNAmethylationTrack(const QString& trackName, const QMap<QString, QVector<QString> >& methState);

};

#endif // GENOMEFRAGMENTCOLLECTION_H
