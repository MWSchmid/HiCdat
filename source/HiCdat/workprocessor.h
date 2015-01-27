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

#ifndef WORKPROCESSOR_H
#define WORKPROCESSOR_H

#include <QtGui>
#ifdef Q_WS_MAC
#include <QtWidgets> // Qt5 instead of 4
#endif
#include <QMap>
#include "helperClasses.h"

#include "genomefragmentcollection.h"


class workProcessor: public QThread
{ Q_OBJECT
private:
    //! some stuff
    bool _allFine;
    QString _errorMessage;
    float _progress;
    bool _isCancelled;

    //! variables used for the thread control
    QQueue<QString> _work; // a queue with the project files that shall be processed
    QWaitCondition _workAdded;
    QMutex _mutex;

    //! the genomefragmentcollection for the addTracks
    genomeFragmentCollection ADDTRACKSgenomeFragments;

public:
    workProcessor();
    ~workProcessor();
    void addWork(QString workString); // will add some work to the queue

protected:
    void run();

private:
    // react to added work
    bool doProcess(QString workString);
    // non-specific gunctions
    quint64 getNumberLinesFromTxt(QString fileName);
    quint64 getNumberLinesFromBam(QString fileName);
    QMap<uint, QString> getIDtoNAMEfromBam(QString fileName);
    //! functions previously implemented in the readprocessor class
    bool MERGEPAIRSdoProcess(QString workString);
    void MERGEPAIRSaddPartner(QString fileName, QMap<std::string, readPair> & readCollection, QMap<uint, QString> &idToName, quint64 numLines);
    void MERGEPAIRSwriteCollection(QTextStream& out, QMap<std::string, readPair> &readCollection);
    //! functions previously implemented in the mapprocessor class
    bool MAPREADPAIRSdoProcess(QString workString);
    bool CREATEFRAGMENTSdoProcess(QString workString);
    bool CREATERCODEdoProcess(QString workString);
    //! functions previously implemented in the addtracksprocessor class
    bool ADDTRACKSdoProcess(QString workString);
    bool ADDTRACKSloadFragments(QString fileName);
    bool ADDTRACKSsaveFragments(QString fileName);
    bool ADDTRACKScloseFragments();
    bool ADDTRACKSaddGenomeAnnotationTrack(QString fileName);
    bool ADDTRACKSaddCountFeatureTrack(QString fileName);
    bool ADDTRACKSaddDensityFeatureTrack(QString fileName);
    bool ADDTRACKSaddDNAmethylationTrack(QString fileName);

public slots:
    void cancelProcessing();

signals:
    void errorMessage(QString message);
    void processStatus(QString status); // signals what is going on
    void processProgress(int progress);
    void idleAgain();
    void workFinished(QString workString);
    //! addTrack specific stuff
    void ADDTRACKSupdateFragmentsView(QString header); // used to update the fragment listWidget - the string should contain all headers
    void ADDTRACKSgenomeAnnotationTrackAdded(QString fileName); // used to remove the track in in the "toProcess" side
    void ADDTRACKScountFeatureTrackAdded(QString fileName); // used to remove the track in in the "toProcess" side
    void ADDTRACKSdensityFeatureTrackAdded(QString fileName); // used to remove the track in in the "toProcess" side
    void ADDTRACKSDNAmethylationTrackAdded(QString fileName); // used to remove the track in in the "toProcess" side
};

#endif // WORKPROCESSOR_H
