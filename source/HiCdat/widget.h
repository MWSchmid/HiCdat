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

#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>
#include "workprocessor.h"

namespace Ui {
class Widget;
}

class Widget : public QWidget
{
    Q_OBJECT
    
public:
    explicit Widget(QWidget *parent = 0);
    ~Widget();
    
private:
    Ui::Widget *ui;
    workProcessor worker;
    int progress;
    QString status;
    QString ADDTRACKSheaderAtBeginning;
    QString ADDTRACKScurrentHeader;

private slots:
    void removeFromQueue(QString workString);
    void updateProgress(int progressIN);
    void updateStatus(QString statusIN);
    //! from merge reads  >> merge tab
    void MERGEsearchReadsInfileR1();
    void MERGEsearchReadsInfileR2();
    void MERGEsearchMatrixOutfile();
    void MERGEsendProcess();
    //! from map >> create fragments tab
    void CREATEFRAGMENTSsearchFastaFile();
    void CREATEFRAGMENTSsearchFragListFile();
    void CREATEFRAGMENTSsendProcess();
    //! from map >> map tab
    void MAPREADSsearchFragListFile();
    void MAPREADSsearchPreMatrixFile();
    void MAPREADSsearchPostMatrixFile();
    void MAPREADSsearchReducedMatrixFile();
    void MAPREADSsendProcess();
    //! from tracks >> add tracks tab
    void ADDTRACKSsetFragmentView(QString header); // used to update the fragment listWidget - the string should contain all headers
    void ADDTRACKSremoveGenomeAnnotationFiles(QString fileName); // used to remove the track in in the "toProcess" side
    void ADDTRACKSremoveCountFeatureFiles(QString fileName); // used to remove the track in in the "toProcess" side
    void ADDTRACKSremoveDensityFeatureFiles(QString fileName); // used to remove the track in in the "toProcess" side
    void ADDTRACKSremoveDNAmethylationFiles(QString fileName); // used to remove the track in in the "toProcess" side
    void ADDTRACKSloadFragments();
    void ADDTRACKSsaveFragments();
    void ADDTRACKScloseFragments();
    void ADDTRACKSaddGenomeAnnotationFiles();
    void ADDTRACKSaddCountFeatureFiles();
    void ADDTRACKSaddDensityFeatureFiles();
    void ADDTRACKSaddDNAmethylationFiles();
    void ADDTRACKSsendProcess();
    //! from map >> create R-code tab
    void CREATERCODEsearchFastaFile();
    void CREATERCODEsearchRfile();
    void CREATERCODEsendProcess();



signals:
    void analysisCanceled();

};

#endif // WIDGET_H
