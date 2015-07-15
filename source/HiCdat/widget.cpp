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

#include "widget.h"
#include "ui_widget.h"

Widget::Widget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Widget)
{
    ui->setupUi(this);
    setWindowTitle(tr("HiCdatPre"));

    //! the header of the in-memory fragments
    this->ADDTRACKSheaderAtBeginning = "";
    this->ADDTRACKScurrentHeader = "";

    //! progress bar
    progress = 0;
    status = "idle";
    ui->progressLabel->setText(status);
    ui->progressBar->setValue(progress);

    //! hide the filter self circle check box
    ui->MAPREADSselfCircleCheck->hide();

    //! connections to the worker
    connect(ui->progressPushButtonCancelCurrentJob, SIGNAL(clicked()), &worker, SLOT(cancelProcessing()));
    connect(&worker, SIGNAL(processStatus(QString)), this, SLOT(updateStatus(QString)));
    connect(&worker, SIGNAL(processProgress(int)), this, SLOT(updateProgress(int)));
    connect(&worker, SIGNAL(errorMessage(QString)), this, SLOT(updateStatus(QString)));
    connect(&worker, SIGNAL(workFinished(QString)), this, SLOT(removeFromQueue(QString)));
    //! connections in merge tab
    connect(ui->MERGEr1_find, SIGNAL(clicked()), this, SLOT(MERGEsearchReadsInfileR1()));
    connect(ui->MERGEr2_find, SIGNAL(clicked()), this, SLOT(MERGEsearchReadsInfileR2()));
    connect(ui->MERGEout_find, SIGNAL(clicked()), this, SLOT(MERGEsearchMatrixOutfile()));
    connect(ui->RUNmergeFiles, SIGNAL(clicked()), this, SLOT(MERGEsendProcess()));
    //! connections in the create fragments tab
    connect(ui->CREATEFRAGMENTSfasta_find, SIGNAL(clicked()), this, SLOT(CREATEFRAGMENTSsearchFastaFile()));
    connect(ui->CREATEFRAGMENTSfragList_find, SIGNAL(clicked()), this, SLOT(CREATEFRAGMENTSsearchFragListFile()));
    connect(ui->RUNCREATEFRAGMENTScreateFragments, SIGNAL(clicked()), this, SLOT(CREATEFRAGMENTSsendProcess()));
    //! connections in the map pairs
    connect(ui->MAPREADSfragList_find, SIGNAL(clicked()), this, SLOT(MAPREADSsearchFragListFile()));
    connect(ui->MAPREADSpreMatrix_find, SIGNAL(clicked()), this, SLOT(MAPREADSsearchPreMatrixFile()));
    connect(ui->MAPREADSpostMatrix_find, SIGNAL(clicked()), this, SLOT(MAPREADSsearchPostMatrixFile()));
    connect(ui->MAPREADSreducedMatrix_find, SIGNAL(clicked()), this, SLOT(MAPREADSsearchReducedMatrixFile()));
    connect(ui->MAPREADSRUNmapReadPairs, SIGNAL(clicked()), this, SLOT(MAPREADSsendProcess()));
    //! connections in the add tracks
    connect(&worker, SIGNAL(ADDTRACKSupdateFragmentsView(QString)), this, SLOT(ADDTRACKSsetFragmentView(QString)));
    connect(&worker, SIGNAL(ADDTRACKSgenomeAnnotationTrackAdded(QString)), this, SLOT(ADDTRACKSremoveGenomeAnnotationFiles(QString)));
    connect(&worker, SIGNAL(ADDTRACKScountFeatureTrackAdded(QString)), this, SLOT(ADDTRACKSremoveCountFeatureFiles(QString)));
    connect(&worker, SIGNAL(ADDTRACKSdensityFeatureTrackAdded(QString)), this, SLOT(ADDTRACKSremoveDensityFeatureFiles(QString)));
    connect(&worker, SIGNAL(ADDTRACKSDNAmethylationTrackAdded(QString)), this, SLOT(ADDTRACKSremoveDNAmethylationFiles(QString)));
    connect(ui->RUNaddAllTracksButton, SIGNAL(clicked()), this, SLOT(ADDTRACKSsendProcess()));
    connect(ui->ADDTRACKSloadFragmentsButton, SIGNAL(clicked()), this, SLOT(ADDTRACKSloadFragments()));
    connect(ui->ADDTRACKSsaveFragmentsButton, SIGNAL(clicked()), this, SLOT(ADDTRACKSsaveFragments()));
    connect(ui->ADDTRACKScloseFragmentsButton, SIGNAL(clicked()), this, SLOT(ADDTRACKScloseFragments()));
    connect(ui->ADDTRACKSgenomeAnnotationFind, SIGNAL(clicked()), this, SLOT(ADDTRACKSaddGenomeAnnotationFiles()));
    connect(ui->ADDTRACKScountFeaturesFind, SIGNAL(clicked()), this, SLOT(ADDTRACKSaddCountFeatureFiles()));
    connect(ui->ADDTRACKSdensityFeaturesFind, SIGNAL(clicked()), this, SLOT(ADDTRACKSaddDensityFeatureFiles()));
    connect(ui->ADDTRACKSDNAmethylationFind, SIGNAL(clicked()), this, SLOT(ADDTRACKSaddDNAmethylationFiles()));
    //! connections in the create R-code tab
    connect(ui->CREATERCODEfasta_find, SIGNAL(clicked()), this, SLOT(CREATERCODEsearchFastaFile()));
    connect(ui->CREATERCODERfile_find, SIGNAL(clicked()), this, SLOT(CREATERCODEsearchRfile()));
    connect(ui->CREATERCODERUNcreateRfile, SIGNAL(clicked()), this, SLOT(CREATERCODEsendProcess()));
    //! set field colors
    QString outputfieldColor = "background-color: rgba(255, 230, 150, 150)";
    QString inputfieldColor = "background-color: rgba(200, 255, 200, 150)";
    ui->MERGEr1_lineEdit->setStyleSheet(inputfieldColor);
    ui->MERGEr2_lineEdit->setStyleSheet(inputfieldColor);
    ui->MERGEout_lineEdit->setStyleSheet(outputfieldColor);
    //ui->MERGEnumReads_spinBox->setStyleSheet(inputfieldColor);
    ui->CREATEFRAGMENTSresPat_lineEdit->setStyleSheet(inputfieldColor);
    ui->CREATEFRAGMENTSfasta_lineEdit->setStyleSheet(inputfieldColor);
    ui->CREATEFRAGMENTSfragList_lineEdit->setStyleSheet(outputfieldColor);
    ui->MAPREADSfragList_lineEdit->setStyleSheet(inputfieldColor);
    ui->MAPREADSpreMatrix_lineEdit->setStyleSheet(inputfieldColor);
    ui->MAPREADSpostMatrix_lineEdit->setStyleSheet(outputfieldColor);
    ui->MAPREADSreducedMatrix_lineEdit->setStyleSheet(outputfieldColor);
    //ui->MAPREADSinwardCloseBox->setStyleSheet(inputfieldColor);
    //ui->MAPREADSoutwardCloseBox->setStyleSheet(inputfieldColor);
    ui->CREATERCODEresPat_lineEdit->setStyleSheet(inputfieldColor);
    ui->CREATERCODEfasta_lineEdit->setStyleSheet(inputfieldColor);
    ui->CREATERCODErFile_lineEdit->setStyleSheet(outputfieldColor);
}

Widget::~Widget()
{
    worker.cancelProcessing();
    worker.addWork("STOPTHREAD");
    worker.wait();
    delete ui;
}

void Widget::removeFromQueue(QString workString) {
    QList<QListWidgetItem *> items = this->ui->queueListWidget->findItems(workString, Qt::MatchExactly);
    foreach (QListWidgetItem * item, items) { delete item; }
    //! std::cerr << QTime::currentTime().toString().toStdString() << std::endl << std::flush;
}

void Widget::updateProgress(int progressIN)
{
    progress = progressIN;
    ui->progressBar->setValue(progress);
    //print_time_and_memory(); //! MEM USAGE
}

void Widget::updateStatus(QString statusIN)
{
    status = statusIN;
    ui->progressLabel->setText(status);
}

//!
//! originally from the merge reads >> merge tab
//!

void Widget::MERGEsearchReadsInfileR1() {
    QString fileName = QFileDialog::getOpenFileName(this,
                                                    tr("select an alignment file"),
                                                    ".",
                                                    tr("compressed SAM files (*.bam)"));
    if (!fileName.isEmpty()) { ui->MERGEr1_lineEdit->setText(fileName); }
}

void Widget::MERGEsearchReadsInfileR2() {
    QString fileName = QFileDialog::getOpenFileName(this,
                                                    tr("select an alignment file"),
                                                    ".",
                                                    tr("compressed SAM files (*.bam)"));
    if (!fileName.isEmpty()) { ui->MERGEr2_lineEdit->setText(fileName); }

}

void Widget::MERGEsearchMatrixOutfile() {
    QString fileName = QFileDialog::getSaveFileName(this,
                                                    tr("select an output file"),
                                                    ".",
                                                    tr("reads matrix (*.txt)"));
    if (!fileName.isEmpty()) { ui->MERGEout_lineEdit->setText(fileName); }

}

void Widget::MERGEsendProcess() {
    QString prefix = "RUNmergeReadPairs";
    QString maxnumReads = ui->MERGEnumReads_spinBox->text();
    QString readsInfileR1 = ui->MERGEr1_lineEdit->text();
    QString readsInfileR2 = ui->MERGEr2_lineEdit->text();
    QString matrixOutfile = ui->MERGEout_lineEdit->text();
    QStringList temp;
    temp << prefix << readsInfileR1 << readsInfileR2 << matrixOutfile << maxnumReads;
    QString messageSent = temp.join("|");
    this->ui->queueListWidget->addItem(messageSent);
    //! std::cerr << QTime::currentTime().toString().toStdString() << std::endl << std::flush;
    this->worker.addWork(messageSent);
}

//!
//! from map >> create fragments tab
//!

void Widget::CREATEFRAGMENTSsearchFastaFile() {
    QString fileName = QFileDialog::getOpenFileName(this,
                                                    tr("select a genome fasta file"),
                                                    ".",
                                                    tr("FASTA files (*.fasta)"));
    if (!fileName.isEmpty()) { ui->CREATEFRAGMENTSfasta_lineEdit->setText(fileName); }
}

void Widget::CREATEFRAGMENTSsearchFragListFile() {
    QString fileName = QFileDialog::getSaveFileName(this,
                                                    tr("select a file for the fragment annotations"),
                                                    ".",
                                                    tr("plain text (*.txt)"));
    if (!fileName.isEmpty()) { ui->CREATEFRAGMENTSfragList_lineEdit->setText(fileName); }
}

void Widget::CREATEFRAGMENTSsendProcess() {
    QString prefix = "RUNcreateFragments";
    QString restrictionPatternOrWinSize = ui->CREATEFRAGMENTSresPat_lineEdit->text();
    QString fastaFile = ui->CREATEFRAGMENTSfasta_lineEdit->text();
    QString fragListFile = ui->CREATEFRAGMENTSfragList_lineEdit->text();
    QStringList temp;
    temp << prefix << restrictionPatternOrWinSize << fastaFile << fragListFile;
    QString messageSent = temp.join("|");
    this->ui->queueListWidget->addItem(messageSent);
    //! std::cerr << QTime::currentTime().toString().toStdString() << std::endl << std::flush;
    this->worker.addWork(messageSent);
}

//!
//! from map >> map tab
//!

void Widget::MAPREADSsearchFragListFile() {
    QString fileName = QFileDialog::getOpenFileName(this,
                                                    tr("select a file for the fragment annotations"),
                                                    ".",
                                                    tr("plain text (*.txt)"));
    if (!fileName.isEmpty()) { ui->MAPREADSfragList_lineEdit->setText(fileName); }
}

void Widget::MAPREADSsearchPreMatrixFile() {
    QString fileName = QFileDialog::getOpenFileName(this,
                                                    tr("select a read matrix file"),
                                                    ".",
                                                    tr("plain text (*.txt)"));
    if (!fileName.isEmpty()) { ui->MAPREADSpreMatrix_lineEdit->setText(fileName); }
}

void Widget::MAPREADSsearchPostMatrixFile() {
    QString fileName = QFileDialog::getSaveFileName(this,
                                                    tr("select a matrix file with mappings added"),
                                                    ".",
                                                    tr("plain text (*.txt)"));
    if (!fileName.isEmpty()) { ui->MAPREADSpostMatrix_lineEdit->setText(fileName); }
}

void Widget::MAPREADSsearchReducedMatrixFile() {
    QString fileName = QFileDialog::getSaveFileName(this,
                                                    tr("select a file for the reduced matrix"),
                                                    ".",
                                                    tr("plain text (*.txt)"));
    if (!fileName.isEmpty()) { ui->MAPREADSreducedMatrix_lineEdit->setText(fileName); }
}

void Widget::MAPREADSsendProcess() {
    QString prefix = "RUNmapReadPairs";
    QString fragListFile = ui->MAPREADSfragList_lineEdit->text();
    QString preMatrixFile = ui->MAPREADSpreMatrix_lineEdit->text();
    QString postMatrixFile = ui->MAPREADSpostMatrix_lineEdit->text();
    QString reducedMatrixFile = ui->MAPREADSreducedMatrix_lineEdit->text();
    QString doReducedMatrix = "n";
    if (ui->MAPREADSreduceMatrixCheck->isChecked()) { doReducedMatrix = "y"; }
    QString useFilterStr = "n";
    QString filterSelfCircleStr = "n";
    QString inwardClose = "1000";
    QString outwardClose = "25000";
    if (ui->MAPREADSuseFilterCheck->isChecked()) {
        useFilterStr = "y";
        if (ui->MAPREADSselfCircleCheck->isChecked()) { filterSelfCircleStr = "y"; } //! NOTE THAT THIS CHECKBOX IS HIDDEN AT THE MOMENT
        inwardClose = ui->MAPREADSinwardCloseBox->text();
        outwardClose = ui->MAPREADSoutwardCloseBox->text();
    }
    QStringList temp;
    temp << prefix << fragListFile << preMatrixFile << postMatrixFile << doReducedMatrix << reducedMatrixFile << useFilterStr << filterSelfCircleStr << inwardClose << outwardClose;
    QString messageSent = temp.join("|");
    this->ui->queueListWidget->addItem(messageSent);
    //std::cerr << QTime::currentTime().toString().toStdString() << " " << messageSent.toStdString() << std::endl << std::flush;
    this->worker.addWork(messageSent);
}

//!
//! from tracks >> add tracks tab
//!

void Widget::ADDTRACKSsetFragmentView(QString header)
{
    if ((ui->ADDTRACKStrackList->count() == 0) && (header != "")) { this->ADDTRACKSheaderAtBeginning = header; }
    else { this->ADDTRACKScurrentHeader = header; }
    ui->ADDTRACKStrackList->clear();
    if (header != "") {
        QStringList headerList = header.split("|");
        this->ui->ADDTRACKStrackList->addItems(headerList);
    }
}

void Widget::ADDTRACKSremoveGenomeAnnotationFiles(QString fileName)
{
    QList<QListWidgetItem *> items = this->ui->ADDTRACKSgenomeAnnotationListWidget->findItems(fileName, Qt::MatchExactly);
    foreach (QListWidgetItem * item, items) { delete item; }
}

void Widget::ADDTRACKSremoveCountFeatureFiles(QString fileName)
{
    QList<QListWidgetItem *> items = this->ui->ADDTRACKScountFeaturesListWidget->findItems(fileName, Qt::MatchExactly);
    foreach (QListWidgetItem * item, items) { delete item; }
}

void Widget::ADDTRACKSremoveDensityFeatureFiles(QString fileName)
{
    QList<QListWidgetItem *> items = this->ui->ADDTRACKSdensityFeaturesListWidget->findItems(fileName, Qt::MatchExactly);
    foreach (QListWidgetItem * item, items) { delete item; }
}

void Widget::ADDTRACKSremoveDNAmethylationFiles(QString fileName)
{
    QList<QListWidgetItem *> items = this->ui->ADDTRACKSDNAmethylationListWidget->findItems(fileName, Qt::MatchExactly);
    foreach (QListWidgetItem * item, items) { delete item; }
}

void Widget::ADDTRACKSloadFragments() {
    QString fileName = QFileDialog::getOpenFileName(this,
                                                    tr("select a fragment list file"),
                                                    ".",
                                                    tr("from HiCdat: create fragments (*.txt)"));
    if (!fileName.isEmpty()) {
        QString command = "RUNaddTracks|LOADFRAGMENTS|" + fileName;
        this->worker.addWork(command);
    }
}

void Widget::ADDTRACKSsaveFragments() {
    QString fileName = QFileDialog::getSaveFileName(this,
                                                    tr("save the fragment list file"),
                                                    ".",
                                                    tr("plain text (*.txt)"));
    if (!fileName.isEmpty()) {
        QString command = "RUNaddTracks|SAVEFRAGMENTS|" + fileName;
        this->worker.addWork(command);
    }
}

void Widget::ADDTRACKScloseFragments() {
    bool doClose = true;
    if ((this->ADDTRACKSheaderAtBeginning != this->ADDTRACKScurrentHeader) && (this->ADDTRACKScurrentHeader != "")) {
        // ask if you want to save it - code from the help page of QMessageBox
        QMessageBox msgBox;
        msgBox.setText("The document has been modified.");
        msgBox.setInformativeText("Do you want to save your changes?");
        msgBox.setStandardButtons(QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
        msgBox.setDefaultButton(QMessageBox::Save);
        int ret = msgBox.exec();
        switch (ret) {
            case QMessageBox::Save:
                this->ADDTRACKSsaveFragments();
                break;
            case QMessageBox::Discard:
                break;
            case QMessageBox::Cancel:
                doClose = false;
                break;
            default:
                break;
        }
    }
    if (doClose) {
        QString command = "RUNaddTracks|CLOSEFRAGMENTS";
        this->worker.addWork(command);
    }
}

void Widget::ADDTRACKSaddGenomeAnnotationFiles()
{
    QStringList fileNames = QFileDialog::getOpenFileNames(this,
                                                    tr("select genome annotation files"),
                                                    ".",
                                                    tr("GFF files (*.gff);;GTF files (*.gtf)"));
    if (!fileNames.isEmpty()) { this->ui->ADDTRACKSgenomeAnnotationListWidget->addItems(fileNames); }
}

void Widget::ADDTRACKSaddCountFeatureFiles()
{
    QStringList fileNames = QFileDialog::getOpenFileNames(this,
                                                    tr("select count feature files"),
                                                    ".",
                                                    tr("BAM files (*.bam)"));//;;GTF files (*.gtf);;GFF files (*.gff)
    if (!fileNames.isEmpty()) { this->ui->ADDTRACKScountFeaturesListWidget->addItems(fileNames); }
}

void Widget::ADDTRACKSaddDensityFeatureFiles()
{
    QStringList fileNames = QFileDialog::getOpenFileNames(this,
                                                    tr("select density feature files"),
                                                    ".",
                                                    tr("BAM files (*.bam)"));//;;GTF files (*.gtf);;GFF files (*.gff)
    if (!fileNames.isEmpty()) { this->ui->ADDTRACKSdensityFeaturesListWidget->addItems(fileNames); }
}

void Widget::ADDTRACKSaddDNAmethylationFiles()
{
    QStringList fileNames = QFileDialog::getOpenFileNames(this,
                                                    tr("select DNA methylation files"),
                                                    ".",
                                                    tr("see tooltip (*.txt)"));
    if (!fileNames.isEmpty()) { this->ui->ADDTRACKSDNAmethylationListWidget->addItems(fileNames); }
}

void Widget::ADDTRACKSsendProcess() {
    QString prefix = "RUNaddTracks";
    QStringList temp;
    QStringList innerTemp;
    QListWidgetItem* item;
    temp << prefix << "PROCESSTRACKS";
    innerTemp.clear();
    for(int i = 0; i < ui->ADDTRACKSgenomeAnnotationListWidget->count(); ++i) {
        item = ui->ADDTRACKSgenomeAnnotationListWidget->item(i);
        innerTemp << item->text();
    }
    temp << innerTemp.join(";");
    innerTemp.clear();
    for(int i = 0; i < ui->ADDTRACKScountFeaturesListWidget->count(); ++i) {
        item = ui->ADDTRACKScountFeaturesListWidget->item(i);
        innerTemp << item->text();
    }
    temp << innerTemp.join(";");
    innerTemp.clear();
    for(int i = 0; i < ui->ADDTRACKSdensityFeaturesListWidget->count(); ++i) {
        item = ui->ADDTRACKSdensityFeaturesListWidget->item(i);
        innerTemp << item->text();
    }
    temp << innerTemp.join(";");
    innerTemp.clear();
    for(int i = 0; i < ui->ADDTRACKSDNAmethylationListWidget->count(); ++i) {
        item = ui->ADDTRACKSDNAmethylationListWidget->item(i);
        innerTemp << item->text();
    }
    temp << innerTemp.join(";");
    QString messageSent = temp.join("|");
    this->ui->queueListWidget->addItem(messageSent);
    //! std::cerr << QTime::currentTime().toString().toStdString() << std::endl << std::flush;
    this->worker.addWork(messageSent);
}

//!
//! from map >> create R-code tab
//!

void Widget::CREATERCODEsearchFastaFile() {
    QString fileName = QFileDialog::getOpenFileName(this,
                                                    tr("select a genome fasta file"),
                                                    ".",
                                                    tr("FASTA files (*.fasta)"));
    if (!fileName.isEmpty()) { ui->CREATERCODEfasta_lineEdit->setText(fileName); }
}

void Widget::CREATERCODEsearchRfile() {
    QString fileName = QFileDialog::getSaveFileName(this,
                                                    tr("select a file for the organism specific R-code"),
                                                    ".",
                                                    tr("R-Script (*.R)"));
    if (!fileName.isEmpty()) { ui->CREATERCODErFile_lineEdit->setText(fileName); }
}

void Widget::CREATERCODEsendProcess() {
    QString prefix = "RUNcreateRcode";
    QString restrictionPatternOrWinSize = ui->CREATERCODEresPat_lineEdit->text();
    QString fastaFile = ui->CREATERCODEfasta_lineEdit->text();
    QString rFile = ui->CREATERCODErFile_lineEdit->text();
    QStringList temp;
    temp << prefix << restrictionPatternOrWinSize << fastaFile << rFile;
    QString messageSent = temp.join("|");
    this->ui->queueListWidget->addItem(messageSent);
    //! std::cerr << QTime::currentTime().toString().toStdString() << std::endl << std::flush;
    this->worker.addWork(messageSent);
}
