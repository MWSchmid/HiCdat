#-------------------------------------------------
#
# Project created by QtCreator 2014-09-24T09:42:08
#
#-------------------------------------------------

QT       += core gui

TARGET = HiCdatPre
TEMPLATE = app


SOURCES += main.cpp \
    widget.cpp \
    workprocessor.cpp \
    annotation_reader.cpp \
    genomefragment.cpp \
    genomefragmentcollection.cpp \
    genomeannotationelement.cpp

HEADERS  += widget.h\
            helperClasses.h \
    workprocessor.h \
    annotation_reader.h \
    genomefragment.h \
    genomefragmentcollection.h \
    genomeannotationelement.h \
    printTimeAndMem.h

FORMS    += widget.ui

# note that the mac specificity is not only due to mac but also due to Qt version 5... instead of 4...
macx {
    QT += widgets
    LIBS += -lZ
    INCLUDEPATH +=  ../zlib-1.2.8 \
}

