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
#ifdef Q_WS_MAC
#include <QtWidgets/QApplication> // Qt5 instead of 4
#endif
#ifndef Q_WS_MAC
#include <QtGui/QApplication> // Qt 4
#endif

int main(int argc, char *argv[])
{
    #ifdef Q_WS_MAC
    QString path = QString(argv[0]);
    QString end = path.split("/").last();
    path.remove(path.length() - end.length(), end.length());
    QDir dir(path);
    dir.cdUp();
    dir.cd("PlugIns");
    QApplication::setLibraryPaths(QStringList(dir.absolutePath()));
    #endif
    QApplication a(argc, argv);
    Widget w;
    w.show();
    
    return a.exec();
}
