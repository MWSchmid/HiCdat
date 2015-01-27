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

#ifndef GENOMEANNOTATIONELEMENT_H
#define GENOMEANNOTATIONELEMENT_H

#include <QtCore>

class genomeAnnotationElement
{
public:
    QString _chrom;
    uint _start;
    uint _end;

    genomeAnnotationElement();
    //~genomeAnnotationElement();

    void init(QString& chrom, uint& start, uint& end);

};

#endif // GENOMEANNOTATIONELEMENT_H
