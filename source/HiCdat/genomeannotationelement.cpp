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

#include "genomeannotationelement.h"

genomeAnnotationElement::genomeAnnotationElement() :
    _chrom(""),
    _start(0),
    _end(0) {}

void genomeAnnotationElement::init(QString& chrom, uint& start, uint& end) {
    this->_chrom = chrom;
    this->_start = start;
    this->_end = end;
}
