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

#ifndef PRINTTIMEANDMEM_H
#define PRINTTIMEANDMEM_H

//! MEM USAGE
#include <cstdlib>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unistd.h>
#include <stdio.h>
#include <map>
#include <ctime>

#include <ios> // for mem usage
#include <string> // for mem usage

void process_mem_usage(double& vm_usage, double& resident_set)
{

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   std::ifstream stat_stream("/proc/self/stat", std::ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   std::string pid, comm, state, ppid, pgrp, session, tty_nr;
   std::string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   std::string utime, stime, cutime, cstime, priority, nice;
   std::string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}

void print_time_and_memory() {
    std::time_t now = std::time(0);
    std::tm* local_time = std::localtime(&now);
    double vm, rss; //! MEM USAGE
    process_mem_usage(vm, rss);
    std::cerr << "#TIME INFO#" << '\t' << std::asctime(local_time) << std::flush;
    std::cerr << "#PROC INFO#" << '\t' << "VM:  " << vm << std::endl << "#PROC INFO#" << '\t' << "RSS: " << rss << std::endl << std::flush;
} // taken from C++ Cookbook and a website - see process_mem_usage


#endif // PRINTTIMEANDMEM_H
