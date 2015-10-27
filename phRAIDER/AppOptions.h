// SeedChain.h is part of phRAIDER.
//
// phRAIDER is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// phRAIDER is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with phRAIDER.  If not, see <http://www.gnu.org/licenses/>.

// Created by Carly Schaeffer, Nathan Figueroa, and John Karro

// This struct stores the options from the command line.
struct AppOptions {
  // Age. Presets for familyarray, excising, overlaps, and tieup
  int age;
  // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- carly verbose
  int verbosity;
  // Minimum repeat length
  uint min;
  // Minimum number of repeats to be significant
  uint count;
  // Whether or not to use family array
  bool family_array;
  // Whether or not to allow for excising
  bool excising;
  // Whether or not to require overlaps
  bool overlaps;
  // Whether or not to allow for modified tie up
  bool tieup;
  // Do cleanup while going
  bool proactive_split;
  // Keep track of previous family
  bool prev_fam;
  // Do we use a skip-back list?
  bool sbl;
  // Do we do a pre-scan for memory reduction?
  bool prescan;
  // File name for masked file (if specified)
  seqan::CharString mask_file;
  // File name for filter_file (if specified)
  seqan::CharString filter_file;
  
  seqan::CharString seed;
  seqan::CharString sequence_file;
  seqan::CharString output_directory;
  
  AppOptions() :
  age(1), verbosity(1), family_array(true), excising(false), overlaps(false), tieup(false), proactive_split(false), prev_fam(false), sbl(false), prescan(false), mask_file(""), filter_file("") {}
};

