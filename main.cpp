// main.cpp is part of phRAIDER.
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

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/gff_io.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <numeric>

#include "SeedChain.h"
#include "Family.h"
#include "LmerVector.h"

using namespace std;

typedef pair<uint, string> idThreshold;

const string OUTPUT_ELEMENTARY_REPEATS_FILENAME = "elements";
const string OUTPUT_REPEAT_FAMILIES_FILENAME = "families";
const string OUTPUT_COMPOSITES_FILENAME = "composites";
const string OUTPUT_SUMMARY_FILENAME = "summary_info";


typedef pair<int,int> Interval;

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------
seqan::ArgumentParser::ParseResult parseCommandLine(AppOptions & options, int argc, char const ** argv) {
  // Setup ArgumentParser.
  seqan::ArgumentParser parser("RAIDER2");
  // Set short description, version, and date.
  setShortDescription(parser, "RAIDER - Rapid Ab Initio Detection of Elementary Repeats");
  setVersion(parser, "2.0");
  setDate(parser, "June 2015");
  
  // Define usage line and long description.
  addUsageLine(parser, "[\\fIOPTIONS\\fP]  \"\\fISEQUENCE_FILE\\fP\"  \"\\fIOUTPUT_DIRECTORY\\fP\"");
  addDescription(
                 parser,
                 "RAIDER2 parses the given sequence file using the supplied mask (spaced seed) to identify de novo repeats. Minimum repeat size and other options can be configured as described below.");
  
  // We require two arguments.
  //addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "MASK"));
  addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "SEQUENCE_FILE"));
  addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "OUTPUT_DIRECTORY"));
  
  addOption(
            parser,
            seqan::ArgParseOption("s", "seed", "Spaced seed/mask to use. Defaults to 111110011111110001111111000000000000011111.",
                                  seqan::ArgParseOption::STRING));

  addOption(
	    parser,
	    seqan::ArgParseOption("mf", "mask_file", "Masked sequence file. Default: none",
				  seqan::ArgParseOption::STRING));

  addOption(
	     parser,
	     seqan::ArgParseOption("ff", "filter_file", "Filter file. Default: none",
				   seqan::ArgParseOption::STRING));
  
  addOption(
            parser,
            seqan::ArgParseOption("m", "min", "Minimum repeat length. Defaults to pattern length.",
                                  seqan::ArgParseOption::INTEGER));
  addOption(
            parser,
            seqan::ArgParseOption("c", "count", "Minimum number of repeats in a family. Defaults to 5.",
                                  seqan::ArgParseOption::INTEGER));
  addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
  addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
  addOption(parser, seqan::ArgParseOption("vv", "verbose+", "Enable extremely verbose output."));
  addOption(
            parser,
            seqan::ArgParseOption("a", "age", "Age of raiderv2. Defaults to 1.",
                                  seqan::ArgParseOption::INTEGER));
  addOption(parser, seqan::ArgParseOption("na", "noarray", "Disable family array (enabled by default)."));
  addOption(parser, seqan::ArgParseOption("e", "excise", "Enable excising (disabled by default)."));
  addOption(parser, seqan::ArgParseOption("no", "overlaps", "Require overlaps (not required by default)."));
  addOption(parser, seqan::ArgParseOption("t", "tieup", "Enable alternate tie up (disabled by default)."));
  addOption(parser, seqan::ArgParseOption("ps", "prosplit", "Enable proactive splitting(disabled by default)."));
  addOption(parser, seqan::ArgParseOption("pf", "prevfam", "Enable pointers to prev family (disabled by default)."));
  addOption(parser, seqan::ArgParseOption("sbl", "skipbacklist", "Enable skip back list (disabled by default)."));
  addOption(parser, seqan::ArgParseOption("p", "prescan", "Enable prescan."));
  
  // Add Examples Section.
  addTextSection(parser, "Examples");
  addListItem(parser, "\\fBraider\\fP \\fB-v\\fB -s \\fI1110110111\\fP \\fIchr23.fasta\\fP \"\\fIchr23_out\\fP\"",
              "Call with mask \"1110110111\" and verbose output.");
  
  // Parse command line.
  seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
  
  // Only extract  options if the program will continue after parseCommandLine()
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res;
  
  // Extract option values.
  if (isSet(parser, "quiet"))
    options.verbosity = 0;
  if (isSet(parser, "verbose"))
    options.verbosity = 2;
  if (isSet(parser, "verbose+"))
    options.verbosity = 3;
  
  //seqan::getArgumentValue(options.seed, parser, 0);
  seqan::getArgumentValue(options.sequence_file, parser, 0);
  seqan::getArgumentValue(options.output_directory, parser, 1);
  
  
  if (isSet(parser, "count"))
    seqan::getOptionValue(options.count, parser, "count");
  else
    options.count = 5;
  
  if (isSet(parser, "seed"))
    seqan::getOptionValue(options.seed, parser, "seed");
  else
    options.seed = "111110011111110001111111000000000000011111";

  if (isSet(parser, "mask_file"))
    seqan::getOptionValue(options.mask_file, parser, "mask_file");
  else
    options.mask_file = "";

  if (isSet(parser, "filter_file"))
    seqan::getOptionValue(options.filter_file, parser, "filter_file");
  else
    options.filter_file = "";
  
  if (isSet(parser, "min"))
    seqan::getOptionValue(options.min, parser, "min");
  else
    options.min = seqan::length(options.seed);
  
  if (isSet(parser, "age")){
    seqan::getOptionValue(options.age, parser, "age");
    options.family_array = options.age > 0 ? true : false;
    options.excising = options.age > 1 ? true : false;
  }
  else
    options.age = 1;
  
  // If any switch is explicitly set, this overrides age
  if (isSet(parser, "noarray")){
    options.family_array = false;
  }
  
  if (isSet(parser, "excise")){
    options.excising = true;
  }
  
  if (isSet(parser, "overlaps")){
    options.overlaps = true;
  }
  
  if (isSet(parser, "tieup")){
    options.tieup = true;
  }
  
  if (isSet(parser, "prosplit")){
    options.proactive_split = true;
  }
  
  if (isSet(parser, "prevfam")){
    options.prev_fam = true;
  }

  if (isSet(parser, "sbl")){
    options.sbl = true;
  }

  if (isSet(parser, "prescan")){
    options.prescan = true;
  }
#ifdef PRE    // Used for the pre-phRAIDER version.
  options.prescan = true;
#endif
  
  // Ensure a trailing /
  if (options.output_directory[seqan::length(options.output_directory) - 1] != '/') {
    seqan::append(options.output_directory, "/");
  }

  return seqan::ArgumentParser::PARSE_OK;
}

/**
 * Given a sequence stream, read in all sequences and concatenate into a master sequence
 */
bool concatenateSequences(seqan::SequenceStream &seqStream, vector<idThreshold> &thresholds,
                          vector<pair<uint,uint>>& seq_coords, unordered_map<string,uint>& chr2start,
			  seqan::Dna5String &outSequence, int verbosity, int L) {
  thresholds.clear();
  seq_coords.clear();
  chr2start.clear();
  seqan::Dna5String buffer;
  for (int i=0; i < L; i++) // Better way to do this?
    seqan::append(buffer, "N");
  
  seqan::CharString id;
  if (seqan::readRecord(id, outSequence, seqStream) != 0) {
    return false;
  }
  if (verbosity > 0) {
    cout << "Preparing " << id << endl;
  }
  uint currentThreshold = seqan::length(outSequence);
  thresholds.push_back(make_pair(currentThreshold + L, toCString(id)));
  seqan::Dna5String other;
  seq_coords.push_back(make_pair<uint,uint>(0, seqan::length(outSequence)));
  chr2start[toCString(id)] = 0;



  while (seqan::readRecord(id, other, seqStream) == 0) {
    if (verbosity > 0) {
      cout << "Preparing " << id << endl;
    }
    seqan::append(outSequence, buffer);
    currentThreshold += seqan::length(other);

    uint start = seqan::length(outSequence);
    uint finish = start + length(other);
    seq_coords.push_back(make_pair(start,finish));
    thresholds.push_back(make_pair(currentThreshold + L, toCString(id)));
    chr2start[toCString(id)] = start;
    seqan::append(outSequence, other);
  }

  return true;
}

/**
 * Given a file name, read all sequences from file and compile into master sequence. Also return
 * list of thresholds: indices into the master sequence paired with sequence IDs to mark where each
 * individual sequence begins.
 */
bool getSequence(seqan::CharString file, seqan::Dna5String &sequence, vector<idThreshold> &thresholds, vector<pair<uint,uint> >& seq_coords, unordered_map<string,uint>& chr2start, int verbosity, int L) {
  if (verbosity > 0) {
    cout << "Loading sequence..." << endl;
  }
  
  seqan::SequenceStream seqStream(seqan::toCString(file));
  if (!seqan::isGood(seqStream) || !concatenateSequences(seqStream, thresholds, seq_coords, chr2start, sequence, verbosity, L)) {
    cout << "Error: unable to open sequence." << endl;
    return false;
  }

  return true;
}

/**
 * Given an index into the master sequence, return the associated sequence ID
 */
seqan::CharString getSeqId(uint index, vector<idThreshold> &thresholds) {
  for (uint i = 0; i < thresholds.size(); i++) {
    if (index < thresholds[i].first) {
      return thresholds[i].second;
    }
  }
  cout << "ERROR: unable to find associated sequence ID" << endl;
  return "Unknown";
}

/**
 * Given an index into the master sequence, find the corresponding original sequence and
 * return the corresponding index into that sequence
 */
uint getSeqIndex(uint index, vector<idThreshold> &thresholds) {
  uint prevThreshold = 0;
  for (uint i = 0; i < thresholds.size(); i++) {
    if (index < thresholds[i].first) {
      return index - prevThreshold;
    }
    prevThreshold = thresholds[i].first;
  }
  cout << "ERROR: unable to find correct threshold" << endl;
  // TODO throw exceptions
  return index;
}

/**
 * Print the command line arguments back to the user.
 */
void printArgs(AppOptions &options) {
  if (options.verbosity > 0) {
    cout << "__ARGUMENTS____________________________________________________________________" << endl
				<< "VERBOSITY\t" << options.verbosity << endl
				<< "MIN_LENGTH\t" << options.min << endl
				<< "MIN_COUNT\t" << options.count << endl
				<< "SPACED_SEED     \t" << options.seed <<endl
				<< "SEQUENCE_FILE\t" << options.sequence_file << endl
				<< "OUTPUT_DIRECTORY\t" << options.output_directory << endl;
  }
}

void writeFamilies(vector<Family*> &families, AppOptions &options) {
  ofstream ofile;
  string familiesFilePath(seqan::toCString(options.output_directory));
  familiesFilePath.append(OUTPUT_REPEAT_FAMILIES_FILENAME);
  
  if (options.verbosity > 0) {
    cout << "Writing families to " << familiesFilePath << endl;
  }
  
  ofile.open(familiesFilePath.c_str());
  
  ofile << "#  fam\tnum_copies\tcopy_length" << endl;
  for (uint i = 0; i < families.size(); i++) {
    Family* fam = families[i];
    if (fam->repeatLength(options.min) >= options.min && fam->size() - fam->excluded.size() >= options.count) {
      ofile << i << "\t" << fam->size() - fam->excluded.size() << "\t" << fam->repeatLength(options.min) << endl;
    }
  }
}

void writeRepeats(vector<Family*> &families, vector<idThreshold> &thresholds, AppOptions &options) {
  uint repCount = 0;
  
  ofstream ofile;
  string repeatsFilePath(seqan::toCString(options.output_directory));
  repeatsFilePath.append(OUTPUT_ELEMENTARY_REPEATS_FILENAME);
  
  if (options.verbosity > 0) {
    cout << "Writing elementary repeats to " << repeatsFilePath << endl;
  }
  
  ofile.open(repeatsFilePath.c_str());
  
  ofile << "#  fam\tele\tdir\tsequence\tstart\tend" << endl;

  for (uint i = 0; i < families.size(); i++) {
    Family* fam = families[i];
 
    if (fam->repeatLength(options.min) >= options.min && fam->size() - fam->excluded.size() >= options.count) {
      int famId = i;
      int repId = repCount;
      repCount += fam->size();
      
      LmerVector* prefix = fam->getPrefix();
      for (uint i = 0; i < prefix->size(); i++) {
	if (fam->excluded.find(i) != fam->excluded.end())
	  continue;
        repId++;
        uint index = (*prefix)[i];
        uint trueIndex = getSeqIndex(index, thresholds);
        uint length = fam->repeatLength(options.min);
        ofile << famId << "\t" << repId << "\t" << "1" << "\t" << getSeqId(index, thresholds) << "\t"
        << trueIndex << "\t" << trueIndex + length << endl;
      }
    }
    
  }
}


void writeSummary(vector<Family*> &families, AppOptions &options) {
  uint maxLen = 0;
  uint maxSize = 0;
  uint maxSecondSize = 0;
  uint repCount = 0;
  uint famCount = 0;
  
  ofstream ofile;
  string summaryFilePath(seqan::toCString(options.output_directory));
  summaryFilePath.append(OUTPUT_SUMMARY_FILENAME);
  
  if (options.verbosity > 0) {
    cout << "Writing summary to " << summaryFilePath << endl;
  }
  
  ofile.open(summaryFilePath.c_str());
  
  for (uint i = 0; i < families.size(); i++) {
    Family* fam = families[i];
    uint fam_size = fam->size() - fam->excluded.size();
    if (fam->repeatLength(options.min) >= options.min && fam_size >= options.count) {
      famCount++;
      repCount += fam_size;
      
      if (fam->size() > maxSize) {
        maxSecondSize = maxSize;
        maxSize = fam_size;
      }
      
      if (fam->repeatLength(options.min) > maxLen) {
        maxLen = fam->repeatLength(options.min);
      }
    }
  }
  
  ofile << "#FAMILIES\t" << famCount << endl;
  ofile << "#REPEATS\t" << repCount << endl;
  ofile << "LONGEST_REPEAT\t" << maxLen << endl;
  ofile << "LONGEST_FAMILY\t" << maxSize << endl;
  ofile << "SECOND_LONGEST\t" << maxSecondSize << endl;
}

void writeResults(vector<Family*> &families, vector<idThreshold> &thresholds, AppOptions &options) {
  writeFamilies(families, options);
  writeRepeats(families, thresholds, options);
  writeSummary(families, options);
}
 

// Create the skip-back list for the seed.  
vector<int> createSBL(vector<seqan::CharString> seeds, bool sbl) {
  int l = accumulate(seeds.cbegin(), seeds.cend(), 0, [](int x, seqan::CharString s) {return max(x, (int)length(s));});
  set<int> S({1, l+1});

  if (not sbl) {
    vector<int> SBL(l+1);
    iota(SBL.begin(), SBL.end(), 1);
    return SBL;
  }

  set<int> V;



  for (auto i=seeds.begin(); i!=seeds.end(); i++)
    for (int j=0; j < (int)length(*i); j++)
      if ((*i)[j]=='0')
	V.insert(j);


  set<int> V2(V);
  V2.insert(-1);

  for (int i=2; i <= l; i++)
    for (set<int>::iterator j=V.begin(); j!=V.end(); j++) {
      if (V2.find(*j - i) != V2.end()) {
	S.insert(i);
	break;
      }
    }

  S.insert(l+1);
  
  vector<int> SBL(S.size());
  copy(S.begin(), S.end(), SBL.begin());
  return SBL;
}

// Remove lmers that are overlapping exons.  Remove families which are not left with a sufficient number of lmers.
// IMPORTANT: New copies of families are not fully definied.  The vectors property is properly coppied.  Nothing else is.
void filter_families(const vector<Family*>& families, unordered_map<string,uint>& chr2start, AppOptions& options, uint L) {
  assert(options.filter_file != "");
  
  vector<Interval> exon_list;
  seqan::GffStream gffIn(toCString(options.filter_file));
  seqan::GffRecord R;
  while (!atEnd(gffIn)) {
    seqan::readRecord(R, gffIn);
    uint offset = chr2start[toCString(R.ref)];
    exon_list.push_back(make_pair<uint,uint>(R.beginPos + offset, R.endPos + offset));
  }
  sort(exon_list.begin(), exon_list.end());

  vector<Family*> new_families;
  for (Family* f : families) {
    LmerVector* l = f->vectors.front();
    uint len = f->repeatLength(L);
    unordered_set<uint> filtered;
    for (uint i=0; i < l->lmers.size(); i++) {
      uint c = l->lmers[i];
      Interval I = make_pair(c, c + len);
      vector<Interval>::iterator p = lower_bound(exon_list.begin(), exon_list.end(), I);
      if ((p != exon_list.end() && p->first < I.second) || (p != exon_list.begin() && (p-1)->second >= I.first))
	f->excluded.insert(i);
    }
  }
}
      
      

  


void mask_sequence(const vector<Family*>& families, seqan::Dna5String& sequence,
		   vector<idThreshold> thresholds, vector<pair<uint,uint>>& seq_coords,
		   seqan::CharString output, int L, uint count) {
  for (Family* f : families) {
    if (f->size() - f->excluded.size() < count)
      continue;
    for (uint i=0; i < f->getPrefix()->lmers.size(); i++) {
      if (f->excluded.find(i) == f->excluded.end()) {
	uint start = f->getPrefix()->lmers[i];
	for (uint i = 0; i < f->repeatLength(L); i++)
	  sequence[start + i] = 'N';
      }
    }
  }

  ofstream out(seqan::toCString(output));

  for (unsigned i=0; i < thresholds.size(); i++) {
    uint start = seq_coords[i].first;
    uint finish = seq_coords[i].second;
    auto chromosome = infix(sequence, start, finish);
    seqan::write(out, chromosome, thresholds[i].second, seqan::Fasta());
  }
}
  

 

int main(int argc, char const ** argv) {
  AppOptions options;
  parseCommandLine(options, argc, argv);
  printArgs(options);

  vector<Family*> families;
  vector<seqan::CharString> seeds;

  seeds.push_back(options.seed);

  
  // Load sequences into master sequence, mark where each sequence ID begins within the master
  vector<idThreshold> thresholds;
  vector<pair<uint,uint>> seq_coords;
  unordered_map<string, uint> chr2start;
  int L = seqan::length(seeds[0]);   // TO DO: What is there are multiple seeds?
  seqan::Dna5String sequence;
  if (getSequence(options.sequence_file, sequence, thresholds, seq_coords, chr2start, options.verbosity, L) == false) {
    // TODO throw exceptions
    return 1;
  }

  if (options.verbosity > 0) {
    cout << "BASE PAIRS\t" << seqan::length(sequence) << endl;
  }
  

  //vector<int> SBL = createSBL(seeds, options.sbl);

  getElementaryFamilies(sequence, seeds, families, options);

  if (options.verbosity > 0) {
    cout << "Writing results elements..." << endl;
  }

  if (options.filter_file != "") 
    // IMPORTANT: This "ruins" the families -- only the vector property is left intact. (Shouldn't 
    // matter, since the the algorithm is done.)  Also a potential source of memory leaks.  (Also
    // shouldn't matter.)
    filter_families(families, chr2start, options, L);

  if (options.mask_file != "") {
    mask_sequence(families, sequence, thresholds, seq_coords, options.mask_file, L, options.count);
  }

  writeResults(families, thresholds, options);
  
  return 0;
}
