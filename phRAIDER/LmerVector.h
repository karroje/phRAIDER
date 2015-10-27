#ifndef LMER_VECTOR_H
#define LMER_VECTOR_H
#include <vector>
#include <stdexcept>
#include <iostream>

// LmerVector.h is part of phRAIDER.
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

using namespace std;

#ifndef uint
typedef unsigned int uint;
#endif

class Family;

class LmerVector {
public:
  
  LmerVector(uint position) {
    family = nullptr;
    push_back(position);
  }
  
  ~LmerVector() {
  }
  
  void push_back(uint val) {
    if(lmers.size() > 0){
      previous = back();
    }
    lmers.push_back(val);
  }
  
  uint operator [](uint i) const { return lmers[i]; }
  uint & operator [](uint i)     { return lmers[i]; }
  
  uint prev()  const { return previous; }
  uint back()  const { return lmers.back(); }
  uint front() const { return lmers.front(); }
  
  uint getOff() const { return off; }
  void setOff(uint i) { off = i;    }
  
  uint size() const { return lmers.size(); }
  
  Family* getFamily() const   { return family; }
  void setFamily(Family *fam) {
    //if (family) {
    //  prevFamily = family;
    //}
    family = fam;
    skipped = false;
  }
  
  bool gotSkipped() const { return skipped; }
  void setSkipped() { skipped = true; }
  //Family* getPrevFamily() const   { return prevFamily; }
  //void setPrevFamily(Family *fam) { prevFamily = fam;  }
  
  friend ostream &operator<<( ostream &output, const LmerVector &L) {
    output  << "(\tFront:\t" << L.front()
    << "\tBack:\t" << L.back() << "\tLocations:\t";
    for (uint loc : L.lmers){
      output << loc << " ";
    }
    output << "\tSize:\t" << L.size()
    << "\tSkipped?\t" << L.gotSkipped()
    << "\tOffset In Family\t" << L.getOff()
    << "\tFamily:\t" << L.getFamily() << "\t)";
    return output;
  }
  
  vector<uint> lmers;
  uint previous;
  uint off;
  Family* family;
  bool skipped;
  //Family* prevFamily;
};

#endif //LMER_VECTOR_H
