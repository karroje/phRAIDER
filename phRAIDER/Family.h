// Family.h is part of phRAIDER.
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

#ifndef FAMILY_H
#define FAMILY_H
#include <vector>
#include "LmerVector.h"
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <unordered_set>

class Family {
public:
  void adopt(LmerVector *v) {
    v->setFamily(this);
    vectors.push_back(v);
    setLast(v);
    setExpectedEnd(v->back() + 1);
  }
  
  void adopt(LmerVector *v, uint L){
    v->setFamily(this);
    if (vectors.size() > 0){
      uint off = v->front() - vectors.back()->front();
      v->setOff(off - 1);
    }
    else{
      lastSkipped = nullptr;
      v->setOff(0);
    }
    vectors.push_back(v);
    setLast(v);
    setExpectedEnd(v->back() + L);
  }
  
  uint size() const {
    if (vectors.size() == 0) {
      return 0;
    }
    return vectors.front()->size();
  }
  
  uint repeatLength(uint L) const {
    return vectors.back()->front() - vectors.front()->front() + L;
  }
  
  vector<LmerVector*>* getLmers() { return &vectors; }
  
  LmerVector* at(uint index) { return vectors.at(index); }
  LmerVector* see(uint index) const { return vectors.at(index); }
  void push_back(LmerVector* v) { vectors.push_back(v); }
  
  LmerVector* findByPos(uint pos){
    uint sumOff = 0;
    if (pos == 0){
      return vectors.front();
    }
    for (LmerVector* v : vectors){
      sumOff += v->getOff();
      if (pos == sumOff){
        return v;
      }
      sumOff++;
    }
    return vectors.front();
  }
    
  LmerVector* getPrefix() const { return vectors.front(); }
  LmerVector* getSuffix() const { return vectors.back();  }
  
  LmerVector* getLast() { return last; }
  uint getLastIndex() const { return last_index; }
  bool lastRepeatComplete() const { return last == getSuffix(); }

  void setLast(LmerVector* v) {
    assert(v->getFamily() == this);
    last = v;
    last_index = v->back();
  }
  
  LmerVector* getOneAfterLast(){
    std::vector<LmerVector*>::iterator start = std::find(vectors.begin(), vectors.end(), last);
    return *(start+1);
  }

  LmerVector* getLastSkipped(){
    return lastSkipped;
  }
  
  uint getExpectedEnd() const { return expected_end; }
  void setExpectedEnd(uint expected) { expected_end = expected; }
  
  void resetOffs(){
    for (uint i = 0; i < vectors.size(); i++){
      LmerVector* v = vectors.at(i);
      if (v != vectors.front()){
        v->setOff(v->front() - vectors.at(i-1)->front());
      }
    }
  }
  
  std::vector<LmerVector*> popSkipped() {
    //std::vector<LmerVector*> ret(skipped);
    //skipped.clear();
    std::vector<LmerVector*> skipped(vectors.size());
    std::remove_copy_if (vectors.begin(),vectors.end(),skipped.begin(), [](LmerVector* v){return !v->gotSkipped();});
    vectors.shrink_to_fit();
    skipped.shrink_to_fit();  // shrink container to new size
    resetOffs();
    return skipped;
  }
  

  
  void setSkippedRange(LmerVector* lastSeen){
    std::vector<LmerVector*>::iterator start = std::find(vectors.begin(), vectors.end(), last);
    std::vector<LmerVector*>::iterator end = std::find(vectors.begin(), vectors.end(), lastSeen);
    std::for_each(start+1, end, [](LmerVector* v){ v->setSkipped(); });
    lastSkipped = *(end-1);
  }
  
  void setRemainingSkipped(){
    setSkippedRange(vectors.back());
    vectors.back()->setSkipped();
    lastSkipped = vectors.back();
  }
  //  std::vector<LmerVector*>::iterator start = std::find(vectors.begin(), vectors.end(), last);
  //  std::move(start+1, vectors.end(), std::back_inserter(skipped));
  //  vectors.erase(start+1, vectors.end()); // no longer part of backbone of this family
  //  vectors.shrink_to_fit();
  //}
  
  
  //void moveSkippedRange(LmerVector* lastSeen){
  //  std::vector<LmerVector*>::iterator start = std::find(vectors.begin(), vectors.end(), last);
  //  std::vector<LmerVector*>::iterator end = std::find(vectors.begin(), vectors.end(), lastSeen);
  //  std::move(start+1, end, std::back_inserter(skipped));
  //  vectors.erase(start+1, end); // no longer part of backbone of this family
  //  vectors.shrink_to_fit();
  //  lastSeen->setOff(lastSeen->front() - last->front());  // adjust offset
  //}
  
  //void removeOne(LmerVector* v){
  //  std::vector<LmerVector*>::iterator start = std::find(vectors.begin(), vectors.end(), v);
  //  (*(start+1))->setOff((*(start+1))->front() - (*(start-1))->front());
  //  vectors.erase(start, start+1); // no longer part of backbone of this family
  //  vectors.shrink_to_fit();
  //}

  
  friend std::ostream &operator<< (std::ostream &output, const Family &F){
    output  << "Size: " << F.vectors.size()
    << " Expected End: " << F.expected_end
    << "\t\tLast:\t" << *(F.last) << std::endl
    << "\t\tPrefix:\t" << *(F.getPrefix()) << std::endl
    << "\t\tSuffix:\t" << *(F.getSuffix()) << std::endl;
    for(uint i = 0; i < F.vectors.size(); i++){
      output << "\t\t\t"<< i << ":\t" << *(F.see(i)) << std::endl;
    }
    return output;
  }
  
  LmerVector* last;
  LmerVector* lastSkipped;
  uint last_index;
  uint expected_end;
  std::vector<LmerVector*> vectors;
  // std::vector<LmerVector*> skipped;

  // Instances which are "removed" in post-processing, by LmetVector index.
  // Removal only occurs in post-processing -- not relevant to RAIDER algorithm.
  unordered_set<uint> excluded;  
  
};

#endif //FAMILY_H
