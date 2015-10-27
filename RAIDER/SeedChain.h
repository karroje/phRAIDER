#ifndef SCANER_SEEDCHAIN
#define SCANER_SEEDCHAIN

#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <stdexcept>
#include <vector>
#include "Family.h"
#include <unordered_map>
#include <ctime>

using namespace std;

typedef unordered_map<size_t, LmerVector*> LmerMap;

int p = getpid();


size_t baseToInt(char b) {
	switch (toupper(b)) {
	case 'A':
		return 0;
	case 'C':
		return 1;
	case 'G':
		return 2;
	case 'T':
		return 3;
	default:
		throw std::out_of_range("Unknown base type");
		return -1;
	}
}

size_t seedToInt(const char* seed, const seqan::CharString &mask) {
	//TODO use bitset or vector<bool> to allower for potentially large bitsets
	size_t result = 0;
	uint L = seqan::length(mask);
	for (uint i = 0; i < L; i++) {
		if (mask[i] == '1') {
			result <<= 2;
			result |= baseToInt(seed[i]);
		}
	}
	return result;
}

bool getNextSeed(seqan::Dna5String &sequence, uint &index, uint seqLength, uint L, char* seed) {
	seed[L - 1] = 0;
	for (uint i = 0; i < L && index + i < seqLength; i++) {
		if (toupper((char) sequence[index + i]) == 'N') {
			index = index + i + 1;
			i = -1;
		} else {
			seed[i] = sequence[index + i];
		}
	}
	return seed[L - 1];
}

bool isNewFamily(Family* fam, LmerVector* v) {
	assert(!v->getFamily());
	if (!fam || fam->size() > 2) {
		return true;
	}
	assert(fam->size() == 2);
	LmerVector* suffix = fam->getSuffix();
	uint s1 = suffix->front();
	uint s2 = suffix->back();
	uint prev = v->front();
	uint index = v->back();
	// Ensure both Lmer pairs are the same distance apart
	if (index - prev != s2 - s1) {
		return true;
	}
	// Ensure first Lmer pair overlaps second Lmer pair
	if (index != s2 + 1) {
		return true;
	}
	return false;
}

LmerVector* getAndInsert(LmerMap &lmers, size_t seed, uint index) {
	LmerVector*& v = lmers[seed];
	if (v == 0) {
		v = new LmerVector(index);
	} else {
		v->push_back(index);
	}
	return v;
}

Family* splitRepeatsByLmer(Family* fam, LmerVector *v, bool keepV) {
	assert(fam == v->getFamily());

	vector<LmerVector*> *lmers = fam->getLmers();
	Family *newFam = new Family;

	uint oldLength = lmers->size();
	uint offset = oldLength;
	for (uint i = 0; i < oldLength; i++) {
		LmerVector *u = lmers->at(i);
		if (u == v) {
			assert(i != 0);
			offset = i;
		}
		if ((!keepV && i >= offset) || i > offset) {
			newFam->adopt(u);
		}
	}
	lmers->resize(keepV ? offset + 1 : offset);
	lmers->shrink_to_fit();

	fam->setLast(fam->getSuffix());
	fam->setExpectedEnd(fam->getSuffix()->back() + 1);

	if (!keepV) {
		newFam->setLast(v);
		newFam->setExpectedEnd(v->back() + newFam->repeatLength(1));
	}

	return newFam;
}

bool repeatExpects(Family *fam, LmerVector* v) {
	const uint index = v->back();
	const uint length = fam->repeatLength(1);
	uint dist_to_end = fam->getExpectedEnd() - index;
	if (index > fam->getExpectedEnd() || dist_to_end > length) {
		return false;
	}
	const int pos = length - dist_to_end;
	if (pos >= 0 && pos < (int) length && fam->getLmers()->at(pos) == v) {
		return true;
	}
	return false;
}

bool fragmentSplit(LmerVector* v, uint L, vector<Family*> &families) {
	Family* fam = v->getFamily();
	assert(fam);
	if (fam->lastRepeatComplete()) {
		if (v != fam->getPrefix()) {
			Family* newFam = splitRepeatsByLmer(fam, v, false);
			families.push_back(newFam);
			return true;
		}
	} else if (v->back() > fam->getLastIndex() + L) {
		Family* newFam = splitRepeatsByLmer(fam, fam->getLast(), true);
		families.push_back(newFam);
		return true;
	}
	return false;
}


uint maxLength(const vector<seqan::CharString> &masks) {
	uint max = 0;
	for (uint i = 0; i < masks.size(); i++) {
		uint len = seqan::length(masks[i]);
		if (len > max) {
			max = len;
		}
	}
	return max;
}

bool isPrefix(LmerVector *v) {
	return v->getFamily()->getPrefix() == v;
}

void tieLooseEnds(vector<Family*> &families) {
  for (auto fam : families) {
    if (!fam->lastRepeatComplete()) {
      Family* newFam = splitRepeatsByLmer(fam, fam->getLast(), true);
      families.push_back(newFam);
    }
  }
}

void getElementaryFamilies(seqan::Dna5String &sequence, vector<seqan::CharString> &masks, vector<Family*> &families) {
	const uint seqLength = seqan::length(sequence);
	const seqan::CharString mask = masks.front();
	const uint L = seqan::length(mask);
	const uint MAX_L = maxLength(masks);

	uint index = -1;
	char unmasked[MAX_L];
	LmerMap lmers;
	Family* fam = nullptr;

	while (getNextSeed(sequence, ++index, seqLength, MAX_L, unmasked)) {
	  
		size_t seed = seedToInt(unmasked, mask);
		LmerVector* v = getAndInsert(lmers, seed, index);

		if (v->size() == 2) {
			if (isNewFamily(fam, v)) {
				fam = new Family();
				families.push_back(fam);
			}
			fam->adopt(v);
		} else if (v->size() > 2 && !fragmentSplit(v, L, families)) {
			Family* fam = v->getFamily();
			if (isPrefix(v)) {
				fam->setLast(v);
				fam->setExpectedEnd(v->back() + fam->repeatLength(1));
			} else if (repeatExpects(fam, v)) {
				fam->setLast(v);
			}
		}
	}
	tieLooseEnds(families);
}

#endif //SCANER_SEEDCHAIN
