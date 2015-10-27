// LmerVector.h is part of phRAIDER.
//
// RAIDER is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RAIDER is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with phRAIDER.  If not, see <http://www.gnu.org/licenses/>.

// Created by Nathan Figueroa, Carly Schaeffer, and John Karro
#ifndef LMER_VECTOR
#define LMER_VECTOR
#include <vector>
#include <stdexcept>

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

	void push_back(uint val) {
		lmers.push_back(val);
	}

	uint operator [](uint i) const {
		return lmers[i];
	}

	uint & operator [](uint i) {
		return lmers[i];
	}

	uint back() {
		return lmers.back();
	}

	uint front() {
		return lmers.front();
	}

	uint size() {
		return lmers.size();
	}

	Family* getFamily() const {
		return family;
	}

	void setFamily(Family *fam) {
		family = fam;
	}

	vector<uint> lmers;
	Family* family;
};

#endif //LMER_VECTOR
