// Family.h is part of phRAIDER.
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
#ifndef SCANER_FAMILY
#define SCANER_FAMILY
#include <vector>
#include "LmerVector.h"
#include <iostream>
#include <assert.h>

using namespace std;

class Family {
public:

	void adopt(LmerVector *v) {
		v->setFamily(this);
		vectors.push_back(v);
		setLast(v);
		setExpectedEnd(v->back() + 1);
	}

	uint size() const {
		if (vectors.size() == 0) {
			return 0;
		}
		return vectors.front()->size();
	}

	uint repeatLength(uint L) const {
		return vectors.size() + L - 1;
	}

	vector<LmerVector*>* getLmers() {
		return &vectors;
	}

	LmerVector* at(uint index) {
		return vectors.at(index);
	}

	void push_back(LmerVector* v) {
		vectors.push_back(v);
	}

	LmerVector* getPrefix() const {
		return vectors.front();
	}

	LmerVector* getSuffix() const {
		return vectors.back();
	}

	LmerVector* getLast() {
		return last;
	}

	uint getLastIndex() {
		return last_index;
	}

	void setLast(LmerVector* v) {
		assert(v->getFamily() == this);
		last = v;
		last_index = v->back();
	}

	bool lastRepeatComplete() const {
		return last == getSuffix();
	}

	void setExpectedEnd(uint expected) {
		expected_end = expected;
	}

	uint getExpectedEnd() {
		return expected_end;
	}


	LmerVector* last;
	uint last_index;
	uint expected_end;
	vector<LmerVector*> vectors;
};

#endif //SCANER_FAMILY
