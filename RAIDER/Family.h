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
