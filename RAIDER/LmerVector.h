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
