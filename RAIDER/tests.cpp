#include "gtest/gtest.h"
#include "SeedChain.h"
#include "Family.h"
#include <vector>
#include <stdexcept>


using namespace std;


TEST(SeedChain, maxLength) {
	vector<seqan::CharString> masks;
	masks.push_back("0101");
	masks.push_back("010");
	EXPECT_EQ(maxLength(masks), 4);
}

TEST(LmerVector, setFamily_getFamily) {
	LmerVector r(0);
	EXPECT_EQ(r.getFamily(), nullptr);
	Family fam;
	fam.adopt(&r);
	r.setFamily(&fam);
	EXPECT_EQ(r.getFamily(), &fam);
}

TEST(LmerVector, access_assign_push_back) {
	LmerVector r(5);
	r.push_back(7);
	EXPECT_EQ(r[0], 5);
	EXPECT_EQ(r[1], 7);
	r[1] = 9;
	EXPECT_EQ(r[1], 9);
}

TEST(Family, size) {
	LmerVector v(0);
	v.push_back(7);

	Family fam;
	fam.adopt(&v);
	EXPECT_EQ(fam.size(), 2);
}


TEST(SeedChain, baseToInt) {
	EXPECT_NO_THROW({
		size_t result = baseToInt('c');
		EXPECT_EQ(result, 1);
		result = baseToInt('G');
		EXPECT_EQ(result, 2);
	});

	ASSERT_THROW(baseToInt('z'), std::out_of_range);
}


TEST(SeedChain, seedToInt) {
	const char* seed = "ACGTAC";
	seqan::CharString mask = "110101";
	size_t iSeed = seedToInt(seed, mask);

	EXPECT_EQ(iSeed, 29);
}

TEST(SeedChain, getNextSeed) {
	seqan::Dna5String sequence = "NNNNNNACGTNCGGTANNNNN";

	const uint seqLength = seqan::length(sequence);
	const uint mLength = 4;
	uint index = 0;

	char seed[5] = {0};
	EXPECT_TRUE(getNextSeed(sequence, index, seqLength, mLength, seed));
	EXPECT_TRUE((seed[0] == 'A' && seed[1] == 'C' && seed[2] == 'G' && seed[3] == 'T'));

	EXPECT_TRUE(getNextSeed(sequence, ++index, seqLength, mLength, seed));
	EXPECT_TRUE((seed[0] == 'C' && seed[1] == 'G' && seed[2] == 'G' && seed[3] == 'T'));

	EXPECT_TRUE(getNextSeed(sequence, ++index, seqLength, mLength, seed));
	EXPECT_TRUE((seed[0] == 'G' && seed[1] == 'G' && seed[2] == 'T' && seed[3] == 'A'));

	EXPECT_FALSE(getNextSeed(sequence, ++index, seqLength, mLength, seed));
}


TEST(SeedChain, getAndInsert) {
	LmerMap lmers;
	LmerVector v(0);
	lmers[7] = &v;

	EXPECT_EQ(getAndInsert(lmers, 7, 10), &v);
	EXPECT_EQ(v.back(), 10);

	LmerVector *u = getAndInsert(lmers, 9, 100);
	EXPECT_EQ(u->front(), 100);
	delete u;
}

TEST(SeedChain, isNewFamily) {
	LmerVector v(10);
	v.push_back(100);

	Family* fam = 0;

	EXPECT_TRUE(isNewFamily(fam, &v));

	fam = new Family;
	fam->adopt(&v);

	LmerVector u(11);
	u.push_back(101);

	EXPECT_FALSE(isNewFamily(fam, &u));
	delete fam;

	Family newFam;
	newFam.adopt(&v);
	v.push_back(101);

	EXPECT_TRUE(isNewFamily(&newFam, &u));
}

TEST(Family, getPrefix) {
	LmerVector v(0);
	v.push_back(7);

	Family fam;
	fam.adopt(&v);
	EXPECT_EQ(fam.getPrefix(), &v);
}

TEST(Family, adopt) {
	LmerVector v(0);
	v.push_back(7);

	Family fam;
	fam.adopt(&v);

	LmerVector u(3);
	u.push_back(9);

	fam.adopt(&u);

	EXPECT_EQ(u.getFamily(), v.getFamily());
	EXPECT_EQ(fam.getLmers()->back(), &u);
	EXPECT_EQ(fam.getExpectedEnd(), 10);
	EXPECT_EQ(fam.getLast(), &u);
}



TEST(SeedChain, splitRepeatsByLmer_keepV) {
	LmerVector v(0);
	v.push_back(10);

	Family fam;
	fam.adopt(&v);

	LmerVector u(1);
	u.push_back(11);
	fam.adopt(&u);

	LmerVector w(2);
	w.push_back(12);
	fam.adopt(&w);

	LmerVector x(3);
	x.push_back(23);
	fam.adopt(&x);

	Family* newFam = splitRepeatsByLmer(&fam, &w, true);
	EXPECT_EQ(v.getFamily(), &fam);
	EXPECT_EQ(u.getFamily(), &fam);
	EXPECT_EQ(w.getFamily(), &fam);
	EXPECT_EQ(x.getFamily(), newFam);

	EXPECT_EQ(fam.getSuffix(), &w);
	EXPECT_EQ(fam.getExpectedEnd(), 13);
	EXPECT_EQ(fam.getLast(), &w);
	delete newFam;
}



TEST(SeedChain, splitRepeatsByLmer_no_keepV) {
	LmerVector v(0);
	v.push_back(10);

	Family fam;
	fam.adopt(&v);

	LmerVector u(1);
	u.push_back(11);
	fam.adopt(&u);

	LmerVector w(2);
	w.push_back(22);
	fam.adopt(&w);

	LmerVector x(3);
	x.push_back(23);
	fam.adopt(&x);

	Family* newFam = splitRepeatsByLmer(&fam, &w, false);
	EXPECT_EQ(v.getFamily(), &fam);
	EXPECT_EQ(u.getFamily(), &fam);
	EXPECT_EQ(w.getFamily(), newFam);
	EXPECT_EQ(x.getFamily(), newFam);

	EXPECT_EQ(fam.getSuffix(), &u);
	EXPECT_EQ(fam.getExpectedEnd(), 12);
	EXPECT_EQ(fam.getLast(), &u);
	delete newFam;
}


TEST(SeedChain, fragmentSplit_missingPrefix) {
	LmerVector v(0);
	LmerVector u(1);

	Family fam;
	fam.adopt(&v);
	fam.adopt(&u);

	vector<Family*> families;
	families.push_back(&fam);
	uint L = 2;

	u.push_back(9);

	EXPECT_EQ(fam.repeatLength(1), 2);
	EXPECT_EQ(families.size(), 1);
	EXPECT_TRUE(fragmentSplit(&u, L, families));
	EXPECT_EQ(fam.repeatLength(1), 1);
	EXPECT_EQ(families.size(), 2);

}


TEST(SeedChain, fragmentSplit_missingSuffix) {
	Family fam;
	vector<Family*> families;
	families.push_back(&fam);

	LmerVector v(0);
	LmerVector u(1);
	LmerVector w(2);

	fam.adopt(&v);
	fam.adopt(&u);
	fam.adopt(&w);

	v.push_back(10);
	u.push_back(11);
	w.push_back(12);

	uint L = 1;
	fam.setLast(&u);

	EXPECT_FALSE(fragmentSplit(&w, L, families));

	EXPECT_EQ(fam.repeatLength(1), 3);
	EXPECT_EQ(families.size(), 1);

	v.push_back(20);
	fam.setLast(&v);
	fam.setExpectedEnd(23);

	w.push_back(22);
	EXPECT_TRUE(fragmentSplit(&w, L, families));

	EXPECT_EQ(fam.repeatLength(1), 1);
	EXPECT_EQ(families.size(), 2);
}


TEST(Family, repeatLength) {
	Family fam;

	LmerVector v(0);
	v.push_back(10);
	fam.adopt(&v);

	LmerVector u(1);
	u.push_back(11);
	fam.adopt(&u);

	EXPECT_EQ(fam.repeatLength(1), 2);
	EXPECT_EQ(fam.repeatLength(5), 6);
}


TEST(Family, repeatExpects) {
	Family fam;

	LmerVector v(0);
	v.push_back(10);
	fam.adopt(&v);

	LmerVector u(1);
	u.push_back(11);
	fam.adopt(&u);

	v.push_back(20);

	EXPECT_FALSE(repeatExpects(&fam, &v));

	fam.setLast(&v);
	fam.setExpectedEnd(22);

	EXPECT_TRUE(repeatExpects(&fam, &v));
	EXPECT_FALSE(repeatExpects(&fam, &u));
	u.push_back(21);
	EXPECT_TRUE(repeatExpects(&fam, &u));
	u.push_back(22);
	EXPECT_FALSE(repeatExpects(&fam, &u));
}

TEST(Family, setLast) {
	Family fam;

	LmerVector v(0);
	v.push_back(10);
	fam.adopt(&v);

	EXPECT_EQ(fam.getLast(), &v);
	EXPECT_EQ(fam.getLastIndex(), 10);

	LmerVector u(1);
	fam.adopt(&u);
	u.push_back(11);
	fam.setLast(&u);
	EXPECT_EQ(fam.getLast(), &u);
	EXPECT_EQ(fam.getLastIndex(), 11);
	u.push_back(20);
	EXPECT_EQ(fam.getLastIndex(), 11);
}


TEST(SeedChain, tieLooseEnds) {
  Family *fam = new Family;
  vector<Family*> families;
  families.push_back(fam);

  LmerVector v(0);
  LmerVector u(1);
  LmerVector w(2);

  fam->adopt(&v);
  fam->adopt(&u);
  fam->adopt(&w);

  v.push_back(10);
  u.push_back(11);
  w.push_back(12);

  fam->setLast(&w);

  tieLooseEnds(families);
  EXPECT_TRUE(families.size() == 1);

  fam->setLast(&v);

  tieLooseEnds(families);
  ASSERT_TRUE(families.size() == 2);

  EXPECT_TRUE(families[0]->getPrefix() == &v);
  EXPECT_TRUE(families[1]->getPrefix() == &u);

  delete families[0];
  delete families[1];
}



TEST(SeedChain, getElementaryFamilies1) {
	seqan::Dna5String sequence = "TAAACTAGGTCACTGTAAACTTGGTCACT";

	vector<seqan::CharString> masks;
	masks.push_back("111");

	vector<Family*> families;
	EXPECT_NO_THROW(getElementaryFamilies(sequence, masks, families));

	ASSERT_EQ(families.size(), 3);

	//ACT
	Family *fam = families[0];
	ASSERT_EQ(fam->size(), 4);
	EXPECT_EQ(fam->getLmers()->size(), 1);
	EXPECT_EQ(fam->getLmers()->front()->front(), 3);

	//TAAAC
	fam = families[1];
	ASSERT_EQ(fam->size(), 2);
	EXPECT_EQ(fam->getLmers()->size(), 3);
	EXPECT_EQ(fam->getLmers()->front()->front(), 0);

	//GGTCAC
	fam = families[2];
	ASSERT_EQ(fam->size(), 2);
	EXPECT_EQ(fam->getLmers()->size(), 4);
	EXPECT_EQ(fam->getLmers()->front()->front(), 7);
}


TEST(SeedChain, getElementaryFamilies2) {
	seqan::Dna5String sequence = "CCACGTACTNACGTNCCACGTAANTACACGTA";

	vector<seqan::CharString> masks;
	masks.push_back("111");

	vector<Family*> families;
	EXPECT_NO_THROW(getElementaryFamilies(sequence, masks, families));

	ASSERT_EQ(families.size(), 5);

	//ACTT
	Family *fam = families[0];
	ASSERT_EQ(fam->size(), 4);
	EXPECT_EQ(fam->getLmers()->size(), 2);
	EXPECT_EQ(fam->getLmers()->front()->front(), 2);

	//CCA
	fam = families[1];
	ASSERT_EQ(fam->size(), 2);
	EXPECT_EQ(fam->getLmers()->size(), 1);
	EXPECT_EQ(fam->getLmers()->front()->front(), 0);

	//GTA
	fam = families[2];
	ASSERT_EQ(fam->size(), 3);
	EXPECT_EQ(fam->getLmers()->size(), 1);
	EXPECT_EQ(fam->getLmers()->front()->front(), 4);

	//TAC
	fam = families[3];
	ASSERT_EQ(fam->size(), 2);
	EXPECT_EQ(fam->getLmers()->size(), 1);
	EXPECT_EQ(fam->getLmers()->front()->front(), 5);

	//CAC
	fam = families[4];
	ASSERT_EQ(fam->size(), 3);
	EXPECT_EQ(fam->getLmers()->size(), 1);
	EXPECT_EQ(fam->getLmers()->front()->front(), 1);
}


TEST(SeedChain, getElementaryFamilies3) {
	seqan::Dna5String sequence = "ACGTACCTGNTACCTGNACGTACNTACCTGNACGTACCTG";

	vector<seqan::CharString> masks;
	masks.push_back("111");

	vector<Family*> families;
	EXPECT_NO_THROW(getElementaryFamilies(sequence, masks, families));

	EXPECT_EQ(families.size(), 3);

	//TAC
	Family *fam = families[0];
	ASSERT_EQ(fam->size(), 5);
	EXPECT_EQ(fam->getLmers()->size(), 1);
	EXPECT_EQ(fam->getLmers()->front()->front(), 3);

	//ACGTA
	fam = families[1];
	ASSERT_EQ(fam->size(), 3);
	EXPECT_EQ(fam->getLmers()->size(), 3);
	EXPECT_EQ(fam->getLmers()->front()->front(), 0);

	//ACCTG
	fam = families[2];
	ASSERT_EQ(fam->size(), 4);
	EXPECT_EQ(fam->getLmers()->size(), 3);
	EXPECT_EQ(fam->getLmers()->front()->front(), 4);
}
