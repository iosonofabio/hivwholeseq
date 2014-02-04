// vim:
/*
 *  Created on:	30 Jan 2014
 *  Author:	Fabio Zanini
 *  Contents:	Look for seed sequence in a reference
 */
#include "seqanpy.h"

int find_seed(std::string refseq, std::string seed, int *scoreout) {

	int score = 0;
	int score_old = -1;
	int pos_seed = -1;

	// this is not like Python, but it should not matter (a little bit less optimized)
	int rl = seed.size();

	// set the end
	std::string::iterator end = refseq.end();
	for(int i=0; i < rl; i++, end--);

	bool found = false;
	int i = 0;
	for (std::string::iterator it=refseq.begin(); (it!=end) and (!found); i++, it++) {
		score = 0;
		for(std::string::iterator itl=it, its=seed.begin(); its != seed.end();
				itl++, its++) {
			if(*itl == *its)
				score ++;
		}
		if(score == rl) {
			pos_seed = i;
			found = true;
		} else if(score > score_old) {
			pos_seed = i;
			score_old = score;
		}
	}

	*scoreout = score;

	return pos_seed;
}
