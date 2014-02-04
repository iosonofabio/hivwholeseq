// vim:
/*
 *  Created on:	26/01/2014
 *  Author:	Fabio Zanini
 *  Contents:	Main header for the python wrapper of the SeqAn library.
 */
// Includes
#include <string>
#include  <iostream>
#include <sstream>

// Functions
int nothing();
int find_seed(std::string refseq, std::string seed, int *scoreout);
