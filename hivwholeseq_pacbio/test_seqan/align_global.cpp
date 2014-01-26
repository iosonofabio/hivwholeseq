#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main(int argc, char *argv[]){
    typedef String <char> TSequence;
    typedef Align<TSequence,ArrayGaps> TAlign;

    typedef StringSet<TSequence> TStringSet;       // container for strings
    typedef StringSet<TSequence, Dependent<> > TDepStringSet;   // dependent string set
    typedef Graph<Alignment<TDepStringSet> > TAlignGraph;   // alignment graph

    int band = -1;
    int score;

    // Command line arguments
    if((argc != 3) and (argc != 4)) {
	    std::cerr << "Wrong arguments! Example calls: ";
	    std::cerr << argv[0] << " ACCT ACAT, ";
	    std::cerr << argv[0] << " ACCT ACAT 100 (for banded alignments).";
	    std::cerr << std::endl;
	    return 1;
    } else if(argc == 4)
	    band = atoi(argv[3]);

    TSequence seq1 = argv[1];
    TSequence seq2 = argv[2];

    //// Global alignment
    //TAlign align;
    //resize(rows(align), 2);
    //assignSource(row(align,0),seq1);
    //assignSource(row(align,1),seq2);

    //if(band >= 0)
    //    score = globalAlignment(align, Score<int,Simple>(0,-1,-1), -band, band);
    //else
    //    score = globalAlignment(align, Score<int,Simple>(0,-1,-1));

    //std::cout << "Score: " << score << std::endl;
    //std::cout << align << std::endl;

    
    // Overlap alignment
    TStringSet sequences;
    appendValue(sequences,seq1);
    appendValue(sequences,seq2);

    TAlignGraph alignG(sequences);
    if(band >= 0)
	    score = globalAlignment(alignG,
		    Score<int,Simple>(1,-1,-1),
		    AlignConfig<true, false, false, true>(),
		    -band, band);
    else
	    score = globalAlignment(alignG,
		    Score<int,Simple>(1,-1,-1),
		    AlignConfig<true, true, true, true>());
    std::cout << "Score: " << score << ::std::endl;
    std::cout << alignG << ::std::endl;


    return 0;
}

