#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main(int argc, char *argv[]){
    typedef String <char> TSequence;
    typedef Align<TSequence,ArrayGaps> TAlign;

    typedef StringSet<TSequence> TStringSet;       // container for strings
    typedef StringSet<TSequence, Dependent<> > TDepStringSet;   // dependent string set
    typedef Graph<Alignment<TDepStringSet> > TAlignGraph;   // alignment graph

    std::cout << "N of args: " << argc << std::endl;
    if(argc != 3) {
	    std::cerr << "Wrong arguments!";
	    std::cerr << "Example call: " << argv[0] << " ACCT ACAT" << std::endl;
	    return 1;
    }

    TSequence seq1 = argv[1];
    TSequence seq2 = argv[2];

    // Global
    std::cout << "GLOBAL" << std::endl;
    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0),seq1);
    assignSource(row(align,1),seq2);

    int score = globalAlignment(align, Score<int,Simple>(0,-1,-1),-100,100);

    std::cout << "Score: " << score << std::endl;
    std::cout << align << std::endl;

    // Overlap
    std::cout << "OVERLAP" << std::endl;
    TStringSet align2;
    appendValue(align2,seq1);
    appendValue(align2,seq2);
    TAlignGraph alignG(align2);

    score = globalAlignment(alignG,
		    Score<int,Simple>(3,-3,-2,-2),
		    AlignConfig<true, false, false, true>(),
		    -100, 100);
    std::cout << "Score: " << score << std::endl;
    std::cout << alignG << std::endl;

    // Local
    std::cout << "LOCAL" << std::endl;
    Align< String<char> > ali;
    resize(rows(ali), 2);
    assignSource(row(ali, 0), argv[1]);
    assignSource(row(ali, 1), argv[2]);

    score = localAlignment(ali, Score<int>(3,-3,-2, -2));
    std::cout << "Score = " << score << std::endl;
    std::cout << ali;
    std::cout << "Aligns Seq1[" << clippedBeginPosition(row(ali, 0)) << ":" << (clippedEndPosition(row(ali, 0))-1) << "]";
    std::cout << " and Seq2[" << clippedBeginPosition(row(ali, 1)) << ":" <<  (clippedEndPosition(row(ali, 1))-1) << "]" << std::endl << std::endl;

    return 0;
}

