#ifndef SEQUENCE _H
#define SEQUENCE _H

#include <iostream>
using namespace std;

class Sequence
{
protected:
    char* seq;
public:
    // constructors and destructor
    Sequence();
    Sequence(int length);
    Sequence(const Sequence& rhs);
    ~Sequence();
    // pure virtual function that should be overridden because every
    // type of sequence has its own details to be printed

    virtual char* getSeq()= 0;
    virtual void setSeq(char*)= 0;
    virtual void Print()= 0;
    virtual void SaveSequenceToFile(char* filename)=0;
    virtual void LoadSequenceFromFile(char* filename)=0;
    // friend function that will find the LCS (longest common
    // subsequence) between 2 sequences of any type, according to
    // polymorphism
    friend  char* Align(Sequence&  s1, Sequence&  s2);
    friend  char* LocalAlign(Sequence&  s1, Sequence&  s2);
   // void lcs( DNA *dna1, DNA *dna2 );
};

#endif // SEQUENCE _H
