#ifndef DNA_H
#define DNA_H

#include "Sequence .h"
#include "RNA.h"
///#include "Protein.h"

#include <algorithm>
#include <iostream>
#include <cctype>
#include <string>
#include <cstring>
using namespace std;

enum DNA_Type {promoter, motif, tail, noncoding};       //enum for the four DNA types

inline ostream& operator<<(ostream &out,const DNA_Type& t)  //oustream operator overloaded for the enum
{
    switch(t)
    {
    case promoter:
        return out<<"(promoter)";
    case motif:
        return out<<"(motif)";
    case tail:
        return out<<"(tail)";
    case noncoding:
        return out<<"(noncoding)";
    default:
        return out<<"(invalid value)";
    }
}

inline istream& operator>>( istream& in, DNA_Type& t )     //istream operator overloaded for the enum
{
    string text;
    while (in >> text)
    {
        if (text == "promoter"){
            t = promoter;
            break;
        }
        else if (text == "motif"){
            t = motif;
            break;
        }
        else if (text == "tail"){
            t = tail;
            break;
        }
        else if (text == "noncoding"){
            t = noncoding;
            break;
        }
        else{
            cout<<"ERROR : Type Not Found !"<<endl;       //check if the type is wrong
        }
    }
    return in;
}

class DNA : public Sequence
{
private:
    DNA_Type type;                //enum type
    DNA * complementary_strand;
    int startIndex;
    int endIndex;
    char *myDNA_Seq;
public:
    // constructors and destructor
    DNA();
    DNA(char* seq1, DNA_Type atype);
    DNA(const DNA& rhs);
    ~DNA();

    void setIndex(int indx1,int indx2);
    void setIndexComplement(int indx1,int indx2);
    void setIndexComplement();
    bool check_ACGT();

    // function to build the second strand/pair of DNA sequence
    // To build a complementary_strand (starting from the startIndex to
    // the endIndex), convert each A to T, each T to A, each C to G, and
    // each G to C. Then reverse the resulting sequence.
    void BuildComplementaryStrand();

    // function printing DNA sequence information to user
    void Print();

    char* getSeq()
    {
        return seq;
    }
    void setSeq(char* seqq)
    {
        seq=seqq;
    }
    // function to convert the DNA sequence to RNA sequence
    // It starts by building the complementary_strand of the current
    // DNA sequence (starting from the startIndex to the endIndex), then,
    // it builds the RNA corresponding to that complementary_strand.

    RNA& ConvertToRNA();

    bool operator==(const DNA& d1);
    bool operator!=(const DNA& d1);
    DNA operator+(DNA& d1);

    friend ostream& operator<<(ostream& out,const DNA& dna)
    {
        cout<<"DNA sequence "<<dna.type<<" has 1 strand\nstrand: "<<
            dna.seq<<"\t index "<<dna.startIndex<<" : "<<dna.endIndex;
        return out;
    }

    friend istream& operator>>(istream& in,DNA& dna)            //istream operator overloaded to cin the DNA
    {
        string temp;
        cout<<"Enter DNA sequence: ";
        in>>temp;
        transform(temp.begin(), temp.end(), temp.begin(), ::toupper);   //transforming all the seq to UPPER
        char *mySeq = new char[temp.length() + 1];     //
        strcpy(mySeq, temp.c_str());                   //  transform from string to array of char
        dna.seq=mySeq;                                 //
        cout<<"Enter DNA Start Index: ";
        in>>dna.startIndex;
        cout<<"Enter DNA End Index: ";
        in>>dna.endIndex;
        cout<<"Enter DNA type: ";
        in>>dna.type;
        //Exception Handling
        try{
            for(int i=0;i<temp.length();i++){
                if(temp[i]!='A'&&temp[i]!='C'&&temp[i]!='G'&&temp[i]!='T'){
                    throw 'E';     // throw character if wrong
                }
            }
            if(dna.startIndex<0 || dna.endIndex<0){
                    throw 0;        //throw integer if wrong
            }
        }
        catch(char ch){
            cout<<"\nInvalid Characters for DNA !\n";
            cout<<"Enter DNA sequence: ";
            in>>temp;
            transform(temp.begin(), temp.end(), temp.begin(), ::toupper);   //transforming all the seq into a UPPER
            char *mySeq = new char[temp.length() + 1];
            strcpy(mySeq, temp.c_str());
            dna.seq=mySeq;
        }
        catch(int n){
            cout<<"\nInvalid Indexes !\n";
            cout<<"Enter DNA Start Index: ";
            in>>dna.startIndex;
            cout<<"Enter DNA End Index: ";
            in>>dna.endIndex;
        }
        return in;
    }
    char* Align(Sequence & s1, Sequence & s2);                  //LCS Alignment
    char* LocalAlign(Sequence&  s1, Sequence&  s2);             //Local Alignment

    //Files Handling
    void SaveSequenceToFile(char* filename);
    void LoadSequenceFromFile(char* filename);
};
#endif // DNA_H
