#ifndef RNA_H
#define RNA_H

#include "Sequence .h"
#include "DNA.h"
#include "Protein.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cctype>
#include <string>
#include <cstring>

using namespace std;



class DNA;
class Protein;

struct Codon
{
    char value[4];
    char AminoAcid;
};

class CodonsTable
{
private:
    Codon codons[64];
    //fstream file;
public:
    CodonsTable()
    {

    }

    ~CodonsTable()
    {

    }
    //to read all the codon from the file
    void LoadCodonsFromFile(string codonsFileName)
    {
        ifstream file(codonsFileName);
        int i=0 ;
        while ( file && !file.eof() )
        {
            file >> codons[i].value;
            file >> codons[i++].AminoAcid;
        }
        file.close();
    }
    //to change the codon's equality
    void setCodon( char AminoAcid2, int index){
        for(int i=0;i<64;++i){
            if(i==index){
                codons[i].AminoAcid=AminoAcid2;
                cout<<codons[i].value<<" : "<<codons[i].AminoAcid;
                break;
            }
        }
    }

    string getAminoAcid(char* valuee,int steps)
    {
        vector<string> RNA_Seq;
        string temp, out = "";

        for ( int i = 0 ; i <= steps-3 ; i+=3 )     //temp-3 because there are 3 non-equaled
        {
            temp="";
            temp += valuee[i];
            temp += valuee[i+1];
            temp += valuee[i+2];
            RNA_Seq.push_back( temp );
        }
        for ( int i = 0 ; i < RNA_Seq.size() ; i++ )
        {
            for ( int j = 0; j < 64 ; j++ )
            {
                if ( RNA_Seq[i] == codons[j].value )
                {
                    out += codons[j].AminoAcid;
                    break;
                }
                else if ( j == 60 )
                    out += ".";
            }
        }
        return out;
    }

    string getAminoAcid(string val)
    {
        vector<string> RNA_Seq;
        string temp, out = "";
        for ( int i = 0 ; i <= val.length()-3 ; i+=3 )
        {
            temp = val.substr(i,3);      //the same as what was in the last function
            RNA_Seq.push_back( temp );
        }
        for ( int i = 0 ; i < RNA_Seq.size() ; i++ )
        {
            for ( int j = 0; j < 64 ; j++ )
            {
                if ( RNA_Seq[i] == codons[j].value )
                {
                    out += codons[j].AminoAcid;
                    break;
                }
                else if ( j == 60 )
                    out += ".";
            }
        }
        return out;
    }
};


enum RNA_Type {mRNA, pre_mRNA, mRNA_exon, mRNA_intron};      //enum types of the RNA

inline ostream& operator<<(ostream &out,const RNA_Type& t)    //inline function th avoid the Error in the cpp.file
{
    switch(t)    //check if the type is right
    {
    case mRNA:
        return out<<"(mRNA)";
    case pre_mRNA:
        return out<<"(pre_mRNA)";
    case mRNA_exon:
        return out<<"(mRNA_exon)";
    case mRNA_intron:
        return out<<"(mRNA_intron)";
    default:
        return out<<"(invalid value)";
    }
}

inline istream& operator>>( istream& in, RNA_Type& t )
{
    string text;
    while (in >> text)     //check if the type is right
    {
        if (text == "mRNA"){
            t = mRNA;
            break;
        }
        else if (text == "pre_mRNA"){
            t = pre_mRNA;
            break;
        }
        else if (text == "mRNA_exon"){
            t = mRNA_exon;
            break;
        }
        else if (text == "mRNA_intron"){
            t = mRNA_intron;
            break;
        }
        else{
            cout<<"ERROR : Type Not Found !";
        }
    }
    return in;
}

class RNA: public Sequence
{
private:
    RNA_Type type;
    int startIndx;
    int endIndx;
    char *myRNA_Seq2;
public:
    // constructors and destructor
    RNA();
    RNA(char* seq, RNA_Type atype);
    RNA(const RNA& rhs);
    ~RNA();
    // function to be overridden to print all the RNA information
    void Print();
    void setIndx(int indx1,int indx2);

    bool check_ACGU();

    char* getSeq()
    {
        return seq;
    }
    void setSeq(char* seqq)
    {
        seq=seqq;
    }


    // function to convert the RNA sequence into protein sequence
    // using the codonsTable object

    Protein ConvertToProtein( CodonsTable & table);


    // function to convert the RNA sequence back to DNA
    DNA& ConvertToDNA();

    bool operator==(const RNA& r1);
    bool operator!=(const RNA& r1);
    RNA operator+(RNA& r1);

    friend ostream& operator<<(ostream& out,const RNA& rna)
    {
        cout<<"RNA sequence "<<rna.type<<" has 1 strands\nstrand1: "<<
            rna.seq<<"\t index "<<rna.startIndx<<" : "<<rna.endIndx;
        return out;
    }

    friend istream& operator>>(istream& in,RNA& rna)    //oustream operator overloaded
    {
        string temp;
        cout<<"Enter RNA sequence: ";
        in>>temp;
        transform(temp.begin(), temp.end(), temp.begin(), ::toupper);
        char *myRNA_Seq = new char[temp.length() + 1];
        strcpy(myRNA_Seq, temp.c_str());
        rna.seq=myRNA_Seq;
        cout<<"Enter RNA Start Index: ";
        in>>rna.startIndx;
        cout<<"Enter RNA End Index: ";
        in>>rna.endIndx;
        cout<<"Enter RNA type: ";
        in>>rna.type;
        //Exception Handling
        try
        {
            for(int i=0; i<temp.length(); i++)
            {
                if(temp[i]!='A'&&temp[i]!='C'&&temp[i]!='G'&&temp[i]!='U')
                {
                    throw 'E';
                }
            }
            if(rna.startIndx<0 || rna.endIndx<0)
            {
                throw 0;
            }
        }
        catch(char ch)
        {
            //retake the RNA
            cout<<"\nInvalid Characters for RNA !\n";
            in>>temp;
            transform(temp.begin(), temp.end(), temp.begin(), ::toupper);
            char *myRNA_Seq = new char[temp.length() + 1];
            strcpy(myRNA_Seq, temp.c_str());
            rna.seq=myRNA_Seq;
        }
        catch(int n)
        {
            //Retake the indexes
            cout<<"\nInvalid Indexes !\n";
            cout<<"Enter RNA Start Index: ";
            in>>rna.startIndx;
            cout<<"Enter RNA End Index: ";
            in>>rna.endIndx;
        }
        return in;
    }

    char* Align(Sequence & s1, Sequence & s2);   //LCS Alignment

    //Files Handling
    void SaveSequenceToFile(char* filename);
    void LoadSequenceFromFile(char* filename);

};

#endif // RNA_H
