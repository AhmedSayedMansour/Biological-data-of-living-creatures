#ifndef PROTEIN_H
#define PROTEIN_H

#include "Sequence .h"
#include "DNA.h"
#include "RNA.h"

#include <algorithm>
#include <iostream>
#include <cctype>
#include <string>
#include <cstring>

using namespace std;

class DNA;
class RNA;
class Sequence;


enum Protein_Type {Hormon, Enzyme, TF, Cellular_Function};

inline ostream& operator<<(ostream &out,const Protein_Type& t)
{
    switch(t)      //finding the type of the protein
    {
    case Hormon:
        return out<<"(Hormon)";
    case Enzyme:
        return out<<"(Enzyme)";
    case TF:
        return out<<"(TF)";
    case Cellular_Function:
        return out<<"(Cellular_Function)";
    default:
        return out<<"(invalid value)";
    }
}

inline istream& operator>>( istream& in, Protein_Type& t )
{
    string text;
    while (in >> text)     //check if the type is right
    {
        if (text == "Hormon"){
            t = Hormon;
            break;
        }
        else if (text == "Enzyme"){
            t = Enzyme;
            break;
        }
        else if (text == "TF"){
            t = TF;
            break;
        }
        else if (text == "Cellular_Function"){
            t = Cellular_Function;
            break;
        }
        else{
            cout<<"Type Not Found !";
        }
    }
    return in;
}

class Protein : public Sequence
{
private:
    Protein_Type type;
    int startIndx;
    int endIndx;
    char *myProtein_Seq2;
public:
    // constructors and destructor
    Protein();
    Protein(char * p,Protein_Type atype);
    Protein(const Protein& rhs);
    ~Protein();

    void setIndx(int indx1,int indx2);

    void Print();   //to print the inf of the

    char* getSeq()
    {
        return seq;
    }
    void setSeq(char* seqq)
    {
        seq=seqq;
    }

    // return an array of DNA sequences that can possibly
    // generate that protein sequence
    DNA GetDNAStrandsEncodingMe( DNA & bigDNA);

    bool operator==(const Protein& p1);
    bool operator!=(const Protein& p1);
    Protein operator+(Protein& p1);

    friend ostream& operator<<(ostream& out,const Protein& p1)
    {
        cout<<"Protein sequence "<<p1.type<<" has 1 strands\nstrand1: "<<
            p1.seq<<"\t index "<<p1.startIndx<<" : "<<p1.endIndx;
        return out;
    }
    friend istream& operator>>(istream& in,Protein& p1)
    {
        string temp;
        cout<<"Enter Protein sequence: ";
        in>>temp;
        transform(temp.begin(), temp.end(), temp.begin(), ::toupper);
        char *myProtein_Seq = new char[temp.length() + 1];
        strcpy(myProtein_Seq, temp.c_str());
        p1.seq=myProtein_Seq;
        cout<<"Enter Protein Start Index: ";
        in>>p1.startIndx;
        cout<<"Enter Protein End Index: ";
        in>>p1.endIndx;
        cout<<"Enter Protein type: ";
        in>>p1.type;
        try
        {
            for(int i=0; i<temp.length(); i++)
            {
                if(temp[i]!='K'&&temp[i]!='N'&&temp[i]!='T'&&temp[i]!='R'&&temp[i]!='S'&&temp[i]!='I'&&
                   temp[i]!='M'&&temp[i]!='I'&&temp[i]!='Q'&&temp[i]!='H'&&temp[i]!='S'&&temp[i]!='P'&&
                   temp[i]!='E'&&temp[i]!='D'&&temp[i]!='A'&&temp[i]!='G'&&temp[i]!='V'&&temp[i]!='Y'&&
                   temp[i]!='S'&&temp[i]!='C'&&temp[i]!='W'&&temp[i]!='L'&&temp[i]!='F')
                {
                    throw 'E';       //checking the error and throwing character if exist
                }
            }
            if(p1.startIndx<0 || p1.endIndx<0)
            {
                throw 0;
            }
        }
        catch(char ch)
        {
            cout<<"\nInvalid Characters for Protein !\n";
            in>>temp;
            transform(temp.begin(), temp.end(), temp.begin(), ::toupper);
            char *myProtein_Seq = new char[temp.length() + 1];
            strcpy(myProtein_Seq, temp.c_str());
            p1.seq=myProtein_Seq;
        }
        catch(int n)
        {
            cout<<"\nInvalid Indexes !\n";
            cout<<"Enter Protein Start Index: ";
            in>>p1.startIndx;
            cout<<"Enter Protein End Index: ";
            in>>p1.endIndx;
        }
        return in;
    }
    char* Align(Sequence & s1, Sequence & s2);  //LCS Alignment

    //files handling
    void SaveSequenceToFile(char* filename);
    void LoadSequenceFromFile(char* filename);

};
#endif // PROTEIN_H
