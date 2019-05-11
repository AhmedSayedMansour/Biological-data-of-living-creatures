#include "DNA.h"
#include "RNA.h"
#include "CodonsTable.h"
#include "Sequence .h"

#include <iostream>
#include <string>
#include <cstring>
#include <sstream>

using namespace std;

//Default Constructor
RNA::RNA()
{
    startIndx=0;
    endIndx=0;
    type=mRNA;
}

RNA::RNA(char* seq, RNA_Type atype)
{
    this->seq=seq;
    type=atype;
}

RNA::~RNA()
{
    ///delete[] myRNA_Seq;
    ///delete[] myRNA_Seq2;
}

RNA::RNA(const RNA& rhs)
{
    startIndx = rhs.startIndx;
    endIndx = rhs.endIndx;
    type = rhs.type;
    seq=rhs.seq;
}

void RNA::Print()    //print the type and the RNA seq
{
    cout<<"RNA sequence "<<type<<" has 1 strand: "<<this->seq;
}

bool RNA::check_ACGU()         //tho check if the RNA seq is right
{
    int steps=abs(startIndx-endIndx);      //the size of the RNA seq
    for(int i=0; i<steps; i++)
    {
        switch(seq[i])
        {
        case 'A':
            continue;
        case 'C':
            continue;
        case 'G':
            continue;
        case 'U':
            continue;
        default:
            return 0;
        }
    }
    return 1;
}

void RNA::setIndx(int indx1,int indx2)
{
    startIndx=indx1;
    endIndx=indx2;
}

DNA& RNA::ConvertToDNA()
{
    int steps=abs(startIndx-endIndx);
    string temp;                       //considered as a container
    for(int i=0; i<steps; i++)
    {
        string t;
        if(seq[i]=='U')
        {
            temp+='T';
        }
        else{
            t=seq[i];
            temp+=t;
        }
    }
    char *myRNA_Seq = new char[temp.length() + 1];
    strcpy(myRNA_Seq, temp.c_str());

    static DNA dna(myRNA_Seq,tail);    //static object just to call it just once then it will be removed
    dna.setIndex(startIndx,endIndx);
    return dna;
}

Protein RNA::ConvertToProtein( CodonsTable & table)
    {

        string FileName;
        cout << "Enter file name to load codons table : ";
        cin >> FileName;
        table.LoadCodonsFromFile( FileName );
        ///cout << "Enter RNA : ";
        ///cin >> rna;
        int steps=abs(startIndx-endIndx);
    //	cout << "AminoAcids are : " <<
    //  table.getAminoAcid(seq ,steps) << endl;
        string t=table.getAminoAcid(seq,steps);

        char *myProtein_Seq = new char[t.length() + 1];
        strcpy(myProtein_Seq, t.c_str());
        Protein p1(myProtein_Seq,Cellular_Function);
        return p1;
    //	cout << "Enter proteins and AminoAcids to find all compinations : " << endl;
    //	cin >> protein;
    //	cin >> AminoAcid;
    //	cout << "All compinations are : " << endl ;
    //	table.proteinsCompinations(protein, AminoAcid);

}

bool RNA::operator==(const RNA& r1)        //to check the equality
{
    int steps=abs(startIndx-endIndx);
    if(steps != abs(r1.startIndx-r1.endIndx))
        return 0;
    for(int i=0; i<steps; i++)
    {
        if(seq[i]!=r1.seq[i])
            return 0;
    }
    return 1;                   //return 1 if not equal
}

bool RNA::operator!=(const RNA& r1)
{
    int steps=abs(startIndx-endIndx);
    if(steps != abs(r1.startIndx-r1.endIndx))
        return 1;
    for(int i=0; i<steps; i++)
    {
        if(seq[i]!=r1.seq[i])
            return 1;
    }
    return 0;
}

RNA RNA::operator+(RNA& r1)  //merge
{
    int steps=abs(startIndx-endIndx);
    string temp;
    for(int i=0; i<steps; i++)
    {
        temp+=seq[i];
    }
    int steps2=abs(r1.startIndx-r1.endIndx);
    for(int i=0; i<steps2; i++)
    {
        temp+=r1.seq[i];
    }
    char *myRNA_Seq2 = new char[temp.length() + 1];
    strcpy(myRNA_Seq2, temp.c_str());

    RNA sum(myRNA_Seq2,type);
    sum.setIndx(startIndx,r1.endIndx);
    return sum;
}

///LCS Alignment  :::

/// The pseudo-code for the algorithm to compute the F matrix of the LCS therefore looks like this:

/*
d <- MismatchScore
for i=0 to length(A)
  F(i,0) <- d*i
for j=0 to length(B)
  F(0,j) <- d*j
for i=1 to length(A)
  for j=1 to length(B)
  {
    Match <- F(i-1,j-1) + S(Ai, Bj)
    Delete <- F(i-1, j) + d
    Insert <- F(i, j-1) + d
    F(i,j) <- max(Match, Insert, Delete)
  }

  */

/// Internal algorithm of the LCS

/*
AlignmentA <- ""
AlignmentB <- ""
i <- length(A)
j <- length(B)
while (i > 0 or j > 0)
{
 if (i > 0 and j > 0 and F(i,j) == F(i-1,j-1) + S(Ai, Bj))
 {
   AlignmentA <- Ai + AlignmentA
   AlignmentB <- Bj + AlignmentB
   i <- i - 1
   j <- j - 1
 }
 else if (i > 0 and F(i,j) == F(i-1,j) + d)
 {
   AlignmentA <- Ai + AlignmentA
   AlignmentB <- "-" + AlignmentB
   i <- i - 1
 }
 else
 {
   AlignmentA <- "-" + AlignmentA
   AlignmentB <- Bj + AlignmentB
   j <- j - 1
 }
}
*/


char* RNA::Align(Sequence & s1, Sequence & s2){

    char *X=s1.getSeq();
    char *Y=s2.getSeq();
    int m = strlen(X);
    int n = strlen(Y);
    int L[m+1][n+1];
    for (int i=0; i<=m; i++)
    {
        for (int j=0; j<=n; j++)
        {
            if (i == 0 || j == 0)
                L[i][j] = 0;
            else if (X[i-1] == Y[j-1])
                L[i][j] = L[i-1][j-1] + 1;
            else
                L[i][j] = max(L[i-1][j], L[i][j-1]);
        }
    }

    int index = L[m][n]; //the index of the box for backtracking
    int l = L[m][n];  //the value of the most right bottom cell

    // Create a character array to store the LCS string
    char lcs[index+1];
    lcs[index] = ' '; // Set the terminating character by space

    // Start from the right-most-bottom-most corner and
    // one by one store characters in LCS[]
    int i = m, j = n;
    while (i > 0 && j > 0)
    {
        // If current character in X[] and Y are same, then
        // current character is part of LCS
        if (X[i-1] == Y[j-1])
        {
            lcs[index-1] = X[i-1]; // Put current character in result
            i--;
            j--;
            index--;     // reduce values of i, j and index
        }

        // If not same, then find the larger of two and
        // go in the direction of larger value
        else if (L[i-1][j] > L[i][j-1])
            i--;
        else
            j--;
    }

    string temp="";
    for(int i=0; i<l; ++i){
        temp+=lcs[i];
    }

    char* myRNA_Seq3 = new char[temp.length() + 1];
    strcpy(myRNA_Seq3, temp.c_str());

    seq =myRNA_Seq3;
    return myRNA_Seq3;
}


void RNA::SaveSequenceToFile(char* filename){
      string text;
      cout<<"Please Enter Your RNA here: ";
      cin.ignore();
      getline(cin ,text);
      fstream myfile ( filename, ios :: out | ios::app);      //opening the file
      if( myfile.fail() )
      {
         cout<<"The file didn't opened !"<<endl;
         return;
      }
      myfile << "RNA\n";
      myfile << text <<"\n";         //printing in the file
      myfile.close();
}
void RNA::LoadSequenceFromFile(char* filename){
    string text,temp;
    fstream myfile ( filename, ios :: in);
    if( myfile.fail() )
    {
        cout<<"The file didn't opened !"<<endl;           //Case of having error in opening the file
        return;
    }
    while (myfile>>temp)
    {

        if(temp=="RNA")
        {
            myfile.ignore();
            getline(myfile,text);                      //Reading from the file
            cout<<text<<endl;
        }
    }
    myfile.close();
}
