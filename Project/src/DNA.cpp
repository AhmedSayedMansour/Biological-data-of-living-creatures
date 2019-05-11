#include "DNA.h"
#include "RNA.h"
#include "Sequence .h"

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <cstring>
#include <sstream>

using namespace std;

DNA::DNA()
{
    startIndex=0;
    endIndex=0;
    type=promoter;       //enum type promoter
    ///complementary_strand=NULL;
}

DNA::DNA(char* seq, DNA_Type atype)
{
    this->seq=seq;
    type=atype;

}

//default constructor
DNA::~DNA()
{
    ///delete[] myDNA_Seq;

}

//copy constructor
DNA::DNA(const DNA& rhs)
{
    startIndex = rhs.startIndex;
    endIndex = rhs.endIndex;
    type = rhs.type;
    seq=rhs.seq;
}

//setting the indeces
void DNA::setIndex(int indx1,int indx2)
{
    startIndex=indx1;
    endIndex=indx2;
}

//checking if there exist a wrong DNA
bool DNA::check_ACGT()
{
    int steps=abs(startIndex-endIndex);
    for(int i=0; i<steps; i++)                     //for loop on all the sequence
    {
        switch(seq[i])
        {                            //checking if all are on all right
        case 'A':
            continue;
        case 'C':
            continue;
        case 'G':
            continue;
        case 'T':
            continue;
        default:
            return 0;
        }
    }
    return 1;
}

void DNA::setIndexComplement(int indx1,int indx2)    //setting the index of the complementary strand
{
    complementary_strand->startIndex=indx2;
    complementary_strand->endIndex=indx1;
}

void DNA::setIndexComplement()                      //complement the index of the complementary strand
{
    complementary_strand->startIndex=endIndex;
    complementary_strand->endIndex=startIndex;
}

void DNA::BuildComplementaryStrand()
{
    int steps=abs(startIndex-endIndex);    //the difference between start and end indexes
    string myComp;

    for(int i=0; i<steps; i++)           //building the complementary strand with a for loop on all the strand
    {
        if(seq[i]=='T')                  //storing the complementary strand in string
            myComp+='A';
        else if(seq[i]=='A')
            myComp+='T';
        else if(seq[i]=='C')
            myComp+='G';
        else if(seq[i]=='G')
            myComp+='C';
    }
    char* myDNA_Seq = new char[myComp.length() + 1];     //creating char with a size of the complemented string
    strcpy(myDNA_Seq, myComp.c_str());                   //coping the string into the character variable
    complementary_strand->seq=myDNA_Seq;                 //
///    delete[]myDNA_Seq;
}

//printing inf of the DNA
void DNA::Print()
{
    cout<<"DNA sequence "<<type<<" has 2 strands\nstrand1: "<<seq
        <<"\t index "<<startIndex<<" : "<<endIndex<<"\nstrand2: "
        << complementary_strand->seq<<"\t index "
        << complementary_strand->startIndex<<" : "
        << complementary_strand->endIndex;
}

RNA& DNA::ConvertToRNA()
{
    string myComp;
    int steps=abs(startIndex-endIndex);

    for(int i=0; i<steps; i++)
    {
        if(seq[i]=='T')
            myComp+='A';
        else if(seq[i]=='A')
            myComp+='U';
        else if(seq[i]=='C')
            myComp+='G';
        else if(seq[i]=='G')
            myComp+='C';
    }
    char* myDNA_Seq = new char[myComp.length() + 1];
    strcpy(myDNA_Seq, myComp.c_str());

    static RNA rna(myDNA_Seq,mRNA);
    rna.setIndx(startIndex,endIndex);
    return rna;
}

bool DNA::operator==(const DNA& d1)       //to check if 2 seq are equal
{
    int steps=abs(startIndex-endIndex);
    if(steps != abs(d1.startIndex-d1.endIndex))
        return 0;

    for(int i=0; i<steps; i++)
    {
        if(seq[i]!=d1.seq[i])
            return 0;
    }
    return 1;
}

bool DNA::operator!=(const DNA& d1)      //to check if 2 seq are not equal
{
    int steps=abs(startIndex-endIndex);
    if(steps != abs(d1.startIndex-d1.endIndex))
        return 1;

    for(int i=0; i<steps; i++)
    {
        cout<<seq[i]<<d1.seq[i];
        if(seq[i]!=d1.seq[i])
            return 1;
    }
    return 0;
}

DNA DNA::operator+(DNA& d1)    //merge
{
    int steps=abs(startIndex-endIndex);
    string temp;                            //to store a new value
    for(int i=0; i<steps; i++)
    {
        temp+=seq[i];
    }
    int steps2=abs(d1.startIndex-d1.endIndex);
    for(int i=0; i<steps2; i++)
    {
        temp+=d1.seq[i];
    }

    myDNA_Seq = new char[temp.length() + 1];
    strcpy(myDNA_Seq, temp.c_str());            //convert from string to character

    DNA sum(myDNA_Seq,type);
    sum.setIndex(startIndex,d1.endIndex);
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


char* DNA::Align(Sequence & s1, Sequence & s2)
{

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

    int index = L[m][n]; //the index of the box for ->  backtracking  <--
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
    for(int i=0; i<l; ++i)
    {
        temp+=lcs[i];
    }

    char* myDNA_Seq3 = new char[temp.length() + 1];
    strcpy(myDNA_Seq3, temp.c_str());

    seq =myDNA_Seq3;
    return myDNA_Seq3;
}


//Local Alignment
//before all here are some algorithms mot fully completed
int align( string a,  string b, int alpha_gap,
           int alpha[26][26], string &a_aligned,
           string &b_aligned)
{
    int  n = a.size();   //size of the first seq
    int  m = b.size();   //size of the first seq

    vector<vector<int> > A(n + 1, vector<int>(m + 1));   //the map of the local alignment

    for (int i = 0; i <= m; ++i)
        A[0][i] = alpha_gap * i;
    for (int i = 0; i <= n; ++i)
        A[i][0] = alpha_gap * i;

    for (int i = 1; i <= n; ++i)
    {
        for (int  j = 1; j <= m; ++j)
        {
            char x_i = a[i-1];
            char y_j = b[j-1];
            A[i][j] = min(A[i-1][j-1] + alpha[x_i - 'A'][y_j - 'A'],min(A[i-1][j] + alpha_gap,
                          A[i][j-1] + alpha_gap));

        }
    }

    // print2DVector(A);

    a_aligned = "";
    b_aligned = "";
    int j = m;
    int i = n;
    for (; i >= 1 && j >= 1; --i)
    {
        char x_i = a[i-1];
        char y_j = b[j-1];
        if (A[i][j] == A[i-1][j-1] + alpha[x_i - 'A'][y_j - 'A'])
        {

///           * I think prepending chars this way to a std::string is very inefficient.
///           * Is there any better way of doing this without using C-style strings?

            a_aligned = x_i + a_aligned;
            b_aligned = y_j + b_aligned;
            --j;
        }
        else if (A[i][j] == A[i-1][j] + alpha_gap)
        {
            a_aligned = x_i + a_aligned;
            b_aligned = '-' + b_aligned;
        }
        else
        {
            a_aligned = '-' + a_aligned;
            b_aligned = y_j + b_aligned;
            --j;
        }
    }

    while (i >= 1 && j < 1)
    {
        a_aligned = a[i-1] + a_aligned;
        b_aligned = '-' + b_aligned;
        --i;
    }
    while (j >= 1 && i < 1)
    {
        a_aligned = '-' + a_aligned;
        b_aligned = b[j-1] + b_aligned;
        --j;
    }

    return A[n][m];
}

char* DNA::LocalAlign(Sequence&  s1, Sequence&  s2)
{
    char *X=s1.getSeq();
    char *Y=s2.getSeq();


    // Penalty for any alphabet matched with a gap
    int gap_penalty = 0;

    /*
     * alpha[i][j] = penalty for matching the ith alphabet with the
     *               jth alphabet.
     * Here: Penalty for matching an alphabet with another one is 1
     *       Penalty for matching an alphabet with itself is 0
     */
    int alpha[26][26];
    for (size_t i = 0; i < 26; ++i)
    {
        for (size_t j = 0; j < 26; ++j)
        {
            if (i == j)
                alpha[i][j] = 0;
            else
                alpha[i][j] = 1;
        }
    }

    // Aligned sequences
    //string a1,a2;
    char* ss1=s1.getSeq();
    char* ss2=s2.getSeq();

    //char* sseq;
    //sseq=bigDNA.getSeq();
    std::string a11(ss1);
    std::string a22(ss2);

//    char* myDNA_Seq4 = new char[temp.length() + 1];
//    strcpy(myDNA_Seq4, temp.c_str());
//
//    char* myDNA_Seq4 = new char[temp.length() + 1];
//    strcpy(myDNA_Seq4, temp.c_str());


    string a2,b2;
    int penalty = align(a11, a22, gap_penalty, alpha, a2, b2);

    cout << "seq 1: " << a11 << endl;
    cout << "seq 2: " << a22 << endl;
    cout << "Needleman-Wunsch Score: " << penalty << endl;
    cout << "Aligned sequences: " << endl;
    cout << a2 << endl;
    cout << b2 << endl;

}

void DNA::SaveSequenceToFile(char* filename)
{
    string text,temp;
    fstream myfile ( filename, ios :: out | ios::app);
    if( myfile.fail() )
    {
        cout<<"The file didn't opened !"<<endl;
        return;
    }
    cout<<"Please Enter Your DNA here: ";
    cin.ignore();
    getline(cin,text);
    myfile << "DNA\n";
    myfile << text <<"\n";
    myfile.close();

}

void DNA::LoadSequenceFromFile(char* filename)
{
    string text,temp;
    fstream myfile ( filename, ios :: in);
    if( myfile.fail() )
    {
        cout<<"The file didn't opened !"<<endl;
        return;
    }
    while (myfile>>temp)
    {
        if(temp=="DNA")
        {
            myfile.ignore();
            getline(myfile,text);
            cout<<text<<endl;
        }
    }
    myfile.close();
}
