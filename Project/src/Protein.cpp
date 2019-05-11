#include "Protein.h"
#include "RNA.h"
#include "DNA.h"
#include "Sequence .h"

Protein::Protein()
{
    type=Cellular_Function;
    startIndx=0;
    endIndx=0;
}

Protein::Protein(char* p,Protein_Type atype)
{
    this->seq=p;
    type = atype;
}

Protein::~Protein()
{

}

Protein::Protein(const Protein& rhs)
{
    startIndx = rhs.startIndx;
    endIndx = rhs.endIndx;
    type = rhs.type;
    seq=rhs.seq;
}

void Protein::setIndx(int indx1,int indx2)
{
    startIndx=indx1;
    endIndx=indx2;
}


void Protein::Print()        //printing the inf of the protein
{
    cout<<"Protein sequence "<<type<<" has 1 strands\nstrand1: "<<
        seq<<"\t index "<<startIndx<<" : "<<endIndx;
}

bool Protein::operator==(const Protein& p1)
{
    int steps=abs(startIndx-endIndx);
    if(steps != abs(p1.startIndx-p1.endIndx))
        return 0;
    for(int i=0; i<steps; i++)
    {
        if(seq[i]!=p1.seq[i])
            return 0;         //if there exist one char not equal
    }
    return 1;
}

bool Protein::operator!=(const Protein& p1)
{
    int steps=abs(startIndx-endIndx);
    if(steps != abs(p1.startIndx-p1.endIndx))
        return 1;
    for(int i=0; i<steps; i++)
    {
        if(seq[i]!=p1.seq[i])
            return 1;
    }
    return 0;
}

Protein Protein::operator+(Protein& p1)
{
    int steps=abs(startIndx-endIndx);
    string temp;
    for(int i=0; i<steps; i++)
    {
        temp+=seq[i];
    }
    int steps2=abs(p1.startIndx-p1.endIndx);
    for(int i=0; i<steps2; i++)
    {
        temp+=p1.seq[i];
    }

    myProtein_Seq2 = new char[temp.length() + 1];
    strcpy(myProtein_Seq2, temp.c_str());
    cout<<myProtein_Seq2<< endl;

    Protein sum;
    sum.seq=myProtein_Seq2;
    sum.setIndx(startIndx,p1.endIndx);

    return sum;
}
DNA Protein::GetDNAStrandsEncodingMe( DNA & bigDNA)
{
    char* sseq;
    sseq=bigDNA.getSeq();
    std::string DNA_Seq(sseq);
    std::string AminoAcid(seq);

    int pSize=DNA_Seq.length(), aSize=AminoAcid.length();
    string temp1, temp2;
    CodonsTable test;
    test.LoadCodonsFromFile( "codon.txt" );

    for(int i=0; i<pSize; i+=aSize*3)
    {
        temp1=DNA_Seq.substr(i, aSize*3);
        temp2=test.getAminoAcid(temp1);
        if (AminoAcid==temp2)
        {
            char *myTemp_Seq = new char[temp1.length() + 1];
            strcpy(myTemp_Seq, temp1.c_str());
            DNA d;
            d.setSeq(myTemp_Seq);
            return d;
        }
    }
}


char* Protein::Align(Sequence & s1, Sequence & s2){

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

    char* myProtein_Seq3 = new char[temp.length() + 1];
    strcpy(myProtein_Seq3, temp.c_str());

    seq =myProtein_Seq3;
    return myProtein_Seq3;
}


void Protein::SaveSequenceToFile(char* filename){
      string text;
      cout<<"Please Enter Your Protein here: ";
      cin.ignore();
      getline(cin ,text);        //set the name of the file
      fstream myfile ( filename, ios :: out | ios::app);
      if( myfile.fail() )
      {
          cout<<"The file didn't opened !"<<endl;     //check if the is already still exist
          return;
      }
      myfile << "Protein\n";
      myfile << text <<"\n";
      myfile.close();
}

void Protein::LoadSequenceFromFile(char* filename){
    string text,temp;
    fstream myfile ( filename, ios :: in);
    if( myfile.fail() )
    {
        cout<<"The file didn't opened !"<<endl;
        return;
    }
    while (myfile>>temp)
    {
        if(temp=="Protein")
        {
            myfile.ignore();
            getline(myfile,text);
            cout<<text<<endl;
        }
    }
    myfile.close();    //finally closing the file
}
