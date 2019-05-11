
/****************************************************************/

/// FCI – Programming 2 – 2018 - Assignment 4
/// Program Name:               project 4.cpp
/// Last Modification Date:     13/12/2018
/// Author1 and ID and Group:   Ahmed Sayed     20170022
/// Author2 and ID and Group:   Eslam Saleh     20170046
/// Author3 and ID and Group:   Eslam Mohammed  20170049
/// Teaching Assistant:         Esraa Salem
/// Groups:                     G1 & G2
/// Version:                    2.4

/****************************************************************/

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "Sequence .h"
#include "CodonsTable.h"
#include "DNA.h"
#include "RNA.h"
#include "Protein.h"

using namespace std;

int main()
{

    /// Menu:-
    cout<<"choose from menue:\n\n"

        //DNA menu
        <<"1- Enter DNA Sequence ?\n"
        <<"2- Print DNA Sequence ?\n"
        <<"3- Compare if 2 DNA Sequences Equals ?\n"
        <<"4- Compare if 2 DNA Sequences Not Equals ?\n"
        <<"5- Sum 2 DNA Sequences ?\n"
        <<"6- Equal DNA Sequence to another One ?\n"
        <<"7- Convert DNA Sequence to RNA Sequence ?\n"

        //RNA menu
        <<"\n8- Enter RNA Sequence ?\n"
        <<"9- Print RNA Sequence ?\n"
        <<"10- Compare if 2 RNA Sequences Equals ?\n"
        <<"11- Compare if 2 RNA Sequences Not Equals ?\n"
        <<"12- Sum 2 RNA Sequences ?\n"
        <<"13- Equal RNA Sequence to another One ?\n"
        <<"14- Convert RNA Sequence to DNA Sequence ?\n"
        <<"15- Convert RNA Sequence to Protein Sequence ?\n"

        //Protein menu
        <<"\n16- Enter Protein Sequence ?\n"
        <<"17- Print Protein Sequence ?\n"
        <<"18- Compare if 2 Protein Sequences Equals ?\n"
        <<"19- Compare if 2 Protein Sequences Not Equals ?\n"
        <<"20- Sum 2 Protein Sequences ?\n"
        <<"21- Equal Protein Sequence to another One ?\n"
        <<"22- Convert Protein Sequence to DNA Sequence ?\n"

        //Some global Methods
        <<"23- Calculate the LCS Alignment for 2 Sequence ?\n"
        <<"24- Calculate the Global Alignment for 2 Sequence ?\n"
        <<"25- Set Codon : \n"

        //Files handling menu
        <<"26- Save a Sequence to File ?\n"
        <<"27- Load a Sequence from File ?\n"

        <<"\nNote: Enter -1 to End Process.\n";
    int result=0;

    while(result!=-1)
    {
        cout<<"\nYour Choice:  ";
        cin>>result;
        switch(result)
        {
        case 1:
        {
            DNA dna1;
            cin>>dna1;
            cout<<"your DNA saved.";
            cout<<dna1;
            continue;
        }
        case 2:
        {
            char* seq1="ACAGCGAT";
            DNA dna1(seq1,tail);
            dna1.setIndex(1,9);
            cout<<dna1;
            continue;
        }
        case 3:
        {
            DNA dna1,dna2;
            cout<<"Enter first DNA:\n";
            cin>>dna1;
            cout<<"Enter second DNA:\n";
            cin>>dna2;
            cout<<((dna1==dna2)?"Yes, they R Equal !":"No, they R not Equal !");
            continue;
        }
        case 4:
        {
            DNA dna1,dna2;
            cout<<"Enter first DNA:\n";
            cin>>dna1;
            cout<<"Enter second DNA:\n";
            cin>>dna2;
            cout<<((dna1!=dna2)?"Yes, they R not Equal !":"No, they R Equal !");
            continue;
        }
        case 5:
        {
            DNA dna1,dna2,dna3;
            cout<<"Enter first DNA:\n";
            cin>>dna1;
            cout<<"Enter second DNA:\n";
            cin>>dna2;
            dna3=dna1+dna2;
            cout<<dna3;
            continue;
        }
        case 6:
        {
            DNA dna1;
            cin>>dna1;
            cout<<"BY Copy Constructor the Second DNA:\n";
            DNA dna2=dna1;
            cout<<dna2;
            continue;
        }
        case 7:
        {
            DNA dna1;
            cin>>dna1;
            RNA rna1;
            rna1=dna1.ConvertToRNA();
            cout<<rna1;
            continue;
        }
        case 8:
        {
            RNA rna1;
            cin>>rna1;
            cout<<"\nyour RNA saved.\n";
            continue;
        }
        case 9:
        {
            RNA rna1;
            cin>>rna1;
            cout<<"\nyour RNA: \n";
            cout<<rna1;
            continue;
        }
        case 10:
        {
            RNA rna1,rna2;
            cout<<"Enter first RNA:\n";
            cin>>rna1;
            cout<<"Enter second RNA:\n";
            cin>>rna2;
            cout<<((rna1==rna2)?"Yes, they R Equal !":"No, they R not Equal !");
            continue;
        }
        case 11:
        {
            RNA rna1,rna2;
            cout<<"Enter first RNA:\n";
            cin>>rna1;
            cout<<"Enter second RNA:\n";
            cin>>rna2;
            cout<<((rna1!=rna2)?"Yes, they R not Equal !":"No, they R Equal !");
            continue;
        }
        case 12:
        {
            RNA rna1,rna2,rna3;
            cout<<"Enter first RNA:\n";
            cin>>rna1;
            cout<<"Enter second RNA:\n";
            cin>>rna2;
            rna3=rna1+rna2;
            cout<<rna3;
            continue;
        }
        case 13:
        {
            RNA rna1;
            cin>>rna1;
            cout<<"BY Copy Constructor the Second RNA:\n";
            RNA rna2=rna1;
            cout<<rna2;
            continue;
        }
        case 14:
        {
            RNA rna1;
            cin>>rna1;
            DNA dna1;
            dna1=rna1.ConvertToDNA();
            cout<<dna1;
            continue;
        }
        case 15:
        {
            RNA rna1;
            cin>>rna1;
            CodonsTable table;
            Protein ptr1;
            ptr1 = rna1.ConvertToProtein(table);
            cout<<ptr1;
            continue;
        }
        case 16:
        {
            Protein prt1;
            cin>>prt1;
            cout<<"\nyour Protein saved.\n";
            continue;
        }
        case 17:
        {
            Protein prt1;
            cin>>prt1;
            cout<<"\nyour Protein: \n";
            cout<<prt1;
            continue;
        }
        case 18:
        {
            Protein prt1,prt2;
            cout<<"Enter first Protein:\n";
            cin>>prt1;
            cout<<"Enter second Protein:\n";
            cin>>prt2;
            cout<<((prt1==prt2)?"Yes, they R Equal !":"No, they R not Equal !");
            continue;
        }
        case 19:
        {
            Protein prt1,prt2;
            cout<<"Enter first Protein:\n";
            cin>>prt1;
            cout<<"Enter second Protein:\n";
            cin>>prt2;
            cout<<((prt1!=prt2)?"Yes, they R not Equal !":"No, they R Equal !");
            continue;
        }
        case 20:
        {
            Protein prt1,prt2,prt3;
            cout<<"Enter first Protein:\n";
            cin>>prt1;
            cout<<"Enter second Protein:\n";
            cin>>prt2;
            prt3=prt1+prt2;
            cout<<prt3;
            continue;
        }
        case 21:
        {
            Protein prt1;
            cin>>prt1;
            cout<<"BY Copy Constructor the Second Protein:\n";
            Protein prt2=prt1;
            cout<<prt2;
            continue;
        }
        case 22:
        {
            Protein prt1;
            cin>>prt1;
            DNA dna1;
            cin>>dna1;
            DNA dna2;
            dna2=prt1.GetDNAStrandsEncodingMe(dna1);
            cout<<dna2;
            continue;
        }
        case 23:
        {
            DNA dna1("AGGTAC",tail);
            DNA dna2("GCTAATC",tail);
            DNA dna3;
            cout << "\nLCS of " << dna1.getSeq() << " and " << dna2.getSeq()<< " is " ;
            cout<<dna3.Align(dna1,dna2);                                         ///LCS alignment

            RNA rna1("AGCU",mRNA);
            RNA rna2("UCGAGAU",mRNA);
            RNA rna3;
            cout << "\nLCS of " << rna1.getSeq() << " and " << rna2.getSeq()<< " is " ;
            cout<<rna3.Align(rna1,rna2);

            Protein prt1("PLVWPY",TF);
            Protein prt2("KPLY",TF);
            Protein prt3;
            cout << "\nLCS of " << prt1.getSeq() << " and " << prt2.getSeq()<< " is " ;
            cout<<prt3.Align(prt1,prt2);

            continue;
        }
        case 24:
        {
            DNA dna1("ATTATCT",tail);
            DNA dna2("TTTCTA",tail);
            DNA dna3;
            cout << "\nLocal Alignment of " << dna1.getSeq() << " and " << dna2.getSeq()<< " is " ;
            cout<<dna3.LocalAlign(dna1,dna2);
            continue;
        }
        case 25:
        {
            CodonsTable table2;
            table2.LoadCodonsFromFile("RNA_codon_table.txt");
            int indx;
            cout<<"Enter your index : ";
            cin>>indx;
            char ch;
            cout<<"Enter your AminoAcid: ";
            cin>>ch;
            table2.setCodon(ch,indx);
            continue;
        }
        case 26:
        {
            string sequence ,FileName;
            cout<<"which sequence you want to save DNA - RNA - Protein ?!\nSequence: ";
            cin>>sequence;
            cout<<"Name of your file:  ";
            cin>>FileName;
            if(sequence=="DNA"){
                DNA dna1;
                char *f = new char[FileName.length() + 1];
                strcpy(f, FileName.c_str());
                dna1.SaveSequenceToFile(f);
            }
            else if(sequence=="RNA"){
                RNA rna1;
                char *f = new char[FileName.length() + 1];
                strcpy(f, FileName.c_str());
                rna1.SaveSequenceToFile(f);
            }
            else if(sequence=="Protein"){
                Protein prt1;
                char *f = new char[FileName.length() + 1];
                strcpy(f, FileName.c_str());
                prt1.SaveSequenceToFile(f);
            }
            else{
                cout<<"Invalid Type !";
            }

            continue;
        }
        case 27:
        {
            string sequence ,FileName;
            cout<<"which sequence you want to load DNA - RNA - Protein ?!\nSequence: ";
            cin>>sequence;
            cout<<"Name of your file:  ";
            cin>>FileName;
            if(sequence=="DNA"){
                DNA dna1;
                char *f = new char[FileName.length() + 1];
                strcpy(f, FileName.c_str());
                dna1.LoadSequenceFromFile(f);
            }
            else if(sequence=="RNA"){
                RNA rna1;
                char *f = new char[FileName.length() + 1];
                strcpy(f, FileName.c_str());
                rna1.LoadSequenceFromFile(f);
            }
            else if(sequence=="Protein"){
                Protein prt1;
                char *f = new char[FileName.length() + 1];
                strcpy(f, FileName.c_str());
                prt1.LoadSequenceFromFile(f);
            }
            else{
                cout<<"Invalid Type !";
            }


            continue;
        }

        default:
            return 0;
        }
    }
    return 0;
}
