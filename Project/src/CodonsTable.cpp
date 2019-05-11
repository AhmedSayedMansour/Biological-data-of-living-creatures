//#include "CodonsTable.h"
//#include "RNA.h"
//#include "Sequence .h"
//
//CodonsTable::CodonsTable()
//{
//    cout<<"HI";
//}
//
//CodonsTable::~CodonsTable()
//{
//    cout<<"HI";
//}
//void CodonsTable::LoadCodonsFromFile(string codonsFileName)
//{
//    file.open(codonsFileName,ios::in);
//    int i=0 ;
//    while ( file && !file.eof() )
//    {
//        file >> codons[i].value;
//        file >> codons[i++].AminoAcid;
//    }
//    file.close();
//}
//string CodonsTable::getAminoAcid(char* valuee ,int steps)
//{
//    vector<string> RNARNA;
//    string temp, out = "";
//
//    for ( int i = 0 ; i <= steps-3 ; i+=3 )
//    {
//        //temp = valuee.substr(i,3);
//        temp="";
//        temp += valuee[i];
//        temp += valuee[i+1];
//        temp += valuee[i+2];
//        RNARNA.push_back( temp );
//    }
//    for ( int i = 0 ; i < RNARNA.size() ; i++ )
//    {
//        for ( int j = 0; j < 61 ; j++ )
//        {
//            if ( RNARNA[i] == codons[j].value )
//            {
//                out += codons[j].AminoAcid;
//                break;
//            }
//            else if ( j == 60 )
//                out += ".";
//        }
//    }
//    return out;
//}
