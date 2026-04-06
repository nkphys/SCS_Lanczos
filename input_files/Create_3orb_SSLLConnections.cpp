#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include "tensor_type.h"

int main(){



    double tx=VALUE_TX;
    double ty=VALUE_TX;
    double tz=VALUE_TZ;


    double Delta=VALUE_DELTA; //Delta X n_{xy}

    double U=VALUE_U;
    double JbyU=VALUE_JBYU;
    double J;
    J=U*JbyU;



    bool H_SQ1_bool, H_SQ2_bool, H_SQ3_bool, H_SL_bool, H_Q1_bool, H_Q2_bool, H_Q3_bool, H_L_bool; 
    H_SQ1_bool=true;H_SQ2_bool=true;H_SQ3_bool=true;
    H_SL_bool=true;
    H_Q1_bool=true;H_Q2_bool=true;H_Q3_bool=true;
    H_L_bool=true;


    int OX, OY, OZ;
    OX=0;OY=1;OZ=2;

    //int L1=2; //X
    //int L2=2; //Y

    int N; 

    //sqrt(8)Xsqrt(8)
   /* 
    N=8;
    Mat_1_int Bonds_Xtype_site1 = {0, 1, 2, 3, 4, 5, 6, 7};
    Mat_1_int Bonds_Xtype_site2 = {6, 7, 1, 0, 2, 3, 5, 4};

    Mat_1_int Bonds_Ytype_site1 = {0, 1, 2, 3, 4, 5, 6, 7};
    Mat_1_int Bonds_Ytype_site2 = {2, 3, 5, 4, 6, 7, 1, 0};
    */
    

   //single bond +X
    /*
    N=2;
    Mat_1_int Bonds_Xtype_site1 = {0};
    Mat_1_int Bonds_Xtype_site2 = {1};

    Mat_1_int Bonds_Ytype_site1;
    Mat_1_int Bonds_Ytype_site2;
    Bonds_Ytype_site1.clear();
    Bonds_Ytype_site2.clear();
    */

  //2x2 system
    /*
    N=4;
    Mat_1_int Bonds_Xtype_site1 = {0, 2};
    Mat_1_int Bonds_Xtype_site2 = {1, 3};
    
    Mat_1_int Bonds_Ytype_site1 = {0, 1};
    Mat_1_int Bonds_Ytype_site2 = {2, 3};
   */

   //4x4 OBC(x) PBC(y)
    /* 
    N=16;
    Mat_1_int Bonds_Ytype_site1 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,  10,  11,  12,  13,  14,  15};
    Mat_1_int Bonds_Ytype_site2 = {1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11,   8,  13,  14,  15,  12};

    Mat_1_int Bonds_Xtype_site1 = {0, 1, 2, 3, 4, 5, 6,  7,  8,  9,  10, 11};
    Mat_1_int Bonds_Xtype_site2 = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
   */

    //6x4 OBC(x) PBC(y)
    N=24;
    Mat_1_int Bonds_Ytype_site1 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,  10,  11,  12,  13,  14,  15, 16, 17, 18 ,19, 20, 21, 22, 23};
    Mat_1_int Bonds_Ytype_site2 = {1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11,   8,  13,  14,  15,  12, 17, 18, 19, 16, 21, 22, 23, 20};

    Mat_1_int Bonds_Xtype_site1 = {0, 1, 2, 3, 4, 5, 6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    Mat_1_int Bonds_Xtype_site2 = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};



    string DeltaConnfile_str = "Delta.dat";
    ofstream DeltaConn_file_stream(DeltaConnfile_str.c_str());

    string H_SQConnfile_str = "H_SQ.dat";
    ofstream H_SQConn_file_stream(H_SQConnfile_str.c_str());

    string H_SLConnfile_str = "H_SL.dat";
    ofstream H_SLConn_file_stream(H_SLConnfile_str.c_str());

    string H_QConnfile_str = "H_Q.dat";
    ofstream H_QConn_file_stream(H_QConnfile_str.c_str());

    string H_LConnfile_str = "H_L.dat";
    ofstream H_LConn_file_stream(H_LConnfile_str.c_str());

    int site_neigh;
    int site;

    Mat_1_string Neigh_type;
    Mat_1_int x_offset, y_offset;
    Mat_1_doub hop_x_or_y;
    Mat_1_string opr_str, opr_str2;
    Neigh_type.push_back("plusX");x_offset.push_back(1);y_offset.push_back(0);hop_x_or_y.push_back(tx);
    Neigh_type.push_back("plusY");x_offset.push_back(0);y_offset.push_back(1);hop_x_or_y.push_back(ty);
    opr_str.push_back("y");opr_str.push_back("x");
    opr_str2.push_back("x");opr_str2.push_back("y");


    string str_temp;



    //Tetragomal spltting \Delta n_{xy}
    for(int i=0;i<N;i++){
        site=i + N;
//        cout<<"1 Sz2 "<<site<<" "<<1.0*Delta<<endl;
	DeltaConn_file_stream<<"1 Sz2 "<<site<<" "<<1.0*Delta<<endl;
    }

    //H_SQ + H_Q
    for(int neigh_type=0;neigh_type<Neigh_type.size();neigh_type++){
        double txy=hop_x_or_y[neigh_type];

        Mat_1_int Bonds_site1, Bonds_site2;
        if(Neigh_type[neigh_type]=="plusX"){
            Bonds_site1 = Bonds_Xtype_site1;
            Bonds_site2 = Bonds_Xtype_site2;
        }
        if(Neigh_type[neigh_type]=="plusY"){
            Bonds_site1 = Bonds_Ytype_site1;
            Bonds_site2 = Bonds_Ytype_site2;
        }

        for(int bond_no=0;bond_no<Bonds_site1.size();bond_no++){
            site=Bonds_site1[bond_no];
            site_neigh=Bonds_site2[bond_no];


            Mat_1_string OPR_type;
            Mat_1_doub Hop_type;

            str_temp = "S"+opr_str[neigh_type] + "2 " + "S"+opr_str[neigh_type] + "2";
            OPR_type.push_back(str_temp);Hop_type.push_back(txy);
            OPR_type.push_back("Sz2 Sz2");Hop_type.push_back(tz);

            Mat_1_string SS_type;
            Mat_1_doub fac;
            SS_type.push_back("Sz Sz");fac.push_back(1.0);
            SS_type.push_back("Sp Sm");fac.push_back(0.5);
            SS_type.push_back("Sm Sp");fac.push_back(0.5);

            if(H_SQ1_bool){
                for(int typeL=0;typeL<2;typeL++){
                    for(int type=0;type<3;type++){
                        H_SQConn_file_stream<<"4 "<<SS_type[type]<<" "<<OPR_type[typeL]<<" ";
                        H_SQConn_file_stream<<site<<" "<<site_neigh<<" "<<site+N<<" "<<site_neigh+N<<" ";
                        H_SQConn_file_stream<<fac[type]*Hop_type[typeL]*Hop_type[typeL]*(1.0/U)*( ((U+J)/(U+2.0*J)) + (2*J/(U-3.0*J)) )<<endl;
                    }
                }
            }


            OPR_type.clear();
            Hop_type.clear();
            str_temp = "S"+opr_str[neigh_type] + "2";
            OPR_type.push_back(str_temp);Hop_type.push_back(txy);
            OPR_type.push_back("Sz2");Hop_type.push_back(tz);

            if(H_SQ2_bool){
                for(int typeL=0;typeL<2;typeL++){
                    for(int type=0;type<3;type++){
                        H_SQConn_file_stream<<"3 "<<SS_type[type]<<" "<<OPR_type[typeL]<<" ";
                        H_SQConn_file_stream<<site<<" "<<site_neigh<<" "<<site+N<<" ";
                        H_SQConn_file_stream<<fac[type]*Hop_type[typeL]*Hop_type[typeL]*(-1.0/U)*((J/(U-3.0*J)) )<<endl;

                        H_SQConn_file_stream<<"3 "<<SS_type[type]<<" "<<OPR_type[typeL]<<" ";
                        H_SQConn_file_stream<<site<<" "<<site_neigh<<" "<<site_neigh+N<<" ";
                        H_SQConn_file_stream<<fac[type]*Hop_type[typeL]*Hop_type[typeL]*(-1.0/U)*((J/(U-3.0*J)) )<<endl;
                    }
                }
            }


            string oprAdd;
            double facAdd;

            if(H_SQ3_bool){
                if(Neigh_type[neigh_type]=="plusX"){
                    oprAdd="i"; facAdd=-1.0;
                }
                else{
                    oprAdd=""; facAdd=1.0;
                }

                str_temp = oprAdd + "Q"+opr_str[neigh_type] + "z"+ " " + oprAdd + "Q"+opr_str[neigh_type] + "z";
                for(int type=0;type<3;type++){
                    H_SQConn_file_stream<<"4 "<<SS_type[type]<<" "<<str_temp<<" ";
                    H_SQConn_file_stream<<site<<" "<<site_neigh<<" "<<site+N<<" "<<site_neigh + N<<" ";
                    H_SQConn_file_stream<<facAdd*fac[type]*txy*tz*(1.0/U)*( 0.5*(J/(U+2.0*J)) + 0.5*((U-J)/(U-3.0*J)) )<<endl;
                }
            }



            //H_Q
            OPR_type.clear();
            Hop_type.clear();
            str_temp = "S"+opr_str[neigh_type] + "2 " + "S"+opr_str[neigh_type] + "2";
            OPR_type.push_back(str_temp);Hop_type.push_back(txy);
            OPR_type.push_back("Sz2 Sz2");Hop_type.push_back(tz);

            if(H_Q1_bool){
                for(int typeL=0;typeL<2;typeL++){
                    H_QConn_file_stream<<"2 "<<OPR_type[typeL]<<" ";
                    H_QConn_file_stream<<site+N<<" "<<site_neigh+N<<" ";
                    H_QConn_file_stream<<Hop_type[typeL]*Hop_type[typeL]*(1.0/U)*( (-1.0*(U+J)/(U+2.0*J)) + (2*(U-J)/(U-3.0*J)) )<<endl;
                }
            }

            OPR_type.clear();
            Hop_type.clear();
            str_temp = "S"+opr_str[neigh_type] + "2";
            OPR_type.push_back(str_temp);Hop_type.push_back(txy);
            OPR_type.push_back("Sz2");Hop_type.push_back(tz);

            if(H_Q2_bool){
                for(int typeL=0;typeL<2;typeL++){
                    H_QConn_file_stream<<"1 "<<OPR_type[typeL]<<" ";
                    H_QConn_file_stream<<site+N<<" ";
                    H_QConn_file_stream<<Hop_type[typeL]*Hop_type[typeL]*(-1.0/U)*(((U-J)/(U-3.0*J)) )<<endl;

                    H_QConn_file_stream<<"1 "<<OPR_type[typeL]<<" ";
                    H_QConn_file_stream<<site_neigh+N<<" ";
                    H_QConn_file_stream<<Hop_type[typeL]*Hop_type[typeL]*(-1.0/U)*(((U-J)/(U-3.0*J)) )<<endl;
                }
            }


            if(H_Q3_bool){
                if(Neigh_type[neigh_type]=="plusX"){
                    oprAdd="i"; facAdd=-1.0;
                }
                else{
                    oprAdd=""; facAdd=1.0;
                }

                str_temp = oprAdd + "Q"+opr_str[neigh_type] + "z"+ " " + oprAdd + "Q"+opr_str[neigh_type] + "z";
                H_QConn_file_stream<<"2 "<<str_temp<<" ";
                H_QConn_file_stream<<site+N<<" "<<site_neigh + N<<" ";
                H_QConn_file_stream<<facAdd*txy*tz*(1.0/U)*( 0.5*(-1.0*J/(U+2.0*J)) + 0.5*((U+J)/(U-3.0*J)) )<<endl;
            }

        }
    }//neigh





    //H_SL + H_L

    for(int neigh_type=0;neigh_type<Neigh_type.size();neigh_type++){
        double txy=hop_x_or_y[neigh_type];


        Mat_1_int Bonds_site1, Bonds_site2;
        if(Neigh_type[neigh_type]=="plusX"){
            Bonds_site1 = Bonds_Xtype_site1;
            Bonds_site2 = Bonds_Xtype_site2;
        }
        if(Neigh_type[neigh_type]=="plusY"){
            Bonds_site1 = Bonds_Ytype_site1;
            Bonds_site2 = Bonds_Ytype_site2;
        }

        for(int bond_no=0;bond_no<Bonds_site1.size();bond_no++){
            site=Bonds_site1[bond_no];
            site_neigh=Bonds_site2[bond_no];


            Mat_1_string OPR_type;
            Mat_1_doub Hop_type;


            string oprAdd;
            double facAdd;
            if(Neigh_type[neigh_type]=="plusY"){
                oprAdd="i"; facAdd=-1.0;
            }
            else{
                oprAdd=""; facAdd=1.0;
            }
            str_temp = oprAdd + "S"+opr_str2[neigh_type] + " " + oprAdd + "S"+opr_str2[neigh_type];
            OPR_type.push_back(str_temp);Hop_type.push_back(txy);

            Mat_1_string SS_type;
            Mat_1_doub fac;
            SS_type.push_back("Sz Sz");fac.push_back(1.0);
            SS_type.push_back("Sp Sm");fac.push_back(0.5);
            SS_type.push_back("Sm Sp");fac.push_back(0.5);




            if(H_SL_bool){
                for(int typeL=0;typeL<1;typeL++){
                    for(int type=0;type<3;type++){
                        H_SLConn_file_stream<<"4 "<<SS_type[type]<<" "<<OPR_type[typeL]<<" ";
                        H_SLConn_file_stream<<site<<" "<<site_neigh<<" "<<site+N<<" "<<site_neigh+N<<" ";
                        H_SLConn_file_stream<<facAdd*fac[type]*Hop_type[typeL]*tz*(0.5/U)*( ((-1.0*J)/(U+2.0*J)) + ((U-J)/(U-3.0*J)) )<<endl;
                    }
                }
            }

            if(H_L_bool){
                for(int typeL=0;typeL<1;typeL++){
                    H_LConn_file_stream<<"2 "<<OPR_type[typeL]<<" ";
                    H_LConn_file_stream<<site+N<<" "<<site_neigh+N<<" ";
                    H_LConn_file_stream<<facAdd*Hop_type[typeL]*tz*(0.5/U)*( ((1.0*J)/(U+2.0*J)) + ((U+J)/(U-3.0*J)) )<<endl;
                }
            }

        }
    }//neigh





    return 0;
}
