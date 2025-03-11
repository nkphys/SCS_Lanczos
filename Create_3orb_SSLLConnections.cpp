#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include "tensor_type.h"

int main(){


    double tx=1.0;
    double ty=1.0;
    double tz=1.4;


    double Delta=-2.0; //Delta X n_{xy}

    double U=20.0;
    double JbyU=0.30;
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
    N=8;
    Mat_1_int Bonds_Xtype_site1 = {0, 1, 2, 3, 4, 5, 6, 7};
    Mat_1_int Bonds_Xtype_site2 = {6, 7, 1, 0, 2, 3, 5, 4};

    Mat_1_int Bonds_Ytype_site1 = {0, 1, 2, 3, 4, 5, 6, 7};
    Mat_1_int Bonds_Ytype_site2 = {2, 3, 5, 4, 6, 7, 1, 0};

    


    string Connfile_str = "ConnFile.dat";
    ofstream Conn_file_stream(Connfile_str.c_str());


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
        cout<<"1 Sz2 "<<site<<" "<<1.0*Delta<<endl;
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
                        cout<<"4 "<<SS_type[type]<<" "<<OPR_type[typeL]<<" ";
                        cout<<site<<" "<<site_neigh<<" "<<site+N<<" "<<site_neigh+N<<" ";
                        cout<<fac[type]*Hop_type[typeL]*Hop_type[typeL]*(1.0/U)*( ((U+J)/(U+2.0*J)) + (2*J/(U-3.0*J)) )<<endl;
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
                        cout<<"3 "<<SS_type[type]<<" "<<OPR_type[typeL]<<" ";
                        cout<<site<<" "<<site_neigh<<" "<<site+N<<" ";
                        cout<<fac[type]*Hop_type[typeL]*Hop_type[typeL]*(-1.0/U)*((J/(U-3.0*J)) )<<endl;

                        cout<<"3 "<<SS_type[type]<<" "<<OPR_type[typeL]<<" ";
                        cout<<site<<" "<<site_neigh<<" "<<site_neigh+N<<" ";
                        cout<<fac[type]*Hop_type[typeL]*Hop_type[typeL]*(-1.0/U)*((J/(U-3.0*J)) )<<endl;
                    }
                }
            }


            string oprAdd;
            double facAdd;

            if(H_SQ3_bool){
                if(Neigh_type[neigh_type]=="plusX"){
                    oprAdd="I"; facAdd=-1.0;
                }
                else{
                    oprAdd=""; facAdd=1.0;
                }

                str_temp = "Q"+opr_str[neigh_type] + "z"+oprAdd+" " + "Q"+opr_str[neigh_type] + "z" + oprAdd;
                for(int type=0;type<3;type++){
                    cout<<"4 "<<SS_type[type]<<" "<<str_temp<<" ";
                    cout<<site<<" "<<site_neigh<<" "<<site+N<<" "<<site_neigh + N<<" ";
                    cout<<facAdd*fac[type]*txy*tz*(1.0/U)*( 0.5*(J/(U+2.0*J)) + 0.5*((U-J)/(U-3.0*J)) )<<endl;
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
                    cout<<"2 "<<OPR_type[typeL]<<" ";
                    cout<<site+N<<" "<<site_neigh+N<<" ";
                    cout<<Hop_type[typeL]*Hop_type[typeL]*(1.0/U)*( (-1.0*(U+J)/(U+2.0*J)) + (2*(U-J)/(U-3.0*J)) )<<endl;
                }
            }

            OPR_type.clear();
            Hop_type.clear();
            str_temp = "S"+opr_str[neigh_type] + "2";
            OPR_type.push_back(str_temp);Hop_type.push_back(txy);
            OPR_type.push_back("Sz2");Hop_type.push_back(tz);

            if(H_Q2_bool){
                for(int typeL=0;typeL<2;typeL++){
                    cout<<"1 "<<OPR_type[typeL]<<" ";
                    cout<<site+N<<" ";
                    cout<<Hop_type[typeL]*Hop_type[typeL]*(-1.0/U)*(((U-J)/(U-3.0*J)) )<<endl;

                    cout<<"1 "<<OPR_type[typeL]<<" ";
                    cout<<site_neigh+N<<" ";
                    cout<<Hop_type[typeL]*Hop_type[typeL]*(-1.0/U)*(((U-J)/(U-3.0*J)) )<<endl;
                }
            }


            if(H_Q3_bool){
                if(Neigh_type[neigh_type]=="plusX"){
                    oprAdd="I"; facAdd=-1.0;
                }
                else{
                    oprAdd=""; facAdd=1.0;
                }

                str_temp = "Q"+opr_str[neigh_type] + "z"+oprAdd+" " + "Q"+opr_str[neigh_type] + "z" + oprAdd;
                cout<<"2 "<<str_temp<<" ";
                cout<<site+N<<" "<<site_neigh + N<<" ";
                cout<<facAdd*txy*tz*(1.0/U)*( 0.5*(-1.0*J/(U+2.0*J)) + 0.5*((U+J)/(U-3.0*J)) )<<endl;
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
                oprAdd="I"; facAdd=-1.0;
            }
            else{
                oprAdd=""; facAdd=1.0;
            }
            str_temp = "S"+opr_str2[neigh_type] + oprAdd + " " + "S"+opr_str2[neigh_type] + oprAdd;
            OPR_type.push_back(str_temp);Hop_type.push_back(txy);

            Mat_1_string SS_type;
            Mat_1_doub fac;
            SS_type.push_back("Sz Sz");fac.push_back(1.0);
            SS_type.push_back("Sp Sm");fac.push_back(0.5);
            SS_type.push_back("Sm Sp");fac.push_back(0.5);




            if(H_SL_bool){
                for(int typeL=0;typeL<1;typeL++){
                    for(int type=0;type<3;type++){
                        cout<<"4 "<<SS_type[type]<<" "<<OPR_type[typeL]<<" ";
                        cout<<site<<" "<<site_neigh<<" "<<site+N<<" "<<site_neigh+N<<" ";
                        cout<<facAdd*fac[type]*Hop_type[typeL]*tz*(0.5/U)*( ((-1.0*J)/(U+2.0*J)) + ((U-J)/(U-3.0*J)) )<<endl;
                    }
                }
            }

            if(H_L_bool){
                for(int typeL=0;typeL<1;typeL++){
                    cout<<"2 "<<OPR_type[typeL]<<" ";
                    cout<<site+N<<" "<<site_neigh+N<<" ";
                    cout<<facAdd*Hop_type[typeL]*tz*(0.5/U)*( ((1.0*J)/(U+2.0*J)) + ((U+J)/(U-3.0*J)) )<<endl;
                }
            }

        }
    }//neigh





    return 0;
}
