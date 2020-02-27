#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include "iconprint.h"
#include "models/Model_3_orb_Hubbard_chain.h"
#include "models/Model_2_orb_Hubbard_chain.h"
#include "models/Model_2_orb_Hubbard_chain_KSector.h"
#include "models/Model_3_orb_Hubbard_chain_two_SzSectors.h"
#include "basis/Basis_3_orb_Hubbard_chain.h"
#include "basis/Basis_2_orb_Hubbard_chain.h"
#include "basis/Basis_2_orb_Hubbard_chain_KSector.h"
#include "basis/Basis_3_orb_Hubbard_chain_two_SzSectors.h"
#include "models/Model_1_orb_Hubbard_chain.h"
#include "models/Model_1_orb_tJ.h"
#include "basis/Basis_1_orb_Hubbard_chain.h"
#include "basis/Basis_1_orb_Hubbard_GC.h"
#include "basis/Basis_1_orb_tJ.h"
#include "basis/Basis_Spins.h"
#include "models/Model_Spins.h"
#include "models/Model_3_orb_Hubbard_chain_GC.h"
#include "models/Model_1_orb_Hubbard_GC.h"
#include "basis/Basis_3_orb_Hubbard_chain_GC.h"
#include "basis/Basis_3_orb_Hubbard_chain_GC_restricted.h"
#include "Lanczos_engine.h"
#include "reading_input.h"

//Remember "cpp" files for templated class over basis need to be included in this code
#include "models/Model_3_orb_Hubbard_chain_GC.cpp"
#include "models/Model_1_orb_Hubbard_GC.cpp"






int main(int argc, char** argv){


    //iconprint();

    string model_name;
    string inp_filename = argv[1];
    bool Do_Dynamics;
    bool Perform_RIXS;
    bool Restricted_Basis;

    reading_input(inp_filename, model_name, Do_Dynamics, Perform_RIXS, Restricted_Basis);
    cout<<"Do_Dynamics ="<<Do_Dynamics<<endl;


    bool DO_FULL_DIAGONALIZATION=false;


    if(model_name=="SpinOnly"){
        MODEL_Spins _MODEL;
        BASIS_Spins _BASIS;

        _MODEL.Read_parameters(_BASIS, inp_filename);
        _BASIS.Construct_basis();


        _MODEL.Add_diagonal_terms(_BASIS);
        _MODEL.Add_non_diagonal_terms(_BASIS);
        _MODEL.Add_connections(_BASIS);


        LANCZOS _LANCZOS;
        _LANCZOS.Dynamics_performed=false;
        _LANCZOS.Read_Lanczos_parameters(inp_filename);


        _LANCZOS.Perform_LANCZOS(_MODEL.Hamil);
        _LANCZOS.Write_full_spectrum();
        Print_vector_in_file(_LANCZOS.Eig_vec,"seed_GS.txt");

        _MODEL.Initialize_one_point_to_calculate_from_file(_BASIS);
        _LANCZOS.Measure_one_point_observables(_MODEL.one_point_obs, _MODEL.One_point_oprts, _BASIS.Length, 0);



        cout<<"Dynamics startedXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;

        if(Do_Dynamics){

            int site_=0;
            _MODEL.Read_parameters_for_dynamics(inp_filename);

            LANCZOS _LANCZOS_Dynamics;
            _LANCZOS_Dynamics.Dynamics_performed=true;
            _LANCZOS_Dynamics.Read_Lanczos_parameters(inp_filename);
            _LANCZOS_Dynamics.Eig_vec=_LANCZOS.Eig_vec;
            _LANCZOS_Dynamics.GS_energy=_LANCZOS.GS_energy;

            _MODEL.Initialize_Opr_for_Dynamics(_BASIS, site_);
            _LANCZOS_Dynamics.Get_Dynamics_seed(_MODEL.Dyn_opr);

            cout<<"size of seed = "<<_LANCZOS_Dynamics.Dynamics_seed.size()<<endl;
            _LANCZOS_Dynamics.omega_sign=1.0;
            _LANCZOS_Dynamics.file_dynamics_out = _LANCZOS_Dynamics.file_dynamics_out +".txt";
            _LANCZOS_Dynamics.Perform_LANCZOS(_MODEL.Hamil);

            //-----------------------

        }


    }

    if (model_name=="2_orb_Hubbard_chain_KSector") {

        MODEL_2_orb_Hubb_chain_KSector _MODEL;
        BASIS_2_orb_Hubb_chain_KSector _BASIS;

        _MODEL.Read_parameters(_BASIS, inp_filename);
        _BASIS.Construct_basis();


        _MODEL.Add_diagonal_terms(_BASIS);
        _MODEL.Add_non_diagonal_terms(_BASIS);
        _MODEL.Add_connections(_BASIS);
        cout<<"Size of Hilbert space = "<<_MODEL.Hamil.nrows<<endl;
        //cout<<scientific<<setprecision(2);
        //Print_Matrix_COO(_MODEL.Hamil);
        cout<<scientific<<setprecision(8);


        LANCZOS _LANCZOS;
        _LANCZOS.Dynamics_performed=false;
        _LANCZOS.Read_Lanczos_parameters(inp_filename);


        _LANCZOS.Perform_LANCZOS(_MODEL.Hamil);
        _LANCZOS.Write_full_spectrum();
        Print_vector_in_file(_LANCZOS.Eig_vec,"seed_GS.txt");


        // Print_file_in_vector(_LANCZOS.Eig_vec,_LANCZOS.seed_file_in,_BASIS.D_Period.size());
        // cout<<"_LANCZOS.Eig_vec done"<<endl;


        assert(true==false);

        double_type overlap_temp;

        //ONLY FOR HALF-FILLING
        if(_BASIS.Nup==_BASIS.Ndn &&
                (_BASIS.Nup==_BASIS.Length)
                ){
            _MODEL.Perform_RVB_Analysis_at_half_filling(_BASIS, _LANCZOS.Eig_vec);
            _MODEL.Perform_Product_state_Analysis_at_half_filling(_BASIS, _LANCZOS.Eig_vec);
            overlap_temp = dot_product(_MODEL.STATE_RVB, _MODEL.STATE_TPS);
            cout<<"<STATE_RVB|STATE_TPS> = "<<overlap_temp<<endl;
        }

        //ONLY FOR 2holes doped on Half-filling, and Sz=0 sector
        if(_BASIS.Nup==_BASIS.Ndn
                &&
                (_BASIS.Nup==(_BASIS.Length -1))
                ){

            _MODEL.Perform_RVB_Analysis_at_2_hole_doped(_BASIS, _LANCZOS.Eig_vec);
            _MODEL.Perform_Product_state_Analysis_at_2_hole_doped(_BASIS, _LANCZOS.Eig_vec);

            overlap_temp = dot_product(_MODEL.STATE_RVB, _MODEL.STATE_TPS);
            cout<<"<STATE_RVB|STATE_TPS> = "<<overlap_temp<<endl;
        }


        //Checking Orbital Symmetry of Ground State
        Mat_1_doub vec1_ ;
        double_type norm1_temp;
        double norm1_, overlap_temp2;

        vec1_ = _MODEL.Act_Orbital_Exchange(_BASIS, _LANCZOS.Eig_vec);
        norm1_temp = dot_product(vec1_,vec1_);
        norm1_ = abs(norm1_temp);
        for(int i=0;i<vec1_.size();i++){
            vec1_[i] = vec1_[i]*sqrt(1.0/norm1_);
        }
        overlap_temp = dot_product(_LANCZOS.Eig_vec,vec1_);
        overlap_temp2 = abs(overlap_temp)*abs(overlap_temp);
        cout <<"|GS> is an eigenstate of Orbital Symmetry Operator with eigenvalue = "<< overlap_temp<<endl;


        vec1_ = _MODEL.Act_Reflection_about_Central_site(_BASIS, _LANCZOS.Eig_vec);
        norm1_temp = dot_product(vec1_,vec1_);
        norm1_ = abs(norm1_temp);
        for(int i=0;i<vec1_.size();i++){
            vec1_[i] = vec1_[i]*sqrt(1.0/norm1_);
        }
        overlap_temp = dot_product(_LANCZOS.Eig_vec,vec1_);
        overlap_temp2 = abs(overlap_temp)*abs(overlap_temp);
        cout <<"|GS> is an eigenstate of Parity Operator [Ref. about plane cutting through center] with eigenvalue = "<< overlap_temp<<endl;




        //Spin Spin corr
        if(true && _BASIS.Momentum_n==0){

            Matrix_COO OPR_;
            double_type sum_;

            Mat_1_string opr_type_;
            opr_type_.push_back("SzSz");
            opr_type_.push_back("SpSm");
            opr_type_.push_back("SmSp");

            Mat_2_doub Corr_orb1_orb2_;
            Corr_orb1_orb2_.resize(_BASIS.Length);
            for(int site1=0;site1<_BASIS.Length;site1++){
                Corr_orb1_orb2_[site1].resize(_BASIS.Length);
            }

            for(int type=0;type<2;type++){
                sum_=zero;

                for(int gamma1=0;gamma1<2;gamma1++){
                    for(int gamma2=gamma1;gamma2<2;gamma2++){
                        cout<<opr_type_[type]<<gamma1<<gamma2<<"[site1][site2]: "<<endl;
                        cout<<"-----------------------------------------------------"<<endl;
                        for(int site1=0;site1<_BASIS.Length;site1++){
                            for(int site2=0;site2<_BASIS.Length;site2++){

                                if(site2>=site1){
                                    OPR_.columns.clear();
                                    OPR_.rows.clear();
                                    OPR_.value.clear();
                                    _MODEL.Initialize_two_point_operator_sites_orbital_specific(opr_type_[type] , OPR_, site1,gamma1, site2,gamma2, _BASIS);
                                    Corr_orb1_orb2_[site1][site2]=_LANCZOS.Measure_observable(OPR_, 0);

                                    cout<< Corr_orb1_orb2_[site1][site2]<<" ";
                                    OPR_.columns.clear();
                                    OPR_.rows.clear();
                                    OPR_.value.clear();

                                    if(site1==site2 && gamma1==gamma2){
                                        sum_ += Corr_orb1_orb2_[site1][site2];
                                    }
                                    else if(site1==site2 && gamma1!=gamma2){
                                        sum_ +=2.0*Corr_orb1_orb2_[site1][site2];
                                    }
                                    else if(site1 !=site2 && gamma1==gamma2){
                                        sum_ +=2.0*Corr_orb1_orb2_[site1][site2];
                                    }
                                    else if(site1 !=site2 && gamma1!=gamma2){
                                        sum_ +=4.0*Corr_orb1_orb2_[site1][site2];
                                    }

                                }
                                else{
                                    cout<< zero<<" ";
                                }

                            }
                            cout<<endl;

                        }
                        cout<<"------------------------------------------------------"<<endl;
                    }
                }
                if(type==0){
                    cout<<"Total S_{z}^2 = "<<sum_<<endl;
                }
                if(type==1){
                    cout<<"Total (SpSm) = "<<sum_<<endl;
                }
                cout<<endl;
                cout<<endl;
            }
        }


    }


    //************************************//

    //=======================================================================================================================================================================================================//
    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //=======================================================================================================================================================================================================//


    if (model_name=="3_orb_Hubbard_chain") {
        MODEL_3_orb_Hubb_chain _MODEL;
        BASIS_3_orb_Hubb_chain _BASIS;

        //---------------------------------------

        _MODEL.Read_parameters(_BASIS, inp_filename);
        _BASIS.Construct_basis();
        _MODEL.Add_diagonal_terms(_BASIS);
        _MODEL.Add_non_diagonal_terms(_BASIS);
        _MODEL.Add_connections(_BASIS);



        cout<<"Size of Hilbert space = "<<_MODEL.Hamil.nrows<<endl;
        cout<<scientific<<setprecision(8);
        //Print_Matrix_COO(_MODEL.Hamil);

        double EG;
        Mat_1_doub vecG;
        if(DO_FULL_DIAGONALIZATION==true){
            Diagonalize(_MODEL.Hamil, EG, vecG);
            cout<<"GS energy from ED(without Lanczos) = "<<EG<<endl;
        }

        LANCZOS _LANCZOS;
        _LANCZOS.Dynamics_performed=false;
        _LANCZOS.Read_Lanczos_parameters(inp_filename);
        _LANCZOS.Perform_LANCZOS(_MODEL.Hamil);
        _LANCZOS.Write_full_spectrum();
        Print_vector_in_file(_LANCZOS.Eig_vec,"seed_GS.txt");



        if(Perform_RIXS){
            BASIS_3_orb_Hubb_chain_two_SzSectors _BASIS_Np1;
            MODEL_3_orb_Hubb_chain_two_SzSectors _MODEL_Np1;

            _MODEL_Np1.Read_parameters(_BASIS_Np1, inp_filename);
            _BASIS_Np1.Construct_basis();

            _MODEL_Np1.Add_diagonal_terms(_BASIS_Np1);
            _MODEL_Np1.Add_non_diagonal_terms(_BASIS_Np1);
            _MODEL_Np1.Add_connections(_BASIS_Np1);

            LANCZOS _LANCZOS_Np1;
            _LANCZOS_Np1.Read_Lanczos_parameters(inp_filename);
            _LANCZOS_Np1.Perform_LANCZOS(_MODEL_Np1.Hamil);


        }
        else{
            _MODEL.Initialize_one_point_to_calculate(_BASIS);
            _LANCZOS.Measure_one_point_observables(_MODEL.one_point_obs, _MODEL.One_point_oprts, _BASIS.Length, 0);


            _MODEL.Initialize_two_point_to_calculate(_BASIS);
            for(int i=0;i<_LANCZOS.states_to_look.size();i++){
                _LANCZOS.Measure_two_point_observables(_MODEL.two_point_obs, _MODEL.Two_point_oprts, _BASIS.Length, i, _MODEL.PBC);
            }

            _MODEL.Initialize_macro_oprs_to_calculate(_BASIS);
            for(int i=0;i<_LANCZOS.states_to_look.size();i++){
                _LANCZOS.Measure_macro_observables(_MODEL.macro_obs, _MODEL.Macro_oprts ,i);
            }

            _MODEL.Calculate_Local_Obs_for_States_to_Look(_LANCZOS,_BASIS);

        }



        if(Do_Dynamics){

            _MODEL.Read_parameters_for_dynamics(inp_filename);

            _MODEL.Initialize_Opr_for_Dynamics(_BASIS);

            LANCZOS _LANCZOS_Dynamics;
            _LANCZOS_Dynamics.Dynamics_performed=true;
            _LANCZOS_Dynamics.Read_Lanczos_parameters(inp_filename);
            _LANCZOS_Dynamics.Eig_vec=_LANCZOS.Eig_vec;
            _LANCZOS_Dynamics.GS_energy=_LANCZOS.GS_energy;
            _LANCZOS_Dynamics.Get_Dynamics_seed(_MODEL.Dyn_opr);

            _LANCZOS_Dynamics.Perform_LANCZOS(_MODEL.Hamil);

        }



    }

    //=======================================================================================================================================================================================================//
    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //=======================================================================================================================================================================================================//





    //=======================================================================================================================================================================================================//
    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //=======================================================================================================================================================================================================//


#ifndef USE_COMPLEX
    if (model_name=="2_orb_Hubbard_chain") {
        bool PERFORM_VARIATIONAL_STATE_ANALYSIS=false;
        bool Read_Half_filling_Variational_states=false;
        bool PERFORM_Nm2_Variational_Calculation=false;
        bool Cheaper_SpinSpincorr=true;



        MODEL_2_orb_Hubb_chain _MODEL;
        BASIS_2_orb_Hubb_chain _BASIS;

        //---------------------------------------

        _MODEL.Read_parameters(_BASIS, inp_filename);
        _BASIS.Construct_basis();

        // _MODEL.Read_parameters_for_variational_state(_BASIS, inp_filename);

        _MODEL.Add_diagonal_terms(_BASIS);
        _MODEL.Add_non_diagonal_terms(_BASIS);
        _MODEL.Add_connections(_BASIS);



        if(_MODEL.Hamil.nrows<400){
            DO_FULL_DIAGONALIZATION=true;
        }
        if(DO_FULL_DIAGONALIZATION==true){
            double EG;
            Mat_1_real Evals_temp;
            Mat_1_doub vecG;
            Diagonalize(_MODEL.Hamil, Evals_temp, vecG);
            cout<<"GS energy from ED(without Lanczos) = "<<Evals_temp[0]<<endl;
            cout<<"All eigenvalues using EG------------------------------"<<endl;
            cout<<"-------------------------------------------------------"<<endl;
            for(int i=0;i<Evals_temp.size();i++){
                cout<<i<<"  "<<Evals_temp[i]<<endl;
            }
            cout<<"-------------------------------------------------------"<<endl;
            cout<<"-------------------------------------------------------"<<endl;

        }





        cout<<"Size of Hilbert space = "<<_MODEL.Hamil.nrows<<endl;
        cout<<scientific<<setprecision(8);
        //Print_Matrix_COO(_MODEL.Hamil);

        double EG;
        Mat_1_doub vecG;

        LANCZOS _LANCZOS;
        _LANCZOS.Dynamics_performed=false;
        _LANCZOS.Read_Lanczos_parameters(inp_filename);
        _LANCZOS.Perform_LANCZOS(_MODEL.Hamil);
        _LANCZOS.Write_full_spectrum();
        Print_vector_in_file(_LANCZOS.Eig_vec,"seed_GS.txt");

        //assert(false);
        cout<<"Orbital symmetry of Ground state for without doping:"<<endl;
        _MODEL.Check_orbital_symmetry(_BASIS, _LANCZOS.Eig_vec);

        _MODEL.Initialize_one_point_to_calculate(_BASIS);
        _LANCZOS.Measure_one_point_observables(_MODEL.one_point_obs, _MODEL.One_point_oprts, _BASIS.Length, 0);


        if(Cheaper_SpinSpincorr==true){

            Matrix_COO OPR_;
            for(int state_=0;state_<_LANCZOS.states_to_look.size();state_++){
                cout<<"Spin-Spin correlations for state = "<<state_<<endl;

                double sum_;

                Mat_1_string opr_type_;
                opr_type_.push_back("SzSz");
                opr_type_.push_back("SpSm");
                opr_type_.push_back("SmSp");

                Mat_2_doub Corr_;
                Corr_.resize(_BASIS.Length);
                for(int site1=0;site1<_BASIS.Length;site1++){
                    Corr_[site1].resize(_BASIS.Length);
                }

                for(int type=0;type<3;type++){
                    sum_=0.0;
                    cout<<opr_type_[type]<<": "<<endl;

                    for(int site1=0;site1<_BASIS.Length;site1++){
                        for(int site2=site1;site2<_BASIS.Length;site2++){
                            OPR_.columns.clear();
                            OPR_.rows.clear();
                            OPR_.value.clear();
                            _MODEL.Initialize_two_point_operator_sites_specific(opr_type_[type] , OPR_, site1, site2, _BASIS);
                            Corr_[site1][site2]=_LANCZOS.Measure_observable(OPR_, state_);
                            if(site1 != site2){
                                Corr_[site2][site1]=Corr_[site1][site2];
                            }
                            OPR_.columns.clear();
                            OPR_.rows.clear();
                            OPR_.value.clear();
                        }
                    }
                    for(int site1=0;site1<_BASIS.Length;site1++){
                        for(int site2=0;site2<_BASIS.Length;site2++){
                            cout<< Corr_[site1][site2]<<" ";
                            sum_ +=Corr_[site1][site2];
                        }
                        cout<<endl;
                    }
                    cout<<"sum = "<<sum_<<endl;
                }

                vector< int >().swap( OPR_.columns );
                vector< int >().swap( OPR_.rows );
                vector< double_type >().swap( OPR_.value );
            }
        }
        else{
            _MODEL.Initialize_two_point_to_calculate(_BASIS);
            for(int i=0;i<_LANCZOS.states_to_look.size();i++){
                _LANCZOS.Measure_two_point_observables(_MODEL.two_point_obs, _MODEL.Two_point_oprts, _BASIS.Length, i, _MODEL.PBC);
            }
        }

        _MODEL.Initialize_macro_oprs_to_calculate(_BASIS);
        for(int i=0;i<_LANCZOS.states_to_look.size();i++){
            _LANCZOS.Measure_macro_observables(_MODEL.macro_obs, _MODEL.Macro_oprts ,i);
        }

        _MODEL.Calculate_Local_Obs_for_States_to_Look(_LANCZOS,_BASIS);



        if(Do_Dynamics){

            _MODEL.Read_parameters_for_dynamics(inp_filename);

            _MODEL.Initialize_Opr_for_Dynamics(_BASIS);

            LANCZOS _LANCZOS_Dynamics;
            _LANCZOS_Dynamics.Dynamics_performed=true;
            _LANCZOS_Dynamics.Read_Lanczos_parameters(inp_filename);
            _LANCZOS_Dynamics.Eig_vec=_LANCZOS.Eig_vec;
            _LANCZOS_Dynamics.GS_energy=_LANCZOS.GS_energy;
            _LANCZOS_Dynamics.Get_Dynamics_seed(_MODEL.Dyn_opr);

            _LANCZOS_Dynamics.Perform_LANCZOS(_MODEL.Hamil);

        }


        if( (_BASIS.Ndn==_BASIS.Ndn)
                &&
                (_BASIS.Nup== _BASIS.Length
                 ))
        {
            if(PERFORM_VARIATIONAL_STATE_ANALYSIS==true){
                _MODEL.Get_Variational_State(_BASIS, _BASIS.Length);
                cout<<"Overlap b/w Ansatz from input and GS = "<<dot_product(_LANCZOS.Eig_vec, _MODEL.State_)<<endl;

                if(!Read_Half_filling_Variational_states){
                    _MODEL.Get_overlap_matrix_for_Anzatz_basis(_BASIS,_LANCZOS.Eig_vec);

                    cout<<scientific<<setprecision(4)<<endl;
                    cout<<"+++++++++OVERLAP MATRIX+++++++++++++"<<endl;
                    Print_Matrix(_MODEL.overlap_matrix_for_Anzatz_basis);
                    cout<<"++++++++++++++++++++++++++++++++++++"<<endl;
                    _MODEL.Create_OS_TS_states_by_reading(_BASIS,_LANCZOS.Eig_vec);
                }
                else{
                    _MODEL.Read_Anzatz_basis(_BASIS,_LANCZOS.Eig_vec);
                    _MODEL.Create_OS_TS_states_by_reading(_BASIS,_LANCZOS.Eig_vec);
                }

                cout<<scientific<<setprecision(10)<<endl;
                _MODEL.Variational_state_optimization(_BASIS,_LANCZOS.Eig_vec);

                cout<<"XXXXXXXXXXXXX S.S for Variational stateXXXXXXXXXXXXXXXX"<<endl;
                // _LANCZOS.Measure_two_point_observables(_MODEL.two_point_obs, _MODEL.Two_point_oprts, _BASIS.Length,  _MODEL.State_ , _MODEL.PBC);
                cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;

                cout<<"Orbital symmetry of Ansatz for half-filling:"<<endl;
                _MODEL.Check_orbital_symmetry(_BASIS, _MODEL.State_);

            }



            //Need N-2 Sector for this part of calculation
            if(PERFORM_Nm2_Variational_Calculation==true){

                BASIS_2_orb_Hubb_chain _BASIS_Nm2;
                MODEL_2_orb_Hubb_chain _MODEL_Nm2;
                _MODEL_Nm2.Read_parameters(_BASIS_Nm2, inp_filename);


                //In Sz=0 sector
                _BASIS_Nm2.Nup = _BASIS_Nm2.Nup -1;
                _BASIS_Nm2.Ndn = _BASIS_Nm2.Ndn -1;
                _BASIS_Nm2.Construct_basis();


                //XXXXXXXXX----GS for N-2 sector---XXXXXXXXXXXXXXXXXXXXXX
                _MODEL_Nm2.Add_diagonal_terms(_BASIS_Nm2);
                _MODEL_Nm2.Add_non_diagonal_terms(_BASIS_Nm2);
                _MODEL_Nm2.Add_connections(_BASIS_Nm2);
                cout<<"Size of Hilbert space = "<<_MODEL_Nm2.Hamil.nrows<<endl;
                cout<<scientific<<setprecision(8);
                //Print_Matrix_COO(_MODEL.Hamil);


                LANCZOS _LANCZOS_Nm2;
                _LANCZOS_Nm2.Dynamics_performed=false;
                _LANCZOS_Nm2.Read_Lanczos_parameters(inp_filename);
                _LANCZOS_Nm2.Perform_LANCZOS(_MODEL_Nm2.Hamil);

                cout<<"Orbital symmetry of Ground state for (half-filling - 2):"<<endl;
                _MODEL.Check_orbital_symmetry(_BASIS_Nm2, _LANCZOS_Nm2.Eig_vec);


                _MODEL_Nm2.Initialize_two_point_to_calculate(_BASIS_Nm2);
                for(int i=0;i<_LANCZOS_Nm2.states_to_look.size();i++){
                    _LANCZOS_Nm2.Measure_two_point_observables(_MODEL_Nm2.two_point_obs, _MODEL_Nm2.Two_point_oprts, _BASIS_Nm2.Length, i, _MODEL_Nm2.PBC);
                }



                //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                cout<<"XXXXX Create_Anzatz_basis_2Holes XXXXXXXXXXXXXXXXXXX"<<endl;
                _MODEL_Nm2.Create_Anzatz_basis_2Holes(_BASIS_Nm2, _LANCZOS_Nm2.Eig_vec);
                assert(false);
                cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
                //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                cout<<"here"<<endl;
                cout<<"XXXXX Optimize_Anzatz_2Holes XXXXXXXXXXXXXXXXXXX"<<endl;
                Mat_1_doub OPT_2HOLE_RVB;
                _MODEL_Nm2.Optimize_Anzatz_basis_2Holes(_BASIS_Nm2, OPT_2HOLE_RVB, "READ_OVERLAPS_NOT");
                cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
                cout<<"<GS(N-2)|OPT_VEC> = "<<dot_product(_LANCZOS_Nm2.Eig_vec, OPT_2HOLE_RVB )<<endl;
                cout<<"here 2"<<endl;



                Mat_1_trio_int Hole_positions;
                Hole_positions.resize(2);
                double prob_;

                Hole_positions[0].orb_ = 0; Hole_positions[0].spin_=0; Hole_positions[0].site_=2;




                cout<<"Normalized probability of holes in |OPT_2_HOLES>, with one hole at:"<<
                      Hole_positions[0].orb_<<"  "<<Hole_positions[0].site_<<"  "<<endl;


                for(int orb_val=0;orb_val<2;orb_val++){
                    for(int site_val=0;site_val<_BASIS.Length;site_val++){
                        Hole_positions[1].orb_ = orb_val; Hole_positions[1].site_=site_val;

                        prob_ = _MODEL_Nm2.Get_Holes_Projected_state_probability(_BASIS_Nm2, Hole_positions, OPT_2HOLE_RVB);
                        cout<< Hole_positions[1].orb_<<"  "<<Hole_positions[1].site_<<"   = "<< prob_<<endl;

                    }
                }
                cout<<endl;

                cout<<"2-holes are fixed at:"<<endl;
                cout<<Hole_positions[0].orb_<<"  "<<Hole_positions[0].site_<<"  "<<endl;

                Hole_positions[1].orb_ = 1; Hole_positions[1].site_=3;
                cout<<Hole_positions[1].orb_<<"  "<<Hole_positions[1].site_<<"  "<<endl;


                Mat_1_doub vec_2_holes_, vec_2_holes_temp_;
                vec_2_holes_ = _MODEL_Nm2.Get_2Holes_Projected_state(_BASIS_Nm2, Hole_positions, OPT_2HOLE_RVB);

                cout<<"XXXXXXXXX SS for |OPT_2_HOLES> XXXXXXXXXXXXXXXXX"<<endl;
                _LANCZOS_Nm2.Measure_two_point_observables(_MODEL_Nm2.two_point_obs, _MODEL_Nm2.Two_point_oprts,
                                                           _BASIS_Nm2.Length, vec_2_holes_ , _MODEL_Nm2.PBC);
                cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;





                cout<<"############FOR N-2 GS####################"<<endl;
                Hole_positions[0].orb_ = 0; Hole_positions[0].spin_=0; Hole_positions[0].site_=2;


                cout<<"Normalized probability of holes in |GS(N-2)>, with one hole at:"<<
                      Hole_positions[0].orb_<<"  "<<Hole_positions[0].site_<<"  "<<endl;


                for(int orb_val=0;orb_val<2;orb_val++){
                    for(int site_val=0;site_val<_BASIS.Length;site_val++){
                        Hole_positions[1].orb_ = orb_val; Hole_positions[1].site_=site_val;

                        prob_ = _MODEL_Nm2.Get_Holes_Projected_state_probability(_BASIS_Nm2, Hole_positions, _LANCZOS_Nm2.Eig_vec);
                        cout<< Hole_positions[1].orb_<<"  "<<Hole_positions[1].site_<<"   = "<< prob_<<endl;

                    }
                }
                cout<<endl;

                cout<<"2-holes are fixed at:"<<endl;
                cout<<Hole_positions[0].orb_<<"  "<<Hole_positions[0].site_<<"  "<<endl;

                Hole_positions[1].orb_ = 1; Hole_positions[1].site_=3;
                cout<<Hole_positions[1].orb_<<"  "<<Hole_positions[1].site_<<"  "<<endl;


                vec_2_holes_ = _MODEL_Nm2.Get_2Holes_Projected_state(_BASIS_Nm2, Hole_positions, _LANCZOS_Nm2.Eig_vec);

                cout<<"XXXXXXXXX SS for |GS(N-2)> XXXXXXXXXXXXXXXXX"<<endl;
                _LANCZOS_Nm2.Measure_two_point_observables(_MODEL_Nm2.two_point_obs, _MODEL_Nm2.Two_point_oprts,
                                                           _BASIS_Nm2.Length, vec_2_holes_ , _MODEL_Nm2.PBC);
                cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;


                cout<<"######################################"<<endl;




                assert(false);


                //here
                double overlap_;
                double norm_;
                double alpha_S, alpha_O;
                Mat_1_doub vec1, vec2, vec_total, vec_total2;
                trio_int trio_0, trio_1;

                vec_total.clear();
                vec_total.resize(_BASIS_Nm2.D_dn_basis.size()*_BASIS_Nm2.D_up_basis.size());
                vec_total2.clear();
                vec_total2.resize(_BASIS_Nm2.D_dn_basis.size()*_BASIS_Nm2.D_up_basis.size());

                //------------------------
                for(int site=0;site<_BASIS.Length;site++){
                    // for(int site=0;site<1;site++){
                    int site_p1= site +1;
                    if(site_p1==_BASIS.Length){
                        site_p1=0;
                    }
                    //----------Spin-Singlet, Orbital Antisymmetric pair Anhilation From Variational state------------//
                    alpha_S =-1.0;alpha_O=-1.0;

                    if(PERFORM_VARIATIONAL_STATE_ANALYSIS==true){


                        //ALONG THE DIAGONAL
                        trio_0.orb_ = 1; trio_0.spin_=1; trio_0.site_=site_p1;
                        trio_1.orb_ = 0; trio_1.spin_=0; trio_1.site_=site;
                        _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 1.0*alpha_S);
                        Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _MODEL.State_, vec1);

                        trio_0.orb_ = 1; trio_0.spin_=0; trio_0.site_=site_p1;
                        trio_1.orb_ = 0; trio_1.spin_=1; trio_1.site_=site;
                        _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 1.0*1.0);
                        Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _MODEL.State_, vec2);
                        Subtract(vec1, -1.0, vec2);


                        trio_0.orb_ = 0; trio_0.spin_=1; trio_0.site_=site_p1;
                        trio_1.orb_ = 1; trio_1.spin_=0; trio_1.site_=site;
                        _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 1.0*alpha_S*alpha_O);
                        Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _MODEL.State_, vec2);
                        Subtract(vec1, -1.0, vec2);

                        trio_0.orb_ = 0; trio_0.spin_=0; trio_0.site_=site_p1;
                        trio_1.orb_ = 1; trio_1.spin_=1; trio_1.site_=site;
                        _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 1.0*alpha_O);
                        Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _MODEL.State_, vec2);
                        Subtract(vec1, -1.0, vec2);


                        //ALONG THE CHAIN
                        trio_0.orb_ = 0; trio_0.spin_=1; trio_0.site_=site_p1;
                        trio_1.orb_ = 0; trio_1.spin_=0; trio_1.site_=site;
                        _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 0.0*alpha_S);
                        Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _MODEL.State_, vec2);
                        Subtract(vec1, -1.0, vec2);


                        trio_0.orb_ = 0; trio_0.spin_=0; trio_0.site_=site_p1;
                        trio_1.orb_ = 0; trio_1.spin_=1; trio_1.site_=site;
                        _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 0.0*1.0);
                        Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _MODEL.State_, vec2);
                        Subtract(vec1, -1.0, vec2);


                        trio_0.orb_ = 1; trio_0.spin_=1; trio_0.site_=site_p1;
                        trio_1.orb_ = 1; trio_1.spin_=0; trio_1.site_=site;
                        _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 0.0*alpha_S*alpha_O);
                        Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _MODEL.State_, vec2);
                        Subtract(vec1, -1.0, vec2);

                        trio_0.orb_ = 1; trio_0.spin_=0; trio_0.site_=site_p1;
                        trio_1.orb_ = 1; trio_1.spin_=1; trio_1.site_=site;
                        _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 0.0*alpha_O);
                        Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _MODEL.State_, vec2);
                        Subtract(vec1, -1.0, vec2);




                        //----------------------------------------//

                        Subtract(vec_total, -1.0, vec1);
                    }





                    //----------Spin-Singlet, Orbital Antisymmetric pair  Anhilation From GS(N)------------//


                    //ALONG THE DIAGONAL
                    trio_0.orb_ = 1; trio_0.spin_=1; trio_0.site_=site_p1;
                    trio_1.orb_ = 0; trio_1.spin_=0; trio_1.site_=site;
                    _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 1.0*alpha_S);
                    Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _LANCZOS.Eig_vec, vec1);

                    trio_0.orb_ = 1; trio_0.spin_=0; trio_0.site_=site_p1;
                    trio_1.orb_ = 0; trio_1.spin_=1; trio_1.site_=site;
                    _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 1.0*1.0);
                    Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _LANCZOS.Eig_vec, vec2);
                    Subtract(vec1, -1.0, vec2);


                    trio_0.orb_ = 0; trio_0.spin_=1; trio_0.site_=site_p1;
                    trio_1.orb_ = 1; trio_1.spin_=0; trio_1.site_=site;
                    _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 1.0*alpha_S*alpha_O);
                    Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _LANCZOS.Eig_vec, vec2);
                    Subtract(vec1, -1.0, vec2);

                    trio_0.orb_ = 0; trio_0.spin_=0; trio_0.site_=site_p1;
                    trio_1.orb_ = 1; trio_1.spin_=1; trio_1.site_=site;
                    _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 1.0*alpha_O);
                    Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _LANCZOS.Eig_vec, vec2);
                    Subtract(vec1, -1.0, vec2);


                    //ALONG THE CHAIN
                    trio_0.orb_ = 0; trio_0.spin_=1; trio_0.site_=site_p1;
                    trio_1.orb_ = 0; trio_1.spin_=0; trio_1.site_=site;
                    _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 0.0*alpha_S);
                    Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _LANCZOS.Eig_vec, vec2);
                    Subtract(vec1, -1.0, vec2);

                    trio_0.orb_ = 0; trio_0.spin_=0; trio_0.site_=site_p1;
                    trio_1.orb_ = 0; trio_1.spin_=1; trio_1.site_=site;
                    _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 0.0*1.0);
                    Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _LANCZOS.Eig_vec, vec2);
                    Subtract(vec1, -1.0, vec2);


                    trio_0.orb_ = 1; trio_0.spin_=1; trio_0.site_=site_p1;
                    trio_1.orb_ = 1; trio_1.spin_=0; trio_1.site_=site;
                    _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 0.0*alpha_S*alpha_O);
                    Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _LANCZOS.Eig_vec, vec2);
                    Subtract(vec1, -1.0, vec2);

                    trio_0.orb_ = 1; trio_0.spin_=0; trio_0.site_=site_p1;
                    trio_1.orb_ = 1; trio_1.spin_=1; trio_1.site_=site;
                    _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 0.0*alpha_O);
                    Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _LANCZOS.Eig_vec, vec2);
                    Subtract(vec1, -1.0, vec2);



                    //ALONG THE RUNG
                    trio_0.orb_ = 1; trio_0.spin_=1; trio_0.site_=site;
                    trio_1.orb_ = 0; trio_1.spin_=0; trio_1.site_=site;
                    _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 0.0*alpha_S);
                    Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _LANCZOS.Eig_vec, vec2);
                    Subtract(vec1, -1.0, vec2);

                    trio_0.orb_ = 1; trio_0.spin_=0; trio_0.site_=site;
                    trio_1.orb_ = 0; trio_1.spin_=1; trio_1.site_=site;
                    _MODEL_Nm2.Get_Pair_Operator_Matrix(_BASIS, _BASIS_Nm2, trio_0, trio_1, 0.0*1.0);
                    Matrix_COO_vector_multiplication("full",_MODEL_Nm2.Pair_Annihilation , _LANCZOS.Eig_vec, vec2);
                    Subtract(vec1, -1.0, vec2);

                    //----------------------------------------//

                    Subtract(vec_total2, -1.0, vec1);


                }
                //-------------------







                cout<<"###############GET 2-hole PROJECTED STATE PHYSICS####################"<<endl;

                cout<<"2-holes are fixed at:"<<endl;
                cout<<Hole_positions[0].orb_<<"  "<<Hole_positions[0].site_<<"  "<<endl;

                Hole_positions[1].orb_ = 1; Hole_positions[1].site_=4;
                cout<<Hole_positions[1].orb_<<"  "<<Hole_positions[1].site_<<"  "<<endl;


                Mat_1_doub vec_2_holes, vec_2_holes_temp;
                vec_2_holes = _MODEL_Nm2.Get_2Holes_Projected_state(_BASIS_Nm2, Hole_positions, _LANCZOS_Nm2.Eig_vec);





                cout<<"Normalized probability of holes in vec_2_holes, with one hole at:"<<
                      Hole_positions[0].orb_<<"  "<<Hole_positions[0].site_<<"  "<<endl;


                vec_2_holes_temp=vec_2_holes;
                for(int orb_val=0;orb_val<2;orb_val++){
                    for(int site_val=0;site_val<_BASIS.Length;site_val++){
                        Hole_positions[1].orb_ = orb_val; Hole_positions[1].site_=site_val;

                        prob_ = _MODEL_Nm2.Get_Holes_Projected_state_probability(_BASIS_Nm2, Hole_positions, vec_2_holes_temp);
                        cout<< Hole_positions[1].orb_<<"  "<<Hole_positions[1].site_<<"   = "<< prob_<<endl;

                    }
                }
                cout<<endl;


                cout<<"XXXXXXXXX SS for H[Delta|Ansatz(N)>] XXXXXXXXXXXXXXXXX"<<endl;
                _LANCZOS_Nm2.Measure_two_point_observables(_MODEL_Nm2.two_point_obs, _MODEL_Nm2.Two_point_oprts,
                                                           _BASIS_Nm2.Length, vec_2_holes , _MODEL_Nm2.PBC);
                cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;

                cout<<endl;
                cout<<"###################################################"<<endl;





                if(PERFORM_VARIATIONAL_STATE_ANALYSIS==true){

                    norm_=dot_product(vec_total, vec_total);
                    cout <<"Norm(Delta|Ansatz(N)>)  = "<<norm_<<endl;
                    overlap_ = dot_product(_LANCZOS_Nm2.Eig_vec, vec_total);
                    cout <<"<GS(N-2)|Delta|Ansatz(N)>  = "<<(overlap_)<<endl;

                    cout<<"Orbital symmetry of Delta|Ansatz(N)>:"<<endl;
                    _MODEL.Check_orbital_symmetry(_BASIS_Nm2, vec_total);

                    cout<<"XXXXXXXXX SS for Delta|Ansatz(N)>XXXXXXXXXXXXXXXXX"<<endl;
                    _LANCZOS_Nm2.Measure_two_point_observables(_MODEL_Nm2.two_point_obs, _MODEL_Nm2.Two_point_oprts, _BASIS_Nm2.Length, vec_total , _MODEL_Nm2.PBC);
                    cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;


                    /*###############GET PROJECTED STATE####################*/

                    cout<<"Normalized probability of holes in Delta|Ansatz(N)>, with one hole at:"<<
                          Hole_positions[0].orb_<<"  "<<Hole_positions[0].site_<<"  "<<endl;


                    for(int orb_val=0;orb_val<2;orb_val++){
                        for(int site_val=0;site_val<_BASIS.Length;site_val++){
                            Hole_positions[1].orb_ = orb_val; Hole_positions[1].site_=site_val;

                            prob_ = _MODEL_Nm2.Get_Holes_Projected_state_probability(_BASIS_Nm2, Hole_positions, vec_total);
                            cout<< Hole_positions[1].orb_<<"  "<<Hole_positions[1].site_<<"   = "<< prob_<<endl;

                        }
                    }
                    cout<<endl;
                    /*###################################################*/


                    cout<<"#################################################################################"<<endl;
                    cout<<"IMPROVING Delta|Ansatz(N)> state by Acting Hamiltonian 1 time:"<<endl;
                    Mat_1_doub vec_total_updated;
                    double norm_updated_;
                    vec_total_updated.clear();
                    vec_total_updated.resize(_BASIS_Nm2.D_dn_basis.size()*_BASIS_Nm2.D_up_basis.size());
                    Matrix_COO_vector_multiplication("U", _MODEL_Nm2.Hamil, vec_total, vec_total_updated);

                    norm_updated_=dot_product(vec_total_updated, vec_total_updated);
                    cout <<"Norm(H[Delta|Ansatz(N)>])  = "<<norm_updated_<<endl;
                    overlap_ = dot_product(_LANCZOS_Nm2.Eig_vec, vec_total_updated);
                    cout <<"<GS(N-2)|H [Delta|Ansatz(N)>]  = "<<(overlap_)<<endl;

                    cout<<"Orbital symmetry of H[Delta|Ansatz(N)>]:"<<endl;
                    _MODEL.Check_orbital_symmetry(_BASIS_Nm2, vec_total_updated);

                    cout<<"XXXXXXXXX SS for H[Delta|Ansatz(N)>] XXXXXXXXXXXXXXXXX"<<endl;
                    _LANCZOS_Nm2.Measure_two_point_observables(_MODEL_Nm2.two_point_obs, _MODEL_Nm2.Two_point_oprts, _BASIS_Nm2.Length, vec_total_updated , _MODEL_Nm2.PBC);
                    cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;


                    /*###############GET PROJECTED STATE####################*/

                    cout<<"Normalized probability of holes in H[Delta|Ansatz(N)>], with one hole at:"<<
                          Hole_positions[0].orb_<<"  "<<Hole_positions[0].site_<<"  "<<endl;


                    for(int orb_val=0;orb_val<2;orb_val++){
                        for(int site_val=0;site_val<_BASIS.Length;site_val++){
                            Hole_positions[1].orb_ = orb_val; Hole_positions[1].site_=site_val;

                            prob_ = _MODEL_Nm2.Get_Holes_Projected_state_probability(_BASIS_Nm2, Hole_positions, vec_total_updated);
                            cout<< Hole_positions[1].orb_<<"  "<<Hole_positions[1].site_<<"   = "<< prob_<<endl;

                        }
                    }
                    cout<<endl;
                    /*###################################################*/


                    cout<<"#################################################################################"<<endl;


                    cout<<"#################################################################################"<<endl;
                    cout<<"IMPROVING: Discarded double occupancies from H[Delta|Ansatz(N)>]:"<<endl;

                    _MODEL.Discard_double_occupancies(_BASIS_Nm2, vec_total_updated);
                    norm_updated_=dot_product(vec_total_updated, vec_total_updated);
                    cout <<"Norm(Dis. H[Delta|Ansatz(N)>])  = "<<norm_updated_<<endl;
                    overlap_ = dot_product(_LANCZOS_Nm2.Eig_vec, vec_total_updated);
                    cout <<"<GS(N-2)|Dis.H [Delta|Ansatz(N)>]  = "<<(overlap_)<<endl;

                    cout<<"Orbital symmetry of Dis.H[Delta|Ansatz(N)>]:"<<endl;
                    _MODEL.Check_orbital_symmetry(_BASIS_Nm2, vec_total_updated);

                    cout<<"XXXXXXXXX SS for Dis.H[Delta|Ansatz(N)>] XXXXXXXXXXXXXXXXX"<<endl;
                    _LANCZOS_Nm2.Measure_two_point_observables(_MODEL_Nm2.two_point_obs, _MODEL_Nm2.Two_point_oprts, _BASIS_Nm2.Length, vec_total_updated , _MODEL_Nm2.PBC);
                    cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;


                    /*###############GET PROJECTED STATE####################*/

                    cout<<"Normalized probability of holes in Dis.H[Delta|Ansatz(N)>], with one hole at:"<<
                          Hole_positions[0].orb_<<"  "<<Hole_positions[0].site_<<"  "<<endl;


                    for(int orb_val=0;orb_val<2;orb_val++){
                        for(int site_val=0;site_val<_BASIS.Length;site_val++){
                            Hole_positions[1].orb_ = orb_val; Hole_positions[1].site_=site_val;

                            prob_ = _MODEL_Nm2.Get_Holes_Projected_state_probability(_BASIS_Nm2, Hole_positions, vec_total_updated);
                            cout<< Hole_positions[1].orb_<<"  "<<Hole_positions[1].site_<<"   = "<< prob_<<endl;

                        }
                    }
                    cout<<endl;
                    /*###################################################*/


                    cout<<"#################################################################################"<<endl;






                    cout<<"#################################################################################"<<endl;
                    cout<<"IMPROVING Delta|Ansatz(N)> using beta|old> + alpha discarded(H|old>):"<<endl;

                    double VALUE_CHECK_OLD, VALUE_CHECK;
                    double beta, beta_final;
                    double alpha, alpha_final;

                    Mat_1_doub vec_temp;
                    vec_temp.clear();
                    vec_total_updated.clear();
                    vec_total_updated.resize(_BASIS_Nm2.D_dn_basis.size()*_BASIS_Nm2.D_up_basis.size());
                    Matrix_COO_vector_multiplication("U", _MODEL_Nm2.Hamil, vec_total, vec_total_updated);

                    _MODEL.Discard_double_occupancies(_BASIS_Nm2, vec_total_updated);

                    VALUE_CHECK_OLD =0.0;



                    for(alpha=1.0;alpha<=1.0;alpha=alpha+0.01){
                        for(beta=0;beta<=0;beta=beta+0.01){

                            vec_temp = vec_total_updated;
                            Sum(vec_temp, alpha, vec_total, beta);


                            norm_updated_=dot_product(vec_temp, vec_temp);
                            overlap_ = dot_product(_LANCZOS_Nm2.Eig_vec, vec_temp);
                            VALUE_CHECK = (overlap_*overlap_)/(norm_updated_);

                            if( VALUE_CHECK > VALUE_CHECK_OLD){
                                VALUE_CHECK_OLD=VALUE_CHECK;
                                alpha_final=alpha;
                                beta_final = beta;
                            }

                            cout<<"alpha = "<<alpha<<", beta = "<<beta<<", value = "<<VALUE_CHECK<<endl;

                        }
                    }

                    cout<<"FINAL ALPHA = "<<alpha_final<<endl;
                    cout<<"FINAL BETA = "<<beta_final<<endl;



                    Sum(vec_total_updated, alpha_final, vec_total, beta_final);

                    norm_updated_=dot_product(vec_total_updated, vec_total_updated);
                    cout <<"Norm(a*Dis.H*Delta|Ansatz(N)>)+ b*Delta|Ansatz(N)>)  = "<<norm_updated_<<endl;
                    overlap_ = dot_product(_LANCZOS_Nm2.Eig_vec, vec_total_updated);
                    cout <<"<GS(N-2)|a*Dis.H[Delta|Ansatz(N)> + b*Delta|Ansatz(N)>]  = "<<(overlap_)<<endl;


                    cout<<"Orbital symmetry of a*Dis.H[Delta|Ansatz(N)> + b*Delta|Ansatz(N)>:"<<endl;
                    _MODEL.Check_orbital_symmetry(_BASIS_Nm2, vec_total_updated);

                    cout<<"XXXXXXXXX SS for a*Dis.H[Delta|Ansatz(N)> + b*Delta|Ansatz(N)> XXXXXXXXXXXXXXXXX"<<endl;
                    _LANCZOS_Nm2.Measure_two_point_observables(_MODEL_Nm2.two_point_obs, _MODEL_Nm2.Two_point_oprts, _BASIS_Nm2.Length, vec_total_updated , _MODEL_Nm2.PBC);
                    cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;


                    /*###############GET PROJECTED STATE####################*/

                    cout<<"Normalized probability of holes in a*Dis.H[Delta|Ansatz(N)> + b*Delta|Ansatz(N)>, with one hole at:"<<
                          Hole_positions[0].orb_<<"  "<<Hole_positions[0].site_<<"  "<<endl;


                    for(int orb_val=0;orb_val<2;orb_val++){
                        for(int site_val=0;site_val<_BASIS.Length;site_val++){
                            Hole_positions[1].orb_ = orb_val; Hole_positions[1].site_=site_val;

                            prob_ = _MODEL_Nm2.Get_Holes_Projected_state_probability(_BASIS_Nm2, Hole_positions, vec_total_updated);
                            cout<< Hole_positions[1].orb_<<"  "<<Hole_positions[1].site_<<"   = "<< prob_<<endl;

                        }
                    }
                    cout<<endl;
                    /*###################################################*/


                    cout<<"#################################################################################"<<endl;



                    cout<<"#################################################################################"<<endl;
                    cout<<"IMPROVING Delta|Ansatz(N)> using internal coefficients, "
                       <<"State(3_coeff) = a0*Holes(1)[Delta|Ansatz(N)>] + a1*Holes(5)[Delta|Ansatz(N)>] + beta*[Delta|Ansatz(N)>] "
                      <<endl;

                    vec_temp.clear();
                    vec_total_updated.clear();
                    vec_total_updated.resize(_BASIS_Nm2.D_dn_basis.size()*_BASIS_Nm2.D_up_basis.size());




                    VALUE_CHECK_OLD=0.0;

                    Mat_1_doub Alpha_;
                    Alpha_.resize(2);

                    Mat_1_doub Alpha_final;
                    Alpha_final.resize(2);

                    Alpha_[0]=0.0;
                    while(Alpha_[0]<=0.0){

                        Alpha_[1]=1.0;
                        while(Alpha_[1]<=1.0){

                            //Alpha_[2]=-1.0;
                            //while(Alpha_[2]<=1.0){

                            //  Alpha_[3]=-1.0;
                            //  while(Alpha_[3]<=1.0){

                            if( true

                                    //(Alpha_[3]==Alpha_[1] )
                                    //&&
                                    //(Alpha_[2]==Alpha_[0] )

                                    ){

                                vec_total_updated=vec_total;
                                _MODEL_Nm2.Create_states_with_hole_hopping(_BASIS_Nm2, Alpha_, vec_total_updated);

                                //cout<<"here 1"<<endl;
                                for(beta=0.0;beta<=4;beta=beta+10000.0){

                                    vec_temp = vec_total_updated;
                                    Sum(vec_temp, 1.0, vec_total, beta);


                                    norm_updated_=dot_product(vec_temp, vec_temp);
                                    // cout <<"Norm(H[Delta|Ansatz(N)>])  = "<<norm_updated_<<endl;
                                    overlap_ = dot_product(_LANCZOS_Nm2.Eig_vec, vec_temp);
                                    // cout <<"<GS(N-2)|H [Delta|Ansatz(N)>]  = "<<(overlap_)<<endl;
                                    VALUE_CHECK = (overlap_*overlap_)/(norm_updated_);

                                    if( VALUE_CHECK > VALUE_CHECK_OLD){
                                        VALUE_CHECK_OLD=VALUE_CHECK;
                                        Alpha_final=Alpha_;
                                        beta_final = beta;
                                    }

                                    for(int i=0;i<Alpha_.size();i++){
                                        cout<<"Alpha_["<<i<<"] = "<<Alpha_[i]<<", ";
                                    }
                                    cout<<"beta = "<<beta<<", value = "<<VALUE_CHECK<<endl;

                                }

                            }
                            //   Alpha_[3]=Alpha_[3]+0.5;
                            // }
                            //  Alpha_[2]=Alpha_[2]+0.5;
                            // }
                            Alpha_[1]=Alpha_[1]+1;
                        }
                        Alpha_[0]=Alpha_[0]+1;
                    }

                    for(int i=0;i<Alpha_final.size();i++){
                        cout<<"FINAL ALPHA["<<i<<"] = "<<Alpha_final[i]<<endl;
                    }
                    cout<<"FINAL BETA = "<<beta_final<<endl;



                    vec_total_updated=vec_total;
                    _MODEL_Nm2.Create_states_with_hole_hopping(_BASIS_Nm2, Alpha_final, vec_total_updated);
                    Sum(vec_total_updated, 1.0, vec_total, beta_final);


                    norm_updated_=dot_product(vec_total_updated, vec_total_updated);
                    cout <<"Norm(State(3_coeff)])  = "<<norm_updated_<<endl;
                    overlap_ = dot_product(_LANCZOS_Nm2.Eig_vec, vec_total_updated);
                    cout <<"<GS(N-2)|State(3_coeff)]  = "<<(overlap_)<<endl;


                    cout<<"Orbital symmetry of State(3_coeff):"<<endl;
                    _MODEL.Check_orbital_symmetry(_BASIS_Nm2, vec_total_updated);

                    cout<<"XXXXXXXXX SS for State(3_coeff) XXXXXXXXXXXXXXXXX"<<endl;
                    _LANCZOS_Nm2.Measure_two_point_observables(_MODEL_Nm2.two_point_obs, _MODEL_Nm2.Two_point_oprts, _BASIS_Nm2.Length, vec_total_updated , _MODEL_Nm2.PBC);
                    cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;


                    /*###############GET PROJECTED STATE####################*/

                    cout<<"Normalized probability of holes in State(3_coeff), with one hole at:"<<
                          Hole_positions[0].orb_<<"  "<<Hole_positions[0].site_<<"  "<<endl;


                    for(int orb_val=0;orb_val<2;orb_val++){
                        for(int site_val=0;site_val<_BASIS.Length;site_val++){
                            Hole_positions[1].orb_ = orb_val; Hole_positions[1].site_=site_val;

                            prob_ = _MODEL_Nm2.Get_Holes_Projected_state_probability(_BASIS_Nm2, Hole_positions, vec_total_updated);
                            cout<< Hole_positions[1].orb_<<"  "<<Hole_positions[1].site_<<"   = "<< prob_<<endl;

                        }
                    }
                    cout<<endl;
                    /*###################################################*/


                    cout<<"#################################################################################"<<endl;
                    cout<<"Getting hole projections for [Delta|Ansatz>] "<<endl;

                    Mat_1_doub State_class_resolved;
                    Mat_1_doub Projected_state;

                    _MODEL_Nm2.Create_STATES_OS_TS_DIAGONAL_HOLES(_BASIS, _BASIS_Nm2, _MODEL.STATES_OS_TS);

                    Mat_1_doub alpha_used;
                    alpha_used.resize(2);
                    alpha_used[0]=1.0;alpha_used[1]=1.0;

                    for(int ci=0;ci<_MODEL_Nm2.STATES_OS_TS_DIAGONAL_HOLES.size();ci++){

                        State_class_resolved=_MODEL_Nm2.STATES_OS_TS_DIAGONAL_HOLES[ci];


                        //Analysis on State with diagonal holes------------
                        vec_total_updated=State_class_resolved;
                        cout<<"Normalized probability of holes in STATES_OS_TS_DIAGONAL_HOLES["<<ci<<"], with one hole at:"<<
                              Hole_positions[0].orb_<<"  "<<Hole_positions[0].site_<<"  "<<endl;


                        for(int orb_val=0;orb_val<2;orb_val++){
                            for(int site_val=0;site_val<_BASIS.Length;site_val++){
                                Hole_positions[1].orb_ = orb_val; Hole_positions[1].site_=site_val;

                                prob_ = _MODEL_Nm2.Get_Holes_Projected_state_probability(_BASIS_Nm2, Hole_positions, vec_total_updated);
                                cout<< Hole_positions[1].orb_<<"  "<<Hole_positions[1].site_<<"   = "<< prob_<<endl;

                            }
                        }

                        cout<<"XXXXXXXXX SS for STATES_OS_TS_DIAGONAL_HOLES["<<ci<<"] XXXXXXXXXXXXXXXXX"<<endl;
                        Projected_state = State_class_resolved;
                        _MODEL_Nm2.Get_Holes_Projected_state(_BASIS_Nm2, Hole_positions, Projected_state);
                        _LANCZOS_Nm2.Measure_two_point_observables(_MODEL_Nm2.two_point_obs, _MODEL_Nm2.Two_point_oprts, _BASIS_Nm2.Length,Projected_state , _MODEL_Nm2.PBC);
                        cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;

                        //----------------------------------------------



                        //Analysis on State with rung and sqrt(5) holes------------

                        vec_total_updated=State_class_resolved;

                        _MODEL_Nm2.Create_states_with_hole_hopping(_BASIS_Nm2, alpha_used, vec_total_updated);


                        cout<<"Normalized probability of holes in (a0*hole[1] + a1*hole[5])|V_"<<ci<<">, with one hole at:"<<
                              Hole_positions[0].orb_<<"  "<<Hole_positions[0].site_<<"  "<<endl;


                        for(int orb_val=0;orb_val<2;orb_val++){
                            for(int site_val=0;site_val<_BASIS.Length;site_val++){
                                Hole_positions[1].orb_ = orb_val; Hole_positions[1].site_=site_val;

                                prob_ = _MODEL_Nm2.Get_Holes_Projected_state_probability(_BASIS_Nm2, Hole_positions, vec_total_updated);
                                cout<< Hole_positions[1].orb_<<"  "<<Hole_positions[1].site_<<"   = "<< prob_<<endl;

                            }
                        }


                        Projected_state = vec_total_updated;
                        _MODEL_Nm2.Get_Holes_Projected_state(_BASIS_Nm2, Hole_positions, Projected_state);

                        cout<<"XXXXXXXXX SS for [one hole fixed at 0,2](a0*hole[1] + a1*hole[5])|V_"<<ci<<"> XXXXXXXXXXXXXXXXX"<<endl;
                        _LANCZOS_Nm2.Measure_two_point_observables(_MODEL_Nm2.two_point_obs, _MODEL_Nm2.Two_point_oprts, _BASIS_Nm2.Length, Projected_state , _MODEL_Nm2.PBC);
                        cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;

                        //----------------------------------------------------


                    }

                    cout<<"#################################################################################"<<endl;


                    cout<<"#################################################################################"<<endl;


                    cout<<"IMPROVING sum_(i)(i=0..7) a^(i)_(0)*Holes(dis=1)[Delta_(D)|V_(i)(N)>] + a^(i)_(1)*Holes(dis=sqrt(5))[Delta_(D)|V_(i)(N)>] + beta(i)[Delta_(D)|V_(i)(N)>]"
                       <<endl;


                    _MODEL_Nm2.Optimize_overlap_Nm2_Variational_State(_BASIS, _BASIS_Nm2, _MODEL.STATES_OS_TS, _LANCZOS_Nm2.Eig_vec, vec_total_updated);


                    cout<<"Orbital symmetry of FULL_VAR[Delta|Ansatz(N)>]:"<<endl;
                    _MODEL.Check_orbital_symmetry(_BASIS_Nm2, vec_total_updated);

                    cout<<"XXXXXXXXX SS for FULL_VAR[Delta|Ansatz(N)>] XXXXXXXXXXXXXXXXX"<<endl;
                    _LANCZOS_Nm2.Measure_two_point_observables(_MODEL_Nm2.two_point_obs, _MODEL_Nm2.Two_point_oprts, _BASIS_Nm2.Length, vec_total_updated , _MODEL_Nm2.PBC);
                    cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;


                    //###############GET PROJECTED STATE####################

                    cout<<"Normalized probability of holes in FULL_VAR[Delta|Ansatz(N)>], with one hole at:"<<
                          Hole_positions[0].orb_<<"  "<<Hole_positions[0].site_<<"  "<<endl;


                    for(int orb_val=0;orb_val<2;orb_val++){
                        for(int site_val=0;site_val<_BASIS.Length;site_val++){
                            Hole_positions[1].orb_ = orb_val; Hole_positions[1].site_=site_val;

                            prob_ = _MODEL_Nm2.Get_Holes_Projected_state_probability(_BASIS_Nm2, Hole_positions, vec_total_updated);
                            cout<< Hole_positions[1].orb_<<"  "<<Hole_positions[1].site_<<"   = "<< prob_<<endl;

                        }
                    }
                    cout<<endl;
                    //###################################################

                    //*/
                    cout<<"#################################################################################"<<endl;


                }



                norm_=dot_product(vec_total2, vec_total2);
                cout <<"Norm(Delta|GS(N)>)  = "<<norm_<<endl;
                overlap_ = dot_product(_LANCZOS_Nm2.Eig_vec, vec_total2);
                cout <<"<GS(N-2)|Delta|GS(N)>  = "<<(overlap_)<<endl;
                cout<<"Orbital symmetry of Delta|GS(N)>:"<<endl;
                _MODEL.Check_orbital_symmetry(_BASIS_Nm2, vec_total2);

                cout<<"XXXXXXXXX SS for Delta|GS(N)>XXXXXXXXXXXXXXXXX"<<endl;
                _LANCZOS_Nm2.Measure_two_point_observables(_MODEL_Nm2.two_point_obs, _MODEL_Nm2.Two_point_oprts, _BASIS_Nm2.Length, vec_total2 , _MODEL_Nm2.PBC);
                cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;




                cout<<"Normalized probability of holes in |GS(N-2)>, with one hole at:"<<
                      Hole_positions[0].orb_<<"  "<<Hole_positions[0].site_<<"  "<<endl;


                for(int orb_val=0;orb_val<2;orb_val++){
                    for(int site_val=0;site_val<_BASIS.Length;site_val++){
                        Hole_positions[1].orb_ = orb_val; Hole_positions[1].site_=site_val;

                        prob_ = _MODEL_Nm2.Get_Holes_Projected_state_probability(_BASIS_Nm2, Hole_positions, _LANCZOS_Nm2.Eig_vec);
                        cout<< Hole_positions[1].orb_<<"  "<<Hole_positions[1].site_<<"   = "<< prob_<<endl;

                    }
                }
                cout<<endl;


                cout<<"Normalized probability of holes in |Delta|GS(N)>, with one hole at:"<<
                      Hole_positions[0].orb_<<"  "<<Hole_positions[0].site_<<"  "<<endl;

                for(int orb_val=0;orb_val<2;orb_val++){
                    for(int site_val=0;site_val<_BASIS.Length;site_val++){
                        Hole_positions[1].orb_ = orb_val; Hole_positions[1].spin_=1; Hole_positions[1].site_=site_val;

                        prob_ = _MODEL_Nm2.Get_Holes_Projected_state_probability(_BASIS_Nm2, Hole_positions, vec_total2);
                        cout<< Hole_positions[1].orb_<<"  "<<Hole_positions[1].site_<<"   = "<< prob_<<endl;

                    }
                }
                cout<<endl;


                cout<<"Normalized probability of holes in |GS(N)>, with one hole at:"<<
                      Hole_positions[0].orb_<<"  "<<Hole_positions[0].site_<<"  "<<endl;
                for(int orb_val=0;orb_val<2;orb_val++){
                    for(int site_val=0;site_val<_BASIS.Length;site_val++){
                        Hole_positions[1].orb_ = orb_val; Hole_positions[1].spin_=1; Hole_positions[1].site_=site_val;

                        prob_ = _MODEL.Get_Holes_Projected_state_probability(_BASIS, Hole_positions, _LANCZOS.Eig_vec);
                        cout<< Hole_positions[1].orb_<<"  "<<Hole_positions[1].site_<<"   = "<< prob_<<endl;

                    }
                }



            }


        }


    }

#endif

    //=======================================================================================================================================================================================================//
    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //=======================================================================================================================================================================================================//




    if (model_name=="3_orb_Hubbard_chain_GC" && Restricted_Basis==false) {


        bool Dynamics_SPDOS = true;
        bool Above_mu = true;
        bool Below_mu= true;

        BASIS_3_orb_Hubb_chain_GC _BASIS;
        MODEL_3_orb_Hubb_chain_GC<BASIS_3_orb_Hubb_chain_GC> _MODEL(_BASIS);


        //---------------------------------------

        _MODEL.Read_parameters(inp_filename);
        cout<<"Parameters done"<<endl;

        //----remove later-----//
        // BASIS_3_orb_Hubb_chain_GC_Restricted _BASIS_Restricted;
        //_MODEL.Read_parameters(inp_filename);
        //_BASIS_Restricted.Construct_basis();
        //---------------------------//


        _BASIS.Construct_basis();

        cout<<"Size of Hilbert space = "<<_BASIS.D_up_basis.size()<<endl;
        cout<<scientific<<setprecision(6);

        cout<<"Diagonal part started"<<endl;
        _MODEL.Add_diagonal_terms();
        cout<<"Diagonal part done"<<endl;

        cout<<"Non Diagonal part started"<<endl;
        _MODEL.Add_non_diagonal_terms();
        cout<<"Non Diagonal part done"<<endl;


        _MODEL.Add_Spin_Orbit_Coupling();

        cout<<"Connections started"<<endl;
        _MODEL.Add_connections();
        cout<<"Connections done"<<endl;


        cout<<"Size of Hilbert space = "<<_MODEL.Hamil.nrows<<endl;
        cout<<scientific<<setprecision(1);
        // Print_Matrix_COO(_MODEL.Hamil);
        cout<<scientific<<setprecision(6);


        LANCZOS _LANCZOS;
        _LANCZOS.Dynamics_performed=false;
        _LANCZOS.Read_Lanczos_parameters(inp_filename);
        _LANCZOS.Perform_LANCZOS(_MODEL.Hamil);
        _LANCZOS.Write_full_spectrum();
        Print_vector_in_file(_LANCZOS.Eig_vec,"seed_GS.txt");

        _MODEL.Initialize_one_point_to_calculate();
        _MODEL.Initialize_two_point_to_calculate();


#ifdef USE_COMPLEX

        if(false){
            //--UPDATING Lanczos Eig vecs with phase---------//
            complex<double> a_;
            //double _PI_ = 3.14159265358979;
            Mat_1_Complex_doub Phases;
            Phases.resize(15);

            //Only works for specific value U=20, H_mag = 0.000001, lambda_SOC = 0.002, Random_seed_generator = 321

            a_=complex <double>(-0.408147,0.00907139);Phases[0]=conj(a_);
            a_=complex <double>(0.432517,0.559401);Phases[1]=conj(a_);
            a_=complex <double>(-0.418029,0.27432);Phases[2]=conj(a_);
            a_=complex <double>(0.206927,0.286673);Phases[3]=conj(a_);
            a_=complex <double>(-0.406994,0.578235);Phases[4]=conj(a_);
            a_=complex <double>(0.0377325,-0.706099);Phases[5]=conj(a_);
            a_=complex <double>(0.205799,0.202436);Phases[6]=conj(a_);
            a_=complex <double>(0.259928,-0.239557);Phases[7]=conj(a_);
            a_=complex <double>(-0.281476,0.648669);Phases[8]=conj(a_);
            a_=complex <double>(0.342784,-0.364005);Phases[9]=conj(a_);
            a_=complex <double>(-0.499813,0.0117241);Phases[10]=conj(a_);
            a_=complex <double>(-0.014609,0.816366);Phases[11]=conj(a_);
            a_=complex <double>(-0.391111,0.31158);Phases[12]=conj(a_);
            a_=complex <double>(0.487954,0.10909);Phases[13]=conj(a_);
            a_=complex <double>(0.113232,0.566138);Phases[14]=conj(a_);

            for(int state_=0;state_<Phases.size();state_++){
                value_multiply_vector(Phases[state_], _LANCZOS.Eig_vecs[state_]);
                Normalize_vec(_LANCZOS.Eig_vecs[state_]);
            }

        }


        //-----------DONE--------------------------------//
#endif


        for(int i=0;i<_LANCZOS.states_to_look.size();i++){
            cout<<"===================FOR STATE NO "<<i<<"============================="<<endl;
            _LANCZOS.Measure_one_point_observables(_MODEL.one_point_obs, _MODEL.One_point_oprts, _BASIS.Length, i);
            // cout<<"Energy = "<<_LANCZOS.Evals_Tri_all[_LANCZOS.Evals_Tri_all.size()-1][i]<<endl;
            //_LANCZOS.Measure_two_point_observables(_MODEL.two_point_obs, _MODEL.Two_point_oprts, _BASIS.Length, i, _MODEL.PBC);
            _LANCZOS.Measure_two_point_observables_smartly(_MODEL.one_point_obs,_MODEL.One_point_oprts, _BASIS.Length, i,"3_orb_Hubbard_chain_GC");
            cout<<"============================================================================"<<endl;

            _LANCZOS.Measure_macro_observables(_MODEL.macro_obs, _MODEL.Macro_oprts, i);
        }

        //_MODEL.Get_Delta_Matrix(_LANCZOS);
        _MODEL.Get_ExcitonCoherence_Length(_LANCZOS.Eig_vecs[0]);

        for(int state_=0;state_<_LANCZOS.Eig_vecs.size();state_++){

            string state_string;
            stringstream state_ss;
            state_ss << state_;
            state_string = state_ss.str();

            string file_name_state = "State"+ state_string + "_.txt";

            Print_vector_in_file(_LANCZOS.Eig_vecs[state_], file_name_state);
        }



        //------------EXACT DIAGONALIZATION---------------------//




        Mat_1_real Eigen_ED;
        Mat_2_doub vecs;
        if(_MODEL.Hamil.nrows>800){
            DO_FULL_DIAGONALIZATION=false;
        }
        if(DO_FULL_DIAGONALIZATION==true){

            string fl_ED_out = "EXACT_RESULTS.txt";
            ofstream file_ED_out(fl_ED_out.c_str());
            cout<<"-----------------------------------------------------------------------"<<endl;
            cout<<"AFTER THIS EXACT DIAGONALIZATION RESULTS ARE WRITTEN IN EXACT_RESULTS.txt-------------------"<<endl;

            file_ED_out<<"//*****************Exact Diagonalization Energies**************//"<<endl;
            file_ED_out<<"//------------------------------------------------------------//"<<endl;
            Diagonalize(_MODEL.Hamil, Eigen_ED, vecs);
            for(int i=0;i<Eigen_ED.size();i++){
                file_ED_out<<i<<"   "<<scientific<<setprecision(6)<<Eigen_ED[i]<<endl;
            }

            Mat_2_doub Dummy_Eig_vecs;
            Dummy_Eig_vecs = _LANCZOS.Eig_vecs;

            _LANCZOS.Eig_vecs.clear();
            _LANCZOS.Eig_vecs.resize(_MODEL.Hamil.nrows);
            for(int s_=0;s_<_MODEL.Hamil.nrows;s_++){
                _LANCZOS.Eig_vecs[s_] = vecs[s_];
            }

                cout<<endl;
            for(int i=0;i<_MODEL.Hamil.nrows;i++){
                cout<<"===================FOR STATE NO "<<i<<"============================="<<endl;
                _LANCZOS.Measure_one_point_observables(_MODEL.one_point_obs, _MODEL.One_point_oprts, _BASIS.Length, i);
                // cout<<"Energy = "<<_LANCZOS.Evals_Tri_all[_LANCZOS.Evals_Tri_all.size()-1][i]<<endl;
                //_LANCZOS.Measure_two_point_observables(_MODEL.two_point_obs, _MODEL.Two_point_oprts, _BASIS.Length, i, _MODEL.PBC);
                _LANCZOS.Measure_two_point_observables_smartly(_MODEL.one_point_obs,_MODEL.One_point_oprts, _BASIS.Length, i,"3_orb_Hubbard_chain_GC");
                cout<<"============================================================================"<<endl;
                _LANCZOS.Measure_macro_observables(_MODEL.macro_obs, _MODEL.Macro_oprts, i);

                cout<<endl<<endl<<endl;
            }

            _LANCZOS.Eig_vecs = Dummy_Eig_vecs;


        }



        cout<<"-----------------------------------------------------------------------"<<endl;
        //------------------------------------------------------//


        // for(int i=0;i<19;i++){
        // _LANCZOS.Measure_two_point_observables(_MODEL.two_point_obs, _MODEL.Two_point_oprts, _BASIS.Length, 0, _MODEL.PBC);

        //_LANCZOS.Measure_two_point_observables_smartly(_MODEL.one_point_obs,_MODEL.One_point_oprts, _BASIS.Length, 0);
        // }


        //_MODEL.Calculate_Local_Obs_for_States_to_Look(_LANCZOS,_BASIS);





        if(Do_Dynamics && (!Dynamics_SPDOS)){

            _MODEL.Read_parameters_for_dynamics(inp_filename);

            _MODEL.Initialize_Opr_for_Dynamics(_LANCZOS);

            LANCZOS _LANCZOS_Dynamics;
            _LANCZOS_Dynamics.Dynamics_performed=true;
            _LANCZOS_Dynamics.Read_Lanczos_parameters(inp_filename);
            _LANCZOS_Dynamics.Eig_vec=_LANCZOS.Eig_vec;
            _LANCZOS_Dynamics.GS_energy=_LANCZOS.GS_energy;
            _LANCZOS_Dynamics.Get_Dynamics_seed(_MODEL.Dyn_opr);

            _LANCZOS_Dynamics.Perform_LANCZOS(_MODEL.Hamil);

        }

        if(Do_Dynamics && Dynamics_SPDOS){

            Mat_1_trio_int TRIO_VEC; Mat_1_doub values_;
            reading_input_dos_trio(inp_filename, TRIO_VEC, values_ );

            if(Below_mu){
                //----------BELOW CHEMICAL POTENTIAL-------
                BASIS_3_orb_Hubb_chain_GC _BASIS_Nm1;
                MODEL_3_orb_Hubb_chain_GC<BASIS_3_orb_Hubb_chain_GC> _MODEL_Nm1(_BASIS_Nm1);


                _MODEL_Nm1.Read_parameters(inp_filename);

                _BASIS_Nm1.N_total = _BASIS_Nm1.N_total - 1;
                _BASIS_Nm1.Construct_basis();

                cout<<"Size of Hilbert space = "<<_BASIS_Nm1.D_up_basis.size()<<endl;
                cout<<scientific<<setprecision(6);

                cout<<"Diagonal part started"<<endl;
                _MODEL_Nm1.Add_diagonal_terms();
                cout<<"Diagonal part done"<<endl;

                cout<<"Non Diagonal part started"<<endl;
                _MODEL_Nm1.Add_non_diagonal_terms();
                cout<<"Non Diagonal part done"<<endl;


                _MODEL_Nm1.Add_Spin_Orbit_Coupling();

                cout<<"Connections started"<<endl;
                _MODEL_Nm1.Add_connections();
                cout<<"Connections done"<<endl;


                cout<<"Size of Hilbert space = "<<_MODEL_Nm1.Hamil.nrows<<endl;
                cout<<scientific<<setprecision(1);
                // Print_Matrix_COO(_MODEL.Hamil);
                cout<<scientific<<setprecision(6);

                _MODEL_Nm1.Read_parameters_for_dynamics(inp_filename);


                LANCZOS _LANCZOS_Dynamics_DOS;
                _LANCZOS_Dynamics_DOS.Dynamics_performed=true;
                _LANCZOS_Dynamics_DOS.Read_Lanczos_parameters(inp_filename);
                _LANCZOS_Dynamics_DOS.Eig_vec=_LANCZOS.Eig_vec;
                _LANCZOS_Dynamics_DOS.GS_energy=_LANCZOS.GS_energy;

                _MODEL.Get_c_on_GS(_LANCZOS_Dynamics_DOS, _BASIS_Nm1,TRIO_VEC, values_ );
                _LANCZOS_Dynamics_DOS.Get_Dynamics_seed(_MODEL.State_c_on_GS);

                _LANCZOS_Dynamics_DOS.omega_sign=-1.0;
                _LANCZOS_Dynamics_DOS.file_dynamics_out = _LANCZOS_Dynamics_DOS.file_dynamics_out + "_below_mu.txt";
                _LANCZOS_Dynamics_DOS.Perform_LANCZOS(_MODEL_Nm1.Hamil);

                //-----------------------
            }

            if(Above_mu){
                //----------ABOVE CHEMICAL POTENTIAL-------
                BASIS_3_orb_Hubb_chain_GC _BASIS_Np1;
                MODEL_3_orb_Hubb_chain_GC<BASIS_3_orb_Hubb_chain_GC> _MODEL_Np1(_BASIS_Np1);

                _MODEL_Np1.Read_parameters(inp_filename);

                _BASIS_Np1.N_total = _BASIS_Np1.N_total + 1;
                _BASIS_Np1.Construct_basis();

                cout<<"Size of Hilbert space = "<<_BASIS_Np1.D_up_basis.size()<<endl;
                cout<<scientific<<setprecision(6);

                cout<<"Diagonal part started"<<endl;
                _MODEL_Np1.Add_diagonal_terms();
                cout<<"Diagonal part done"<<endl;

                cout<<"Non Diagonal part started"<<endl;
                _MODEL_Np1.Add_non_diagonal_terms();
                cout<<"Non Diagonal part done"<<endl;


                _MODEL_Np1.Add_Spin_Orbit_Coupling();

                cout<<"Connections started"<<endl;
                _MODEL_Np1.Add_connections();
                cout<<"Connections done"<<endl;


                cout<<"Size of Hilbert space = "<<_MODEL_Np1.Hamil.nrows<<endl;
                cout<<scientific<<setprecision(1);
                // Print_Matrix_COO(_MODEL.Hamil);
                cout<<scientific<<setprecision(6);

                _MODEL_Np1.Read_parameters_for_dynamics(inp_filename);

                LANCZOS _LANCZOS_Dynamics_DOS2;
                _LANCZOS_Dynamics_DOS2.Dynamics_performed=true;
                _LANCZOS_Dynamics_DOS2.Read_Lanczos_parameters(inp_filename);
                _LANCZOS_Dynamics_DOS2.Eig_vec=_LANCZOS.Eig_vec;
                _LANCZOS_Dynamics_DOS2.GS_energy=_LANCZOS.GS_energy;

                _MODEL.Get_cdagger_on_GS(_LANCZOS_Dynamics_DOS2, _BASIS_Np1,TRIO_VEC, values_ );
                _LANCZOS_Dynamics_DOS2.Get_Dynamics_seed(_MODEL.State_cdagger_on_GS);

                _LANCZOS_Dynamics_DOS2.omega_sign=1.0;
                _LANCZOS_Dynamics_DOS2.file_dynamics_out = _LANCZOS_Dynamics_DOS2.file_dynamics_out + "_above_mu.txt";
                _LANCZOS_Dynamics_DOS2.Perform_LANCZOS(_MODEL_Np1.Hamil);
                //----------------------------------------
            }

        }



    }



    //=======================================================================================================================================================================================================//
    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //=======================================================================================================================================================================================================//




    if (model_name=="3_orb_Hubbard_chain_GC" && Restricted_Basis==true) {

        BASIS_3_orb_Hubb_chain_GC_Restricted _BASIS;
        MODEL_3_orb_Hubb_chain_GC<BASIS_3_orb_Hubb_chain_GC_Restricted> _MODEL(_BASIS);


        //---------------------------------------

        _MODEL.Read_parameters(inp_filename);
        _BASIS.Construct_basis();


        cout<<"Size of Hilbert space = "<<_BASIS.D_up_basis.size()<<endl;
        cout<<scientific<<setprecision(6);

        cout<<"Diagonal part started"<<endl;
        _MODEL.Add_diagonal_terms();
        cout<<"Diagonal part done"<<endl;

        cout<<"Non Diagonal part started"<<endl;
        _MODEL.Add_non_diagonal_terms();
        cout<<"Non Diagonal part done"<<endl;


        _MODEL.Add_Spin_Orbit_Coupling();

        cout<<"Connections started"<<endl;
        _MODEL.Add_connections();
        cout<<"Connections done"<<endl;




        cout<<"Size of Hilbert space = "<<_MODEL.Hamil.nrows<<endl;
        cout<<scientific<<setprecision(1);
        // Print_Matrix_COO(_MODEL.Hamil);
        cout<<scientific<<setprecision(6);


        LANCZOS _LANCZOS;
        _LANCZOS.Dynamics_performed=false;
        _LANCZOS.Read_Lanczos_parameters(inp_filename);
        _LANCZOS.Perform_LANCZOS(_MODEL.Hamil);
        _LANCZOS.Write_full_spectrum();
        Print_vector_in_file(_LANCZOS.Eig_vec,"seed_GS.txt");

        _MODEL.Initialize_one_point_to_calculate();
        _MODEL.Initialize_two_point_to_calculate();

        for(int i=0;i<_LANCZOS.states_to_look.size();i++){
            cout<<"===================FOR STATE NO "<<i<<"============================="<<endl;
            _LANCZOS.Measure_one_point_observables(_MODEL.one_point_obs, _MODEL.One_point_oprts, _BASIS.Length, i);
            cout<<"Energy = "<<_LANCZOS.Evals_Tri_all[_LANCZOS.Evals_Tri_all.size()-1][i]<<endl;
            _LANCZOS.Measure_two_point_observables(_MODEL.two_point_obs, _MODEL.Two_point_oprts, _BASIS.Length, i, _MODEL.PBC);
            _LANCZOS.Measure_two_point_observables_smartly(_MODEL.one_point_obs,_MODEL.One_point_oprts, _BASIS.Length, i, "3_orb_Hubbard_chain_GC");

            cout<<"============================================================================"<<endl;
        }


        // for(int i=0;i<19;i++){
        // _LANCZOS.Measure_two_point_observables(_MODEL.two_point_obs, _MODEL.Two_point_oprts, _BASIS.Length, 0, _MODEL.PBC);

        //_LANCZOS.Measure_two_point_observables_smartly(_MODEL.one_point_obs,_MODEL.One_point_oprts, _BASIS.Length, 0);
        // }


        //_MODEL.Calculate_Local_Obs_for_States_to_Look(_LANCZOS,_BASIS);





        if(Do_Dynamics){

            _MODEL.Read_parameters_for_dynamics(inp_filename);

            _MODEL.Initialize_Opr_for_Dynamics(_LANCZOS);

            LANCZOS _LANCZOS_Dynamics;
            _LANCZOS_Dynamics.Dynamics_performed=true;
            _LANCZOS_Dynamics.Read_Lanczos_parameters(inp_filename);
            _LANCZOS_Dynamics.Eig_vec=_LANCZOS.Eig_vec;
            _LANCZOS_Dynamics.GS_energy=_LANCZOS.GS_energy;
            _LANCZOS_Dynamics.Get_Dynamics_seed(_MODEL.Dyn_opr);

            _LANCZOS_Dynamics.Perform_LANCZOS(_MODEL.Hamil);

        }



        Mat_1_real Eigen_ED;
        Mat_1_doub vecG;
        if(_MODEL.Hamil.nrows<700){
            DO_FULL_DIAGONALIZATION=true;
        }
        if(DO_FULL_DIAGONALIZATION==true){
            cout<<"//*****************Exact Diagonalization Energies**************//"<<endl;
            cout<<"//------------------------------------------------------------//"<<endl;
            Diagonalize(_MODEL.Hamil, Eigen_ED, vecG);
            for(int i=0;i<Eigen_ED.size();i++){
                cout<<i<<"   "<<Eigen_ED[i]<<endl;
            }
        }



    }



    //=======================================================================================================================================================================================================//
    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //=======================================================================================================================================================================================================//





#ifndef USE_COMPLEX
    if (model_name=="1_orb_Hubbard_chain") {

        Mat_1_doub Null_double_vec;
        Null_double_vec.resize(1);
        Mat_1_int Null_int_vec;
        Null_int_vec.resize(1);


        bool Dynamics_SPDOS = true;
        bool Above_mu = true;
        bool Below_mu= true;
        bool Cheaper_observables=true;

        MODEL_1_orb_Hubb_chain _MODEL;
        BASIS_1_orb_Hubb_chain _BASIS;

        //---------------------------------------

        _MODEL.Read_parameters(_BASIS,inp_filename);
        _BASIS.Construct_basis();
        cout<<"Basis contructed"<<endl;
        _MODEL.Add_diagonal_terms(_BASIS);
        cout<<"Diagonal terms added"<<endl;
        _MODEL.Add_connections(_BASIS);
        cout<<"Connections added"<<endl;


        cout<<"Size of Hilbert space = "<<_MODEL.Hamil.nrows<<endl;
        cout<<"Sparsity = "<<(1.0*_MODEL.Hamil.value.size())/(1.0*_MODEL.Hamil.nrows*_MODEL.Hamil.nrows)<<endl;
        cout<<scientific<<setprecision(20);
        //  Print_Matrix_COO(_MODEL.Hamil);

        LANCZOS _LANCZOS;
        _LANCZOS.Dynamics_performed=false;

        double EG;
        Mat_1_doub vecG;

        if(_MODEL.Hamil.nrows<400){
            DO_FULL_DIAGONALIZATION=true;
        }
        if(DO_FULL_DIAGONALIZATION==true){
            double EG;
            Mat_1_real Evals_temp;
            Mat_1_doub vecG;
            Diagonalize(_MODEL.Hamil, Evals_temp, vecG);
            cout<<"GS energy from ED(without Lanczos) = "<<Evals_temp[0]<<endl;
            cout<<"All eigenvalues using ED------------------------------"<<endl;
            cout<<"-------------------------------------------------------"<<endl;
            for(int i=0;i<Evals_temp.size();i++){
                cout<<i<<"  "<<Evals_temp[i]<<endl;
            }
            cout<<"-------------------------------------------------------"<<endl;
            cout<<"-------------------------------------------------------"<<endl;
        }


        _LANCZOS.Read_Lanczos_parameters(inp_filename);
        _LANCZOS.Perform_LANCZOS(_MODEL.Hamil);
        _LANCZOS.Write_full_spectrum();


        if(Cheaper_observables){

            for(int state_no=0;state_no<_LANCZOS.states_to_look.size();state_no++){

                cout <<"**********FOR STATE NO = "<<state_no<<"*********************"<<endl;
                cout<<"***************************************************"<<endl;

                Matrix_COO OPR_;
                Mat_1_string opr_type_;
                /*
            Matrix_COO OPR_;
            double sum_;

            Mat_1_string opr_type_;
            opr_type_.push_back("SzSz");
            opr_type_.push_back("SpSm");
            opr_type_.push_back("SmSp");

            Mat_2_doub Corr_;
            Corr_.resize(_BASIS.Length);
            for(int site1=0;site1<_BASIS.Length;site1++){
                Corr_[site1].resize(_BASIS.Length);
            }

            for(int type=0;type<3;type++){
             sum_=0.0;
                cout<<opr_type_[type]<<": "<<endl;

             for(int site1=0;site1<_BASIS.Length;site1++){
                 for(int site2=site1;site2<_BASIS.Length;site2++){
                     OPR_.columns.clear();
                     OPR_.rows.clear();
                     OPR_.value.clear();
                     _MODEL.Initialize_two_point_operator_sites_specific(opr_type_[type] , OPR_, site1, site2, _BASIS);
                     Corr_[site1][site2]=_LANCZOS.Measure_observable(OPR_, 0);
                     if(site1 != site2){
                    Corr_[site2][site1]=Corr_[site1][site2];
                     }
                     OPR_.columns.clear();
                     OPR_.rows.clear();
                     OPR_.value.clear();
                 }
             }
             for(int site1=0;site1<_BASIS.Length;site1++){
                 for(int site2=0;site2<_BASIS.Length;site2++){
                cout<< Corr_[site1][site2]<<" ";
                sum_ +=Corr_[site1][site2];
                 }
                cout<<endl;
             }
            cout<<"sum = "<<sum_<<endl;
            }
             */

                //one_point_observables
                double value_one_point;
                cout<<"One point observables:"<<endl;
                opr_type_.clear();
                opr_type_.push_back("n_up");
                opr_type_.push_back("n_dn");
                for(int type=0;type<opr_type_.size();type++){
                    for(int site=0;site<_BASIS.Length;site++){
                        OPR_.columns.clear();
                        OPR_.rows.clear();
                        OPR_.value.clear();
                        _MODEL.Initialize_one_point_operator_site_specific(opr_type_[type] , OPR_, site, _BASIS);
                        value_one_point=_LANCZOS.Measure_observable(OPR_, state_no);
                        cout<<opr_type_[type]<<"["<<site<<"]="<<value_one_point<<endl;
                    }
                }

                //two_point_observables
                double value_two_point;
                cout<<"Two point observables:"<<endl;
                opr_type_.clear();
                opr_type_.push_back("SzSz");
                opr_type_.push_back("SpSm");
                opr_type_.push_back("SmSp");
                double Sum_;
                for(int type=0;type<opr_type_.size();type++){
                    cout<<"________________________________"<<endl;
                    cout<<opr_type_[type]<<":"<<endl;
                    cout<<"````````````````````````````````"<<endl;
                    Sum_=0.0;
                    for(int site1=0;site1<_BASIS.Length;site1++){
                        for(int site2=0;site2<_BASIS.Length;site2++){
                            OPR_.columns.clear();
                            OPR_.rows.clear();
                            OPR_.value.clear();
                            _MODEL.Initialize_two_point_operator_sites_specific(opr_type_[type] , OPR_, site1, site2, _BASIS);
                            value_two_point=_LANCZOS.Measure_observable(OPR_, state_no);
                            cout<<value_two_point<<"  ";
                            Sum_+=value_two_point;
                        }
                        cout<<endl;
                    }
                    cout<<"Sum="<<Sum_<<endl;

                }

                vector< int >().swap( OPR_.columns );
                vector< int >().swap( OPR_.rows );
                vector< double_type >().swap( OPR_.value );
            }
        }
        else{
            _MODEL.Initialize_one_point_to_calculate(_BASIS);
            _LANCZOS.Measure_one_point_observables(_MODEL.one_point_obs, _MODEL.One_point_oprts, _BASIS.Length, 0);

            _MODEL.Initialize_two_point_to_calculate(_BASIS);
            _LANCZOS.Measure_two_point_observables(_MODEL.two_point_obs, _MODEL.Two_point_oprts, _BASIS.Length, 0, _MODEL.PBC);

            _LANCZOS.Measure_four_point_observables(_MODEL.Two_point_oprts, _MODEL.Four_point_sites_set, 0);
        }


        if(Do_Dynamics && !Dynamics_SPDOS){

            _MODEL.Read_parameters_for_dynamics(inp_filename);

            _MODEL.Initialize_Opr_for_Dynamics(_BASIS);


            LANCZOS _LANCZOS_Dynamics;

            _LANCZOS_Dynamics.Dynamics_performed=true;
            _LANCZOS_Dynamics.Read_Lanczos_parameters(inp_filename);
            _LANCZOS_Dynamics.Eig_vec=_LANCZOS.Eig_vec;
            _LANCZOS_Dynamics.GS_energy=_LANCZOS.GS_energy;
            _LANCZOS_Dynamics.Get_Dynamics_seed(_MODEL.Dyn_opr);

            _LANCZOS_Dynamics.Perform_LANCZOS(_MODEL.Hamil);

        }
        if(Do_Dynamics && Dynamics_SPDOS){


            vector< int >().swap( _MODEL.Hamil.columns );
            vector< int >().swap( _MODEL.Hamil.rows );
            vector< double_type >().swap( _MODEL.Hamil.value );

            Mat_1_trio_int TRIO_VEC; Mat_1_doub values_;
            reading_input_dos_trio(inp_filename, TRIO_VEC, values_ );

            for(int s_=0;s_<TRIO_VEC.size();s_++){
                if(TRIO_VEC[0].spin_ != TRIO_VEC[s_].spin_){
                    cout<<"ERROR:: spins of all fermionic operators must be same for dynamics"<<endl;
                    assert(TRIO_VEC[0].spin_ == TRIO_VEC[s_].spin_);
                }
            }

            if(Below_mu){
                //----------BELOW CHEMICAL POTENTIAL-------
                BASIS_1_orb_Hubb_chain _BASIS_Nm1;
                MODEL_1_orb_Hubb_chain _MODEL_Nm1;

                _MODEL_Nm1.Read_parameters(_BASIS_Nm1,inp_filename);
                if(TRIO_VEC[0].spin_==0){
                    _BASIS_Nm1.Nup = _BASIS.Nup-1;
                    _BASIS_Nm1.Ndn = _BASIS.Ndn;
                }
                else{
                    assert(TRIO_VEC[0].spin_==1);
                    _BASIS_Nm1.Nup = _BASIS.Nup;
                    _BASIS_Nm1.Ndn = _BASIS.Ndn-1;
                }
                _BASIS_Nm1.Construct_basis();
                cout<<"Basis contructed"<<endl;

                cout<<"Diagonal part started"<<endl;
                _MODEL_Nm1.Add_diagonal_terms(_BASIS_Nm1);
                cout<<"Diagonal part done"<<endl;

                cout<<"Connections started"<<endl;
                _MODEL_Nm1.Add_connections(_BASIS_Nm1);
                cout<<"Connections done"<<endl;


                cout<<"Size of Hilbert space = "<<_MODEL_Nm1.Hamil.nrows<<endl;
                cout<<scientific<<setprecision(1);
                // Print_Matrix_COO(_MODEL.Hamil);
                cout<<scientific<<setprecision(6);

                _MODEL_Nm1.Read_parameters_for_dynamics(inp_filename);


                LANCZOS _LANCZOS_Dynamics_DOS;
                _LANCZOS_Dynamics_DOS.Dynamics_performed=true;
                _LANCZOS_Dynamics_DOS.Read_Lanczos_parameters(inp_filename);
                _LANCZOS_Dynamics_DOS.GS_energy=_LANCZOS.GS_energy;

                _MODEL.Get_c_on_GS(_LANCZOS, _BASIS_Nm1, _BASIS, TRIO_VEC, values_);
                _LANCZOS_Dynamics_DOS.Get_Dynamics_seed(_MODEL.State_c_on_GS);


                _LANCZOS_Dynamics_DOS.omega_sign=-1.0;
                _LANCZOS_Dynamics_DOS.file_dynamics_out = _LANCZOS_Dynamics_DOS.file_dynamics_out + "_below_mu.txt";

                // assert(false);
                _LANCZOS_Dynamics_DOS.Perform_LANCZOS(_MODEL_Nm1.Hamil);


                vector< int >().swap( _MODEL_Nm1.Hamil.columns );
                vector< int >().swap( _MODEL_Nm1.Hamil.rows );
                vector< double_type >().swap( _MODEL_Nm1.Hamil.value );


                //-----------------------
            }

            if(Above_mu){
                //----------ABOVE CHEMICAL POTENTIAL-------
                BASIS_1_orb_Hubb_chain _BASIS_Np1;
                MODEL_1_orb_Hubb_chain _MODEL_Np1;

                _MODEL_Np1.Read_parameters(_BASIS_Np1,inp_filename);

                if(TRIO_VEC[0].spin_==0){
                    _BASIS_Np1.Nup = _BASIS.Nup+1;
                    _BASIS_Np1.Ndn = _BASIS.Ndn;
                }
                else{
                    assert(TRIO_VEC[0].spin_==1);
                    _BASIS_Np1.Nup = _BASIS.Nup;
                    _BASIS_Np1.Ndn = _BASIS.Ndn+1;
                }

                _BASIS_Np1.Construct_basis();

                cout<<"Size of Hilbert space = "<<_BASIS_Np1.D_up_basis.size()<<endl;
                cout<<scientific<<setprecision(6);

                cout<<"Diagonal part started"<<endl;
                _MODEL_Np1.Add_diagonal_terms(_BASIS_Np1);
                cout<<"Diagonal part done"<<endl;

                cout<<"Connections started"<<endl;
                _MODEL_Np1.Add_connections(_BASIS_Np1);
                cout<<"Connections done"<<endl;

                cout<<"Size of Hilbert space = "<<_MODEL_Np1.Hamil.nrows<<endl;
                cout<<scientific<<setprecision(1);
                // Print_Matrix_COO(_MODEL.Hamil);
                cout<<scientific<<setprecision(6);

                _MODEL_Np1.Read_parameters_for_dynamics(inp_filename);

                LANCZOS _LANCZOS_Dynamics_DOS2;
                _LANCZOS_Dynamics_DOS2.Dynamics_performed=true;
                _LANCZOS_Dynamics_DOS2.Read_Lanczos_parameters(inp_filename);
                _LANCZOS_Dynamics_DOS2.GS_energy=_LANCZOS.GS_energy;

                _MODEL.Get_cdagger_on_GS(_LANCZOS, _BASIS_Np1, _BASIS, TRIO_VEC, values_ );
                _LANCZOS_Dynamics_DOS2.Get_Dynamics_seed(_MODEL.State_cdagger_on_GS);

                _LANCZOS_Dynamics_DOS2.omega_sign=1.0;
                _LANCZOS_Dynamics_DOS2.file_dynamics_out = _LANCZOS_Dynamics_DOS2.file_dynamics_out + "_above_mu.txt";
                _LANCZOS_Dynamics_DOS2.Perform_LANCZOS(_MODEL_Np1.Hamil);
                //----------------------------------------

                vector< int >().swap( _MODEL_Np1.Hamil.columns );
                vector< int >().swap( _MODEL_Np1.Hamil.rows );
                vector< double_type >().swap( _MODEL_Np1.Hamil.value );
            }

        }

    }
#endif


    //=======================================================================================================================================================================================================//
    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //=======================================================================================================================================================================================================//



#ifndef USE_COMPLEX
    if (model_name=="1_orb_tJ") {


        MODEL_1_orb_tJ _MODEL;
        BASIS_1_orb_tJ _BASIS;

        //---------------------------------------

        _MODEL.Read_parameters(_BASIS,inp_filename);
        _BASIS.Construct_basis();
        cout<<"Basis contructed"<<endl;
        _MODEL.Add_diagonal_terms(_BASIS, "GS");
        cout<<"Diagonal terms added"<<endl;
        _MODEL.Add_connections(_BASIS, "GS");
        cout<<"Connections added"<<endl;



        cout<<"Size of Hilbert space = "<<_MODEL.Hamil.nrows<<endl;
        cout<<"Sparsity = "<<(1.0*_MODEL.Hamil.value.size())/(1.0*_MODEL.Hamil.nrows*_MODEL.Hamil.nrows)<<endl;
        cout<<scientific<<setprecision(20);
        //Print_Matrix_COO(_MODEL.Hamil);

        LANCZOS _LANCZOS;
        _LANCZOS.Dynamics_performed=false;

        double EG;
        Mat_1_doub vecG;
        if(DO_FULL_DIAGONALIZATION==true){
            Diagonalize(_MODEL.Hamil, EG, vecG);
            cout<<"GS enegry from ED(without Lanczos) = "<<EG<<endl;
        }

        _LANCZOS.Read_Lanczos_parameters(inp_filename);
        _LANCZOS.Perform_LANCZOS(_MODEL.Hamil);
        _LANCZOS.Write_full_spectrum();


        _MODEL.Initialize_one_point_to_calculate_from_file(_BASIS);
        _LANCZOS.Measure_one_point_observables(_MODEL.one_point_obs, _MODEL.One_point_oprts, _BASIS.Length, 0);

        Matrix_COO OPR_;
        Mat_1_string opr_type_;
        //two_point_observables
        double value_two_point;
        cout<<"Two point observables:"<<endl;
        opr_type_.clear();
        opr_type_.push_back("SzSz");
        opr_type_.push_back("SpSm");
        opr_type_.push_back("SmSp");
        double Sum_;
        for(int type=0;type<opr_type_.size();type++){
            cout<<"________________________________"<<endl;
            cout<<opr_type_[type]<<":"<<endl;
            cout<<"````````````````````````````````"<<endl;
            Sum_=0.0;
            for(int site1=0;site1<_BASIS.Length;site1++){
                for(int site2=0;site2<_BASIS.Length;site2++){
                    OPR_.columns.clear();
                    OPR_.rows.clear();
                    OPR_.value.clear();
                    _MODEL.Initialize_two_point_operator_sites_specific(opr_type_[type] , OPR_, site1, site2, _BASIS);
                    value_two_point=_LANCZOS.Measure_observable(OPR_, 0);
                    cout<<value_two_point<<"  ";
                    Sum_+=value_two_point;

                    vector< int >().swap( OPR_.columns );
                    vector< int >().swap( OPR_.rows );
                    vector< double_type >().swap( OPR_.value );
                }
                cout<<endl;
            }
            cout<<"Sum="<<Sum_<<endl;

        }


        //cout<<"//------------------------------------------------TIME EVOLUTION STARTED------------------------------------------------//"<<endl;
        /*
        double A0 = 20.0;
        double tau_p = 2.0;
        double sigma =0.5;
        double omega = 2.0;

        double T_min=0;
        double T_max=8;
        double dt=0.0001;
        double tau;
        int T_slices = int(((T_max - T_min)/dt) + 0.5);

        double Energy_;

        _LANCZOS.Evolving_State = _LANCZOS.Eig_vecs[0];
        for(int slice_no=1;slice_no<T_slices;slice_no++){
            tau=T_min + dt*slice_no;

            _MODEL.Add_diagonal_terms(_BASIS, "time_evolution");
            //cout<<"Diagonal terms added for time slice = "<<slice_no<<endl;

            _MODEL.Ax = A0*exp(-1.0*(((tau - tau_p)*(tau - tau_p) )/(2.0*sigma*sigma)))*cos(omega*(tau - tau_p));
            _MODEL.Ay = _MODEL.Ax;

            _MODEL.Add_connections(_BASIS, "time_evolution");
            //cout<<"Connections added  for time slice = "<<slice_no<<endl;

            _LANCZOS.Time_evolution_type1(_MODEL.Hamil, dt, Energy_);

            cout << tau<< "  "<<Energy_<<"   "<<_MODEL.Ax<<"  #Energy"<<endl;

            _LANCZOS.Measure_one_point_observables(_MODEL.one_point_obs, _MODEL.One_point_oprts, _BASIS.Length, -121);
            _LANCZOS.Measure_KE(_MODEL.H_KE, -121);
            _LANCZOS.Measure_Total_Energy(_MODEL.H_Total, -121);
            _LANCZOS.Measure_two_point_observables_smartly(_MODEL.one_point_obs,_MODEL.One_point_oprts, _BASIS.Length, -121, "1_orb_tJ");

        }

        */


        cout<<"Dynamics startedXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;

        if(Do_Dynamics){
            bool SpSm_Dyn=true;

            if(SpSm_Dyn){
                int site_=0;

                MODEL_1_orb_tJ _MODEL_Szm1;
                BASIS_1_orb_tJ _BASIS_Szm1;


                _MODEL_Szm1.Read_parameters(_BASIS_Szm1,inp_filename);
                _BASIS_Szm1.Nup = _BASIS.Nup - 1;
                _BASIS_Szm1.Ndn = _BASIS.Ndn + 1;
                _BASIS_Szm1.Construct_basis();
                _MODEL_Szm1.Add_diagonal_terms(_BASIS_Szm1, "GS");
                _MODEL_Szm1.Add_connections(_BASIS_Szm1, "GS");

                cout<<"Size of Hilbert space = "<<_MODEL_Szm1.Hamil.nrows<<endl;
                cout<<scientific<<setprecision(1);
                // Print_Matrix_COO(_MODEL.Hamil);
                cout<<scientific<<setprecision(6);

                _MODEL_Szm1.Read_parameters_for_dynamics(inp_filename);


                LANCZOS _LANCZOS_Dynamics_SpSm;
                _LANCZOS_Dynamics_SpSm.Dynamics_performed=true;
                _LANCZOS_Dynamics_SpSm.Read_Lanczos_parameters(inp_filename);
                _LANCZOS_Dynamics_SpSm.Eig_vec=_LANCZOS.Eig_vec;
                _LANCZOS_Dynamics_SpSm.GS_energy=_LANCZOS.GS_energy;

                _MODEL.Get_Sm_on_GS(_LANCZOS.Eig_vec, _BASIS_Szm1, _BASIS, site_);
                _LANCZOS_Dynamics_SpSm.Get_Dynamics_seed(_MODEL.State_Sm_on_GS);

                cout<<"size of seed = "<<_LANCZOS_Dynamics_SpSm.Dynamics_seed.size()<<endl;
                _LANCZOS_Dynamics_SpSm.omega_sign=1.0;
                _LANCZOS_Dynamics_SpSm.file_dynamics_out = _LANCZOS_Dynamics_SpSm.file_dynamics_out + "SpSm.txt";
                _LANCZOS_Dynamics_SpSm.Perform_LANCZOS(_MODEL_Szm1.Hamil);

                //-----------------------
            }

        }


    }

#endif

    //=======================================================================================================================================================================================================//
    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //=======================================================================================================================================================================================================//


    if (model_name=="1_orb_Hubbard_GC") {


        BASIS_1_orb_Hubbard_GC _BASIS;
        MODEL_1_orb_Hubbard_GC<BASIS_1_orb_Hubbard_GC> _MODEL(_BASIS);


        _MODEL.Read_parameters(inp_filename);
        _BASIS.Construct_basis();
        cout<<"Basis contructed"<<endl;
        _MODEL.Add_diagonal_terms();
        cout<<"Diagonal terms added"<<endl;
        _MODEL.Add_connections();
        cout<<"Connections added"<<endl;



        cout<<"Size of Hilbert space = "<<_MODEL.Hamil.nrows<<endl;
        cout<<"Sparsity = "<<(1.0*_MODEL.Hamil.value.size())/(1.0*_MODEL.Hamil.nrows*_MODEL.Hamil.nrows)<<endl;
        cout<<scientific<<setprecision(16);
        //Print_Matrix_COO(_MODEL.Hamil);

        LANCZOS _LANCZOS;
        _LANCZOS.Dynamics_performed=false;



        if(DO_FULL_DIAGONALIZATION==true && (_MODEL.Hamil.nrows<400)){

            Mat_1_real Evals_temp;
            Mat_1_doub vecG;
            Diagonalize(_MODEL.Hamil, Evals_temp, vecG);
            cout<<"GS energy from ED(without Lanczos) = "<<Evals_temp[0]<<endl;
            cout<<"All eigenvalues using EG------------------------------"<<endl;
            cout<<"-------------------------------------------------------"<<endl;
            for(int i=0;i<Evals_temp.size();i++){
                cout<<i<<"  "<<Evals_temp[i]<<endl;
            }
            cout<<"-------------------------------------------------------"<<endl;
            cout<<"-------------------------------------------------------"<<endl;
        }

        _LANCZOS.Read_Lanczos_parameters(inp_filename);
        _LANCZOS.Perform_LANCZOS(_MODEL.Hamil);
        _LANCZOS.Write_full_spectrum();

        cout<<scientific<<setprecision(6);
        _MODEL.Calculate_one_point_observables(_LANCZOS.Eig_vec);
        _MODEL.Calculate_two_point_observables(_LANCZOS.Eig_vec);





        bool Dynamics_SPDOS = true;
        bool Above_mu = true;
        bool Below_mu= true;


        if(Do_Dynamics && Dynamics_SPDOS){


            vector< int >().swap( _MODEL.Hamil.columns );
            vector< int >().swap( _MODEL.Hamil.rows );
            vector< double_type >().swap( _MODEL.Hamil.value );

            Mat_1_trio_int TRIO_VEC;
            Mat_1_doub values_;
            reading_input_dos_trio(inp_filename, TRIO_VEC, values_ );


            if(Below_mu){
                //----------BELOW CHEMICAL POTENTIAL-------
                BASIS_1_orb_Hubbard_GC _BASIS_Nm1;
                MODEL_1_orb_Hubbard_GC<BASIS_1_orb_Hubbard_GC> _MODEL_Nm1(_BASIS_Nm1);


                _MODEL_Nm1.Read_parameters(inp_filename);

                _BASIS_Nm1.N_total = _BASIS.N_total-1;
                _BASIS_Nm1.Construct_basis();
                cout<<"Basis contructed"<<endl;

                _MODEL_Nm1.Add_diagonal_terms();
                cout<<"Diagonal terms added"<<endl;

                _MODEL_Nm1.Add_connections();
                cout<<"Connections added"<<endl;


                cout<<"Size of Hilbert space = "<<_MODEL_Nm1.Hamil.nrows<<endl;
                cout<<scientific<<setprecision(6);

                _MODEL_Nm1.Read_parameters_for_dynamics(inp_filename);


                LANCZOS _LANCZOS_Dynamics_DOS;
                _LANCZOS_Dynamics_DOS.Dynamics_performed=true;
                _LANCZOS_Dynamics_DOS.Read_Lanczos_parameters(inp_filename);
                _LANCZOS_Dynamics_DOS.GS_energy=_LANCZOS.GS_energy;

                _MODEL.Get_c_on_GS(_LANCZOS, _BASIS_Nm1, TRIO_VEC, values_);
                _LANCZOS_Dynamics_DOS.Get_Dynamics_seed(_MODEL.State_c_on_GS);


                _LANCZOS_Dynamics_DOS.omega_sign=-1.0;
                _LANCZOS_Dynamics_DOS.file_dynamics_out = _LANCZOS_Dynamics_DOS.file_dynamics_out + "_below_mu.txt";

                // assert(false);
                _LANCZOS_Dynamics_DOS.Perform_LANCZOS(_MODEL_Nm1.Hamil);


                vector< int >().swap( _MODEL_Nm1.Hamil.columns );
                vector< int >().swap( _MODEL_Nm1.Hamil.rows );
                vector< double_type >().swap( _MODEL_Nm1.Hamil.value );


                //-----------------------
            }

            if(Above_mu){
                //----------ABOVE CHEMICAL POTENTIAL-------
                BASIS_1_orb_Hubbard_GC _BASIS_Np1;
                MODEL_1_orb_Hubbard_GC<BASIS_1_orb_Hubbard_GC> _MODEL_Np1(_BASIS_Np1);

                _MODEL_Np1.Read_parameters(inp_filename);

                _BASIS_Np1.N_total = _BASIS.N_total+1;
                _BASIS_Np1.Construct_basis();

                cout<<"Size of Hilbert space = "<<_BASIS_Np1.D_up_basis.size()<<endl;
                cout<<scientific<<setprecision(6);

                cout<<"Diagonal part started"<<endl;
                _MODEL_Np1.Add_diagonal_terms();
                cout<<"Diagonal part done"<<endl;

                cout<<"Connections started"<<endl;
                _MODEL_Np1.Add_connections();
                cout<<"Connections done"<<endl;

                cout<<"Size of Hilbert space = "<<_MODEL_Np1.Hamil.nrows<<endl;
                cout<<scientific<<setprecision(1);
                // Print_Matrix_COO(_MODEL.Hamil);
                cout<<scientific<<setprecision(6);

                _MODEL_Np1.Read_parameters_for_dynamics(inp_filename);

                LANCZOS _LANCZOS_Dynamics_DOS2;
                _LANCZOS_Dynamics_DOS2.Dynamics_performed=true;
                _LANCZOS_Dynamics_DOS2.Read_Lanczos_parameters(inp_filename);
                _LANCZOS_Dynamics_DOS2.GS_energy=_LANCZOS.GS_energy;

                _MODEL.Get_cdagger_on_GS(_LANCZOS, _BASIS_Np1, TRIO_VEC, values_ );
                _LANCZOS_Dynamics_DOS2.Get_Dynamics_seed(_MODEL.State_cdagger_on_GS);

                _LANCZOS_Dynamics_DOS2.omega_sign=1.0;
                _LANCZOS_Dynamics_DOS2.file_dynamics_out = _LANCZOS_Dynamics_DOS2.file_dynamics_out + "_above_mu.txt";
                _LANCZOS_Dynamics_DOS2.Perform_LANCZOS(_MODEL_Np1.Hamil);
                //----------------------------------------

                vector< int >().swap( _MODEL_Np1.Hamil.columns );
                vector< int >().swap( _MODEL_Np1.Hamil.rows );
                vector< double_type >().swap( _MODEL_Np1.Hamil.value );
            }

        }

    }



    //=======================================================================================================================================================================================================//
    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //=======================================================================================================================================================================================================//


    return 0;
}



