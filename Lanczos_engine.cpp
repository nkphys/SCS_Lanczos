#include "Lanczos_engine.h"
#include <stdlib.h>
#define PI 3.14159265
using namespace std;
//#define USE_COMPLEX
//#ifdef USE_COMPLEX


template <typename Basis_type, typename Model_type>
void LANCZOS<Basis_type, Model_type>::Perform_LANCZOS(Matrix_COO &Hamil){


    // Print_Matrix_COO(Hamil);
    B2.clear();A.clear();
    red_eig_vec.clear();Norms.clear();
    red_eig_vecs.clear();

    int seed_lanczos=Random_seed_value;
    int lanc_iter=0;
    double eps=Lanc_Error;
    int max_lanczos_states=min(max_steps,Hamil.nrows);
    double diff_E;
    double temp2, E0, E0_old, temp1;
    double temp2_type_double;
    double_type temp1_type_double,  temp3_type_double, temp4_type_double;

    int Target_state=0;
    Mat_1_doub Kvector_n,Kvector_nm1,Kvector_np1; //[element]
    Krylov_space_vecs.clear();
    double_type val_temp;
    complex<double> Recursive_GF, Recursive_GF_old;
    double diff_GF_imag;
    double norm_dynamics;
    double_type joker_double;
    Check_Ghosts=true;

    if(get_overlap_with_basis==true){
        need_few_eig_vecs=true;
        Get_Full_Spectrum=true;
    }

    if(need_few_eig_vecs==true){
        few_=states_to_look.size();
    }
    else{
        few_=1;
        states_to_look.clear();
        states_to_look.push_back(0);
    }
    Mat_1_int dummy_states_to_look;

    srand(seed_lanczos);
    double tmpnrm_type_double,tmpnrm_type_double2;
    double tmpnrm,tmpnrm2;


    //--------- CREATING SEED-----------------------------------//
    if(Dynamics_performed==false){

        if(Read_the_seed==false){

            if(Get_SeedVec_from_another_routine){
                Kvector_n = Seed_used;

                tmpnrm_type_double=Norm(Kvector_n);
                tmpnrm=sqrt(tmpnrm_type_double);
                for(int j=0;j<Hamil.nrows;j++){
                    Kvector_n[j] = (Kvector_n[j]/(tmpnrm));
                }
            }
            else{
                for(int j=0;j<Hamil.nrows;j++){
#ifdef USE_COMPLEX
                    temp1_type_double.real((rand()%RAND_MAX));
                    temp1_type_double.imag((rand()%RAND_MAX));
#endif
#ifndef USE_COMPLEX
                    temp1_type_double = (rand()%RAND_MAX);
#endif
                    temp1_type_double=(temp1_type_double)*(1.0/(RAND_MAX*1.0));
                    Kvector_n.push_back(temp1_type_double);
                }

                tmpnrm_type_double=Norm(Kvector_n);
                tmpnrm=sqrt(tmpnrm_type_double);

                for(int j=0;j<Hamil.nrows;j++){
                    Kvector_n[j] = (Kvector_n[j]/(tmpnrm));
                }
            }


        }
        else{
            ifstream infile_seed(seed_file_in.c_str());
            cout<<"Seed for Lanczos is being read from "<<seed_file_in<<endl;
            for(int j=0;j<Hamil.nrows;j++){

#ifdef USE_COMPLEX
                infile_seed>>tmpnrm;joker_double.real(tmpnrm);infile_seed>>tmpnrm;joker_double.imag(tmpnrm);
#endif
#ifndef USE_COMPLEX
                infile_seed>>joker_double;
#endif
                Kvector_n.push_back(joker_double);
            }
        }
    }
    else{
        Kvector_n=Dynamics_seed;
        tmpnrm_type_double=Norm(Kvector_n);
        tmpnrm=sqrt(tmpnrm_type_double);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for(int j=0;j<Hamil.nrows;j++){
            Kvector_n[j] = (Kvector_n[j]/(tmpnrm));
        }
        norm_dynamics=tmpnrm_type_double;
    }


    if(Save_the_Seed){
        Saved_Seed = Kvector_n;
    }
    //-------------SEED CREATED----------------------------//


    cout<<"In Lanczos milestone 1"<<endl;

    E0_old=0;
    diff_E=1.0;
    diff_GF_imag=1.0;
    Recursive_GF_old.imag(0);Recursive_GF_old.real(0);

    Evals_Tri_all.clear();

    bool do_next=true;
    while( do_next ){

        bool to_reorth;
        if(lanczos_reorthogonalization==true){
            if(lanc_iter>0 && lanc_iter%1==0){
                to_reorth=true;
            }
            else{
                to_reorth=false;
            }
        }
        else{
            to_reorth=false;
        }

        if(save_all_Krylov_space_vecs==true){
            Krylov_space_vecs.push_back(Kvector_n);
            //cout<<"Overlap b/w 0th and last Krylov vector : "<<dot_product(Krylov_space_vecs[0],Krylov_space_vecs[lanc_iter])<<endl;
        }

        Evals_Tri_all.resize(lanc_iter+1);
        clock_t Lanc_time = clock();

        temp1 = sqrt(Norm(Kvector_n));//*1.0e-10;
        Norms.push_back(temp1);
        if(lanc_iter==0){B2.push_back(0);}
        else{
            B2.push_back(tmpnrm*tmpnrm);
        }

        clock_t oprt_SB_time = clock();

        if(Saving_Hamil){
            Matrix_COO_vector_multiplication("U", Hamil, Kvector_n, Kvector_np1); // saved in K_vector_np1
        }
        else{
            model.Act_Hamil(basis, Kvector_n, Kvector_np1);
        }

        cout<<"Time to operate Hamiltonian : "<<double( clock() - oprt_SB_time ) / (double)CLOCKS_PER_SEC<<endl;//cout<<"here"<<endl;

        temp4_type_double = dot_product(Kvector_n, Kvector_np1);

        A.push_back(temp4_type_double);

#ifdef USE_COMPLEX
        assert(A[lanc_iter].imag()<0.0000001);
#endif

        Subtract(Kvector_np1, A[lanc_iter], Kvector_n);
        if(lanc_iter!=0){
            Subtract(Kvector_np1, sqrt(B2[lanc_iter]), Kvector_nm1);
        }

        if(save_all_Krylov_space_vecs==true){
            cout<<"Overlap b/w 0th and Kvector_np1 : "<<dot_product(Krylov_space_vecs[0],Kvector_np1)<<endl;
        }


        if(lanc_iter!=0){
            if(to_reorth==true){
                Mat_1_doub temp_dot;
                temp_dot.resize(Krylov_space_vecs.size());

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
                for(int i=0;i<Krylov_space_vecs.size();i++){
                    temp_dot[i]=dot_product(Kvector_np1,Krylov_space_vecs[i]);
                }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(val_temp)
#endif
                for(int l=0;l<Kvector_n.size();l++){
                    val_temp=zero;
                    for(int i=0;i<Krylov_space_vecs.size();i++){
                        val_temp = val_temp + Krylov_space_vecs[i][l]*temp_dot[i];
                    }
                    Kvector_np1[l] = Kvector_np1[l] - val_temp;
                }
                temp_dot.clear();

                cout<<"Reorthogonalization performed"<<endl;
            }
        }

        //Normalizaton of Knp1, added by myself, not included in std. Lanczos
        tmpnrm_type_double=Norm(Kvector_np1);
        tmpnrm=sqrt(tmpnrm_type_double);

        if((tmpnrm)>0.0000000000000000001){
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        for(int i=0;i<Kvector_np1.size();i++){
            Kvector_np1[i] = (Kvector_np1[i]/(tmpnrm));//*1.0e-10;
        }

    }
        if(save_all_Krylov_space_vecs==true){
            cout<<"Overlap b/w 0th and Normalized Kvector_np1 : "<<dot_product(Krylov_space_vecs[0],Kvector_np1)<<endl;
        }

        if(Get_Full_Spectrum==false){
            Diagonalize(A,B2,E0,red_eig_vec);
        }
        else{
            Evals_Tri_all[lanc_iter].resize(A.size());

            dummy_states_to_look.clear();
            for(int i_temp=0;i_temp<states_to_look.size();i_temp++){
                if(states_to_look[i_temp]<=lanc_iter){
                    dummy_states_to_look.push_back(states_to_look[i_temp]);
                }
            }

            Diagonalize(A,B2,Evals_Tri_all,red_eig_vecs,lanc_iter,dummy_states_to_look.size(),dummy_states_to_look);
            red_eig_vec=red_eig_vecs[0];
            E0=Evals_Tri_all[lanc_iter][0];
        }


        diff_E=	fabs(E0-E0_old);
        if(!Dynamics_performed){
            cout<<"iter = "<<lanc_iter<<" diff_E = "<<diff_E<<" E0 = "<<E0<<" E0_old = "<<E0_old<<endl;
            GS_energy=E0;
        }
        if(lanc_iter==0){if(Norm(Kvector_np1)==zero){diff_E = 0;}}
        //cout<<"Energy for lanc_iter("<<lanc_iter<<") is "<<E0<<endl;

        if(get_overlap_bw_eigstates_and_K_vecs ==true){
            cout<<"Overlap b/w low energy eigenvectors with Krylov vector-0:"<<endl;
            for(int no=0;no<red_eig_vecs.size();no++){
                cout<<"<K0|Psi"<<no<<"> = " << red_eig_vecs[no][0]<<endl;
            }
        }

        E0_old=E0;


        //Kvector_nm1.clear();
        Kvector_nm1=Kvector_n;
        //Kvector_n.clear();
        Kvector_n=Kvector_np1;

        if(lanc_iter<=Target_state){diff_E=1.0;}//doing altleast 2 iterations of Lanczos

        lanc_iter=lanc_iter+1;


        if(lanc_iter==Hamil.ncols){diff_E=0;}
        cout<<"Time for 1 LAnczos iter : "<<double( clock() - Lanc_time ) / (double)CLOCKS_PER_SEC<<endl;

        if(!Dynamics_performed){
            do_next=(diff_E>eps && lanc_iter<max_lanczos_states);}
        else{
            Calculate_recursive_GF(A,B2,Recursive_GF,check_omega,eta,GS_energy);
            diff_GF_imag=fabs(Recursive_GF.imag() - Recursive_GF_old.imag());
            Recursive_GF_old=Recursive_GF;
            do_next=(diff_GF_imag>eps && lanc_iter<max_lanczos_states);
            cout<<"iter = "<<lanc_iter<<" diff_GF_imag= "<<diff_GF_imag<<"  Recursive_GF.imag() = "<<Recursive_GF.imag()<<" Recursive_GF_old.imag() = "<<Recursive_GF_old.imag()<<endl;

        }

    }

    cout<<"Perform_LANCZOS: "<<"NO. of iterations required to get convergence in LANCZOS(pass 1) = "<<lanc_iter<<endl;
    cout<<"Perform_LANCZOS: "<<"Energy(GS)"<<" is "<<scientific<<setprecision(20)<< "Energy = "<<E0<<"   "<<"Lanczos_error = "<<diff_E<<endl;
    cout<<"Perform_LANCZOS: "<<"LANCZOS(pass 2) STARTING FOR Eigenvectors, "<<", Size of Matrix(SB) = "<<Hamil.ncols<<endl<<endl;

    Lanc_iter_done=lanc_iter;
    Kvector_n.clear();	Kvector_nm1.clear(); Kvector_np1.clear();




    if(Dynamics_performed==false){
        if(Eig_vecs_required){

            if(Get_Full_Spectrum == false){
                red_eig_vecs.clear();
                red_eig_vecs.resize(1);
                red_eig_vecs[0]=red_eig_vec;
            }

            int temp_Target_state;
            Eig_vecs.clear();
            Eig_vecs.resize(states_to_look.size());
            for(int Ts=0;Ts<states_to_look.size();Ts++){

                temp_Target_state=states_to_look[Ts];
                if(Ts==0){
                    assert(temp_Target_state==0);
                }
                srand(seed_lanczos);

                //------Creating Seed--------------------//
                if(Read_the_seed==false){

                    if(Get_SeedVec_from_another_routine){
                        Kvector_n = Seed_used;
                        tmpnrm_type_double=Norm(Kvector_n);
                        tmpnrm=sqrt(tmpnrm_type_double);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
                        for(int j=0;j<Hamil.nrows;j++){
                            Kvector_n[j] = (Kvector_n[j]/(tmpnrm));
                        }
                    }
                    else{
                        for(int j=0;j<Hamil.nrows;j++){
#ifdef USE_COMPLEX
                            temp1_type_double.real((rand()%RAND_MAX));
                            temp1_type_double.imag((rand()%RAND_MAX));
#endif
#ifndef USE_COMPLEX
                            temp1_type_double = (rand()%RAND_MAX);
#endif
                            temp1_type_double=(temp1_type_double)*(1.0/(RAND_MAX*1.0));
                            Kvector_n.push_back(temp1_type_double);
                        }

                        tmpnrm_type_double=Norm(Kvector_n);
                        tmpnrm=sqrt(tmpnrm_type_double);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
                        for(int j=0;j<Hamil.nrows;j++){
                            Kvector_n[j] = (Kvector_n[j]/(tmpnrm));
                        }

                    }
                }
                else{
                    ifstream infile_seed(seed_file_in.c_str());
                    for(int j=0;j<Hamil.nrows;j++){
#ifdef USE_COMPLEX
                        infile_seed>>tmpnrm;joker_double.real(tmpnrm);infile_seed>>tmpnrm;joker_double.imag(tmpnrm);
#endif
#ifndef USE_COMPLEX
                        infile_seed>>joker_double;
#endif
                        Kvector_n.push_back(joker_double);
                    }
                }
                //---------Seed Created----------------------//

                Eig_vecs[Ts].clear();
                for(int j=0;j<Hamil.nrows;j++){
                    Eig_vecs[Ts].push_back(0);
                }

                if(save_all_Krylov_space_vecs==true){
                    for(int lanc_iter2=0;lanc_iter2<lanc_iter;lanc_iter2=lanc_iter2+1){
                        Subtract(Eig_vecs[Ts], (-1.0*(red_eig_vecs[Ts][lanc_iter2])), Krylov_space_vecs[lanc_iter2]);
                    }

                }
                else{
                    for(int lanc_iter2=0;lanc_iter2<lanc_iter;lanc_iter2=lanc_iter2+1){
                        Subtract(Eig_vecs[Ts], (-1.0*(red_eig_vecs[Ts][lanc_iter2])), Kvector_n);

                        // saved in K_vector_np1
                        if(Saving_Hamil){
                            Matrix_COO_vector_multiplication("U", Hamil, Kvector_n, Kvector_np1); // saved in K_vector_np1
                        }
                        else{
                            model.Act_Hamil(basis, Kvector_n, Kvector_np1);
                        }

                        Subtract(Kvector_np1, A[lanc_iter2], Kvector_n);
                        if(lanc_iter2!=0){
                            Subtract(Kvector_np1, sqrt(B2[lanc_iter2]), Kvector_nm1);
                        }

                        //Normalizaton of Knp1 , not included in std. Lanczos
                        tmpnrm_type_double = Norm(Kvector_np1); //new
                        tmpnrm=sqrt(tmpnrm_type_double);
                        for(int i=0;i<Kvector_np1.size();i++){
                            Kvector_np1[i] = (Kvector_np1[i]/(tmpnrm));
                        }
                        Kvector_nm1=Kvector_n;
                        Kvector_n=Kvector_np1;
                    }
                }

                double norm_ev=Norm(Eig_vecs[Ts]);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
                for(int j=0;j<Hamil.nrows;j++){
                    Eig_vecs[Ts][j] = (Eig_vecs[Ts][j]/(sqrt(norm_ev)));
                }

                if(Get_insitu_FTLM_overlaps){
                    for(int ip=0;ip<Mat_elements.size();ip++){
                        Mat_elements[ip][Ts] = dot_product(Eig_vecs[Ts], Vecs_FTLM[ip]);
                    }
                    if(Ts>0){
                        vector < double_type >().swap(Eig_vecs[Ts]);
                    }
                }

                Kvector_n.clear();Kvector_nm1.clear();Kvector_np1.clear();
            }

            Eig_vec=Eig_vecs[0];

            if(get_overlap_with_basis==true){
                Overlaps_bare.clear();
                Overlaps_bare.resize(Eig_vecs.size());
                int bi_new;

                Overlaps.clear();
                Overlaps.resize(Eig_vecs.size());
                for(int Ts=0;Ts<Eig_vecs.size();Ts++){
                    Overlaps[Ts].resize(Eig_vecs[0].size());

                    for(int bi=0;bi<Eig_vecs[0].size();bi++){

                        Overlaps[Ts][bi].first=abs( Eig_vecs[Ts][bi])*abs( Eig_vecs[Ts][bi]) ;
                        Overlaps[Ts][bi].second=bi;
                    }
                    cout<<"sorting for state = "<<states_to_look[Ts]<<endl;
                    sort(Overlaps[Ts].begin(), Overlaps[Ts].end(), comp_greater_pair_double_int);
                }

                for(int Ts=0;Ts<Eig_vecs.size();Ts++){
                    Overlaps_bare[Ts].resize(Eig_vecs[0].size());
                    for(int bi=0;bi<Eig_vecs[0].size();bi++){
                        bi_new = Overlaps[Ts][bi].second;
                        Overlaps_bare[Ts][bi].first=Eig_vecs[Ts][bi_new];
                        Overlaps_bare[Ts][bi].second=bi_new;
                    }
                }

                ofstream outfile_overlap(overlap_out_file.c_str());
                for(int bi=0;bi<Eig_vecs[0].size();bi++){
                    outfile_overlap<<bi;
                    for(int Ts=0;Ts<Eig_vecs.size();Ts++){
                        outfile_overlap<<"  "<<Overlaps_bare[Ts][bi].first<<"   "<<Overlaps[Ts][bi].first;
                    }
                    for(int Ts=0;Ts<Eig_vecs.size();Ts++){
                        outfile_overlap<<"  "<<Overlaps_bare[Ts][bi].second;
                    }
                    outfile_overlap<<endl;
                }




                if(Check_Ghosts==true){
                    double_type value_local;
                    string out3 = "ghost.txt";
                    ofstream file_out3(out3.c_str());

                    for(int Ts1=0;Ts1<Eig_vecs.size();Ts1++){

                        for(int Ts2=0;Ts2<Eig_vecs.size();Ts2++){
                            value_local = dot_product(Eig_vecs[Ts1],Eig_vecs[Ts2]);
                            if(abs(value_local)<0.00000001){
                                file_out3<<zero<<"  ";
                            }
                            else{
                                file_out3<<value_local<<"  ";}

                        }
                        file_out3<<endl;

                    }
                }


            }


            /*
bool comp(int i, int j) { return i > j; }
sort(numbers.begin(), numbers.end(), comp);
 */

        }
    }

    else{
        string out = file_dynamics_out;
        ofstream file_out(out.c_str());
        double omega =omega_min;
        double sum_response=0.0;
        while(omega<=omega_max){
            omega = omega_sign*omega;
            Calculate_recursive_GF(A,B2,Recursive_GF,omega,eta,GS_energy);
            omega = omega_sign*omega;
            file_out<<scientific<<setprecision(10)<<omega<<"   "<<(-1.0/PI)*Numerator_Dynamics*(Recursive_GF.imag())<<"  "<<(-1.0/PI)*Numerator_Dynamics*(Recursive_GF.real())<<endl;
            sum_response +=(-1.0/PI)*Numerator_Dynamics*(Recursive_GF.imag());
            omega=omega+d_omega;
        }
        file_out<<"# sum response = "<<sum_response<<endl;
    }





    cout<<"----Lanczos milestone 2-----"<<endl;
    //cout<<TimeEvoPerformed<<endl;

    if(TimeEvoPerformed){
        Vec_new_TimeEvo.resize(Hamil.nrows);
    for(int i=0;i<Hamil.nrows;i++){
    Vec_new_TimeEvo[i]=0.0;
    }
    }

    for(int n=0;n<Krylov_space_vecs.size();n++){

        if(TimeEvoPerformed){
#ifdef USE_COMPLEX
        for(int i=0;i<Hamil.nrows;i++){ //sum over basis size
        for(int l=0;l<M_TimeEvo;l++){ //sum over eigenvectors

        Vec_new_TimeEvo[i] +=red_eig_vecs[l][n]*
                      Krylov_space_vecs[n][i]*
                     conjugate(red_eig_vecs[l][0])*
                     exp(-1.0*iota_comp* (Evals_Tri_all[Lanc_iter_done-1][l])*dt_TimeEvo);

        }
        }
#else
     cout<<"USE_COMPLEX FOR TIME EVOLUTION"<<endl;
     assert(false);
#endif
        }

        vector < double_type >().swap(Krylov_space_vecs[n]);
    }
    Krylov_space_vecs.clear();

    cout<<"----Lanczos milstone 3-----"<<endl;


}

template <typename Basis_type, typename Model_type>
void LANCZOS<Basis_type, Model_type>::Clear(){
    vector< double >().swap( B2 );
    vector< double >().swap( Norms );

    vector< double_type >().swap( A );
    vector< double_type >().swap( red_eig_vec );

    for(int i=0;i>red_eig_vecs.size();i++){
        vector< double_type >().swap( red_eig_vecs[i] );
    }
    red_eig_vecs.clear();


    for(int i=0;i<Krylov_space_vecs.size();i++){
        vector < double_type >().swap(Krylov_space_vecs[i]);
    }
    Krylov_space_vecs.clear();


    vector < double_type >().swap(Eig_vec);

    for(int i=0;i<Eig_vecs.size();i++){
        vector < double_type >().swap(Eig_vecs[i]);
    }
    Eig_vecs.clear();

    for(int i=0;i<Vecs_FTLM.size();i++){
        vector < double_type >().swap(Vecs_FTLM[i]);
    }
    Vecs_FTLM.clear();


}

template <typename Basis_type, typename Model_type>
void LANCZOS<Basis_type, Model_type>::Measure_macro_observables(Mat_1_string macro_obs, Hamiltonian_1_COO &Macro_oprts, int state_no){

    cout<<"Macro obs Measurement is started after completing Lanczos algorithm"<<endl;
    cout<<"Macro obs Measurement is started after completing Lanczos algorithm fo state_no = "<<states_to_look[state_no]<<endl;
    //Krylov_space_vecs



    double_type value;
    Mat_1_doub temp_vec;
    double_type sum_;


    sum_=zero;
    for(int opr_no=0;opr_no<macro_obs.size();opr_no++){
        Matrix_COO_vector_multiplication("U", Macro_oprts[opr_no], Eig_vecs[state_no], temp_vec);
        value = dot_product(Eig_vecs[state_no],temp_vec);
        cout<<macro_obs[opr_no]<<" =  "<<value<<endl;
        sum_ += value;
    }
    cout<<"TOTAL OF MACRO-OBSERVABLES : "<<sum_<<endl;
}



template <typename Basis_type, typename Model_type>
void LANCZOS<Basis_type, Model_type>::Measure_KE(Matrix_COO &Macro_oprt, int state_no){


    double_type value;
    Mat_1_doub temp_vec,vector_used;

    if(state_no==-121){
        vector_used=Evolving_State;
    }
    else{
        vector_used=Eig_vecs[state_no];
    }

    Matrix_COO_vector_multiplication("U", Macro_oprt, vector_used, temp_vec);
    value = dot_product(vector_used,temp_vec);
    cout<<value<<"#KE"<<endl;

}


template <typename Basis_type, typename Model_type>
void LANCZOS<Basis_type, Model_type>::Measure_Total_Energy(Matrix_COO &Macro_oprt, int state_no){


    double_type value;
    Mat_1_doub temp_vec,vector_used;

    if(state_no==-121){
        vector_used=Evolving_State;
    }
    else{
        vector_used=Eig_vecs[state_no];
    }

    Matrix_COO_vector_multiplication("U", Macro_oprt, vector_used, temp_vec);
    value = dot_product(vector_used,temp_vec);
    cout<<value<<"#TE"<<endl;

}




template <typename Basis_type, typename Model_type>
void LANCZOS<Basis_type, Model_type>::Time_evolution_type1(Matrix_COO &Hamil, double dt_, double &Energy_){

#ifdef USE_COMPLEX
    time_evolution_order=3;
    complex<double> iota_(0.0,1.0);
    Mat_1_doub temp_vec,temp_vec2,temp_vec3;

    if(time_evolution_order==1){

        Matrix_COO_vector_multiplication("U", Hamil, Evolving_State, temp_vec);
        Subtract(Evolving_State, iota_*dt_*(-1.0), temp_vec);

    }

    if(time_evolution_order==2){

        Matrix_COO_vector_multiplication("U", Hamil, Evolving_State, temp_vec);
        Matrix_COO_vector_multiplication("U", Hamil, temp_vec, temp_vec2);
        Subtract(Evolving_State, iota_*dt_*(-1.0), temp_vec);
        Subtract(Evolving_State, 0.5*dt_*dt_*(-1.0), temp_vec2);

    }

    if(time_evolution_order==3){

        Matrix_COO_vector_multiplication("U", Hamil, Evolving_State, temp_vec);
        Matrix_COO_vector_multiplication("U", Hamil, temp_vec, temp_vec2);
        Matrix_COO_vector_multiplication("U", Hamil, temp_vec2, temp_vec3);
        Subtract(Evolving_State, iota_*dt_*(-1.0), temp_vec);
        Subtract(Evolving_State, 0.5*dt_*dt_*(-1.0)*one, temp_vec2);
        Subtract(Evolving_State, iota_*(1/6.0)*dt_*dt_*dt_*(-1.0), temp_vec3);

    }

    //Normalizing it every time-----------------
    double norm_ev=dot_product(Evolving_State, Evolving_State).real();
    for(int j=0;j<Hamil.nrows;j++){
        Evolving_State[j] = (Evolving_State[j]/(sqrt(norm_ev)));
    }
    //---------------------------------

    Matrix_COO_vector_multiplication("U", Hamil, Evolving_State, temp_vec);
    Energy_=dot_product(Evolving_State,temp_vec).real();

    temp_vec.clear();
    temp_vec2.clear();
#endif

#ifndef USE_COMPLEX
    cout<<"TIME EVOLUTION WORKS ONLY IF USE_COMPLEX IS USED"<<endl;
    assert(false);
#endif
}

template <typename Basis_type, typename Model_type>
void LANCZOS<Basis_type, Model_type>::Measure_one_point_observables(Mat_1_string one_point_obs, Hamiltonian_2_COO &One_point_oprts, int T_sites, int state_no){

    cout<<"One-point Measurement is started after completing Lanczos algorithm"<<endl;

    double_type value;
    Mat_1_doub temp_vec;
    Mat_1_doub vector_used;

    if(state_no==-121){
        vector_used=Evolving_State;
    }
    else{
        vector_used=Eig_vecs[state_no];
    }

    One_point_observables_values.resize(One_point_oprts.size());
    for(int opr_no=0;opr_no<One_point_oprts.size();opr_no++){
        One_point_observables_values[opr_no].resize(T_sites);

        for(int site=0;site<T_sites;site++){
            Matrix_COO_vector_multiplication("cx", One_point_oprts[opr_no][site], vector_used, temp_vec);
            value = dot_product(temp_vec,vector_used);
            cout<<one_point_obs[opr_no]<<"["<<site<<"]  "<<value<<endl;
            One_point_observables_values[opr_no][site]=value;
        }
    }

    vector_used.clear();

    cout<<"one point obs done"<<endl;

}


template <typename Basis_type, typename Model_type>
void LANCZOS<Basis_type, Model_type>::Measure_two_point_observables_smartly(Mat_1_string one_point_obs, Hamiltonian_2_COO &One_point_oprts, int T_sites, int state_no, string _model_){
    cout<<"Two-point Measurement is done smartly for state no = "<<states_to_look[state_no]<<endl;

    double_type value;
    Mat_1_doub temp_vec, temp_vec2;
    Mat_3_doub corr;
    Matrix_COO TEMP_COO;
    Mat_1_doub vector_used;

    if(state_no==-121){
        vector_used=Evolving_State;
    }
    else{
        vector_used=Eig_vecs[state_no];
    }


    corr.resize(1);

    for(int i=0;i<1;i++){
        corr[i].resize(T_sites);
        for(int j=0;j<T_sites;j++){
            corr[i][j].resize(T_sites);
        }
    }

    for(int opr_no=0;opr_no<one_point_obs.size();opr_no++){
        cout<<"("<<one_point_obs[opr_no]<<"(i)|gs>)^\\dagger "<<one_point_obs[opr_no]<<"(j)|gs>"<<endl;
        for(int site=0;site<T_sites;site++){
            for(int site2=0;site2<T_sites;site2++){
                if(site2>=site){
                    value=zero;
                    //TEMP_COO = Dagger(Sz[site2]);
                    TEMP_COO = (One_point_oprts[opr_no][site]);
                    Matrix_COO_vector_multiplication("Full", One_point_oprts[opr_no][site2], vector_used, temp_vec);
                    Matrix_COO_vector_multiplication("Full", TEMP_COO, vector_used, temp_vec2);
                    value = dot_product(temp_vec2,temp_vec);
                    corr[0][site][site2]=value;
                    corr[0][site2][site]=value;
                    cout<<value<<"  ";
                }
                else{
                    cout<<zero<<"  ";
                }
            }
            cout<<endl;
        }

        double_type sum_all;
        sum_all=zero;
        for(int site=0;site<T_sites;site++){
            for(int site2=0;site2<T_sites;site2++){
                sum_all = sum_all + corr[0][site][site2];
            }
        }
        cout<<"sum_all = "<<sum_all<<endl;

    }


#ifdef USE_COMPLEX
    if(_model_=="3_orb_Hubbard_chain_GC"){

        if(false){
            cout<<"2point OBSERVABLES specifically for 3_orb_Hubbard_chain_GC"<<endl;

            //        Mat_1_int oprs_index;
            //        Mat_1_string oprs_string;
            //        oprs_string.push_back();


            for(int opr_no=14;opr_no<=32;opr_no++){
                cout<<"("<<one_point_obs[opr_no]<<"(i)|gs>)^\\dagger "<<one_point_obs[opr_no]<<"(j)|gs>"<<endl;
                for(int site=0;site<T_sites;site++){
                    for(int site2=0;site2<T_sites;site2++){
                        if(site2>=site){
                            value=zero;
                            //TEMP_COO = Dagger(Sz[site2]);
                            TEMP_COO = (One_point_oprts[opr_no][site]);
                            Matrix_COO_vector_multiplication("Full", One_point_oprts[opr_no][site2], Eig_vecs[state_no], temp_vec);
                            Matrix_COO_vector_multiplication("Full", TEMP_COO, Eig_vecs[state_no], temp_vec2);
                            value = dot_product(temp_vec2,temp_vec);
                            corr[0][site][site2]=value;
                            corr[0][site2][site]=value;
                            cout<<value<<"  ";
                        }
                        else{
                            cout<<zero<<"  ";
                        }
                    }
                    cout<<endl;
                }


                double_type sum_all=zero;
                for(int site=0;site<T_sites;site++){
                    for(int site2=0;site2<T_sites;site2++){
                        sum_all = sum_all + corr[0][site][site2];
                    }
                }
                cout<<"sum_all = "<<sum_all<<endl;
            }
        }


    }
#endif



}


template <typename Basis_type, typename Model_type>
double_type LANCZOS<Basis_type, Model_type>::Measure_observable(Matrix_COO &OPR_, int state_no){

    double_type value;
    string mult_type="not_U";
    Mat_1_doub temp_vec;

    Matrix_COO_vector_multiplication(mult_type, OPR_, Eig_vecs[state_no], temp_vec);
    value = dot_product(Eig_vecs[state_no],temp_vec);

    return value;

}


template <typename Basis_type, typename Model_type>
double_type LANCZOS<Basis_type, Model_type>::Measure_observable(Matrix_COO &OPR_, int state_no, string mult_type){

    double_type value;
    Mat_1_doub temp_vec;

    Matrix_COO_vector_multiplication(mult_type, OPR_, Eig_vecs[state_no], temp_vec);
    value = dot_product(Eig_vecs[state_no],temp_vec);

    return value;

}

template <typename Basis_type, typename Model_type>
void LANCZOS<Basis_type, Model_type>::Measure_four_point_observables(Hamiltonian_3_COO &Two_point_oprts, Mat_1_tetra_int sites_set, int state_no){

    cout<<"Four-point Measurement is started after completing Lanczos algorithm for state_no = "<<states_to_look[state_no]<<endl;
    cout<<"Only <GS|c_{iup}^{dagger}c_{jup}^{dagger}c_{lup}c_{kup}|GS> for 'j' ne 'l' is calculated"<<endl;

    Mat_1_doub temp_vec1,temp_vec2;
    double_type value;
    Matrix_COO Matrix_COO_temp;

    //<c_{first}^{\dagger}c_{second}^{\dagger}c_{third}c_{fourth}>
    // = (for second \ne third)
    //-1.0*<c_{first}^{\dagger}c_{third} c_{second}^{\dagger}c_{fourth}>

    for(int set=0;set<sites_set.size();set++){

        if(sites_set[set].second != sites_set[set].third){

            if(sites_set[set].fourth >= sites_set[set].second){
                Matrix_COO_vector_multiplication("x", Two_point_oprts[2][sites_set[set].second][sites_set[set].fourth],
                        Eig_vecs[state_no], temp_vec2);
            }
            else{
                Matrix_COO_temp=Dagger(Two_point_oprts[2][sites_set[set].fourth][sites_set[set].second]);
                Matrix_COO_vector_multiplication("x",Matrix_COO_temp ,
                                                 Eig_vecs[state_no], temp_vec2);

                Matrix_COO_temp.columns.clear();
                Matrix_COO_temp.rows.clear();
                Matrix_COO_temp.value.clear();

            }

            if(sites_set[set].first >= sites_set[set].third){
                Matrix_COO_vector_multiplication("x", Two_point_oprts[2][sites_set[set].third][sites_set[set].first],
                        Eig_vecs[state_no], temp_vec1);
            }
            else{
                Matrix_COO_temp=Dagger(Two_point_oprts[2][sites_set[set].first][sites_set[set].third]);
                Matrix_COO_vector_multiplication("x", Matrix_COO_temp,
                                                 Eig_vecs[state_no], temp_vec1);

                Matrix_COO_temp.columns.clear();
                Matrix_COO_temp.rows.clear();
                Matrix_COO_temp.value.clear();

            }

            value = dot_product(temp_vec2,temp_vec1);

            cout<< "<GS|cdup["<<sites_set[set].first<<"]"<<"cdup["<<sites_set[set].second<<"]"
                <<"cup["<<sites_set[set].third<<"]"<<"cup["<<sites_set[set].fourth<<"]|GS> = "
               << (value*(one*(-1.0)))<<endl;


        }

        else{
            cout<<"Do not use j==l"<<endl;
        }
    }

}

template <typename Basis_type, typename Model_type>
void LANCZOS<Basis_type, Model_type>::Measure_two_point_observables(Mat_1_string two_point_obs, Hamiltonian_3_COO &Two_point_oprts, int T_sites,  int state_no, bool PBC_check){

    cout<<"Two-point Measurement is started after completing Lanczos algorithm for state_no = "<<states_to_look[state_no]<<endl;

    double_type value,value1,value2;
    Mat_1_doub temp_vec;
    string mult_type;
    Mat_3_doub corr;
    corr.resize(two_point_obs.size());

    for(int i=0;i<two_point_obs.size();i++){
        corr[i].resize(T_sites);
        for(int j=0;j<T_sites;j++){
            corr[i][j].resize(T_sites);
        }
    }


    for(int opr_no=0;opr_no<two_point_obs.size();opr_no++){

        cout<<two_point_obs[opr_no]<<":"<<endl;
        if(two_point_obs[opr_no]=="SzSz"){
            mult_type="U";
        }
        else{
            mult_type="not_U";
        }
        for(int site=0;site<T_sites;site++){
            for(int site2=0;site2<T_sites;site2++){

                if(site2>=site){
                    value=zero;


                    Matrix_COO_vector_multiplication(mult_type, Two_point_oprts[opr_no][site][site2], Eig_vecs[state_no], temp_vec);


                    value = dot_product(temp_vec,Eig_vecs[state_no]);


                    corr[opr_no][site][site2]=value;
                    corr[opr_no][site2][site]=value;

                    cout<<value<<"  ";

                }
                else{
                    cout<<zero<<"  ";
                }
            }
            cout<<endl;

        }

        double_type sum_all;
        sum_all=zero;
        for(int site=0;site<T_sites;site++){
            for(int site2=0;site2<T_sites;site2++){
                sum_all = sum_all + corr[opr_no][site][site2];
            }
        }
        cout<<"sum_all = "<<sum_all<<endl;




    }


    cout<<"Now printing SS correlation in k space:"<<endl;

    cout<<"k                              "<<"sin(i*PI*k)*sin(j*PI*k),       cos((i-j)*PI*k),               cos(i*PI*k)*cos(j*PI*k)"<<endl;
    for(int n=1;n<=T_sites;n++){

        double k;
        if(PBC_check==false){
            k=(n*1.0)/(1.0*(T_sites+1));}
        else{
            k=(2*(n-1)*1.0)/(T_sites*1.0);
        }

        value=zero;
        value1=zero;
        value2=zero;
        for(int i=0;i<T_sites;i++){
            for(int j=0;j<T_sites;j++){
                value = value + (3.0/T_sites)*sin(i*PI*k)*sin(j*PI*k)*corr[0][i][j];

                value1 = value1 + (3.0/T_sites)*cos((i-j)*PI*k)*corr[0][i][j];

                value2 = value2 + (3.0/T_sites)*cos(i*PI*k)*cos(j*PI*k)*corr[0][i][j];

            }

        }
        cout<<k<<"  "<<value<<"  "<<value1<<"   "<<value2<<endl;

    }




}

template <typename Basis_type, typename Model_type>
void LANCZOS<Basis_type, Model_type>::Write_full_spectrum(){

    ofstream file_out(file_eval_spectrum.c_str());

    file_out<<"#----> Eval spectrum for each Lanczos iteration"<<endl;

    for(int i=0;i<Evals_Tri_all.size();i++){
        for(int j=0;j<Evals_Tri_all[i].size();j++){
            file_out<<scientific<<setprecision(10)<<i<<"   "<<Evals_Tri_all[i][j]<<endl;

        }
        file_out<<endl;

    }

}

template <typename Basis_type, typename Model_type>
void LANCZOS<Basis_type, Model_type>::Read_Lanczos_parameters(string filename){

    omega_sign =1.0;


    string random_seed_generator_, Random_Seed_Generator_ = "Random_seed_generator = ";
    string eta_, Eta_ = "eta = ";
    string check_omega_, Check_Omega_ = "check_omega = ";
    string omega_max_, Omega_Max_ = "omega_max = ";
    string omega_min_, Omega_Min_ = "omega_min = ";
    string d_omega_, D_Omega_ = "d_omega = ";

    string lanc_error_, Lanc_Error_ = "Lanc_Error = ";
    string max_steps_, Max_Steps_ = "max_steps = ";
    string few__, Few__ = "few = ";
    string no_basis_to_check_, No_Basis_To_Check_= "calculate_local_obs_for_states_to_look_for_first_n_basis = ";
    string get_full_spectrum_, Get_Full_Spectrum_ = "Get_Full_Spectrum = ";
    string read_the_seed_, Read_The_Seed_ = "Read_the_seed = ";
    string File_Eval_Epectrum_ = "file_eval_spectrum = ";
    string Seed_File_ = "seed_file = ";
    string Loc_Obs_In_Basis_Of_States_ = "file_out_local_obs_for_states_to_look = ";
    string File_Dynamics_Out_ = "file_dynamics_out = ";
    string save_all_krylov_space_vecs_, Save_All_Krylov_Space_Vecs_ = "save_all_Krylov_space_vecs = ";
    string get_overlap_bw_eigstates_and_k_vecs_, Get_Overlap_Bw_Eigstates_And_K_Vecs_ = "get_overlap_bw_eigstates_and_K_vecs = ";
    string calculate_local_obs_for_states_to_look_, Calculate_Local_Obs_For_States_To_Look_ = "calculate_local_obs_for_states_to_look = ";
    string lanczos_reorthogonalization_, Lanczos_Reorthogonalization_ = "lanczos_reorthogonalization = ";
    string saving_hamil_, Saving_Hamil_ = "Saving_Hamiltonian_and_Oprs = ";
    string need_few_eig_vecs_ ,Need_Few_Eig_Vecs_ = "need_few_eig_vecs = ";
    string get_overlap_with_basis_, Get_Overlap_With_Basis_ = "Get_overlap_with_basis = ";
    string states_to_look_, States_To_Look_ = "States_to_look = ";
    string Overlap_Out_File_ = "overlap_out_file = ";

    string temperature_range_ftlm_, Temperature_Range_FTLM_ = "Temperature_Range = ";
    string total_random_states_ftlm_, Total_Random_States_FTLM_ = "Number_of_Random_States = ";
    string m_ftlm_, M_FTLM_ = "Number_of_LowKrylov_States_used = ";
    string energy_offset_ftlm_, Energy_Offset_FTLM_ = "Energy_Offset_FTLM = ";
    string temperature_LTLM_dynamics_, Temperature_LTLM_Dynamics_ = "Temperature_LTLM_Dynamics = ";

    int offset;
    string line;
    ifstream inputfile(filename.c_str());


    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);


            if ((offset = line.find(Energy_Offset_FTLM_, 0)) != string::npos) {
                energy_offset_ftlm_ = line.substr (offset + Energy_Offset_FTLM_.length());		}


            if ((offset = line.find(Temperature_LTLM_Dynamics_, 0)) != string::npos) {
                temperature_LTLM_dynamics_ = line.substr (offset + Temperature_LTLM_Dynamics_.length());		}

            if ((offset = line.find(M_FTLM_, 0)) != string::npos) {
                m_ftlm_ = line.substr (offset + M_FTLM_.length());		}

            if ((offset = line.find(Total_Random_States_FTLM_, 0)) != string::npos) {
                total_random_states_ftlm_ = line.substr (offset + Total_Random_States_FTLM_.length());		}

            if ((offset = line.find(Temperature_Range_FTLM_ , 0)) != string::npos) {
                temperature_range_ftlm_ = line.substr (offset + Temperature_Range_FTLM_ .length());		}

            if ((offset = line.find(Eta_, 0)) != string::npos) {
                eta_ = line.substr (offset + Eta_.length());		}

            if ((offset = line.find(Check_Omega_, 0)) != string::npos) {
                check_omega_ = line.substr (offset + Check_Omega_.length());		}

            if ((offset = line.find(Omega_Max_, 0)) != string::npos) {
                omega_max_= line.substr (offset + Omega_Max_.length());		}

            if ((offset = line.find(Omega_Min_, 0)) != string::npos) {
                omega_min_ = line.substr (offset + Omega_Min_.length());		}

            if ((offset = line.find(D_Omega_, 0)) != string::npos) {
                d_omega_= line.substr (offset + D_Omega_.length());		}

            if ((offset = line.find(Lanc_Error_, 0)) != string::npos) {
                lanc_error_ = line.substr (offset + Lanc_Error_.length());		}

            if ((offset = line.find(Max_Steps_, 0)) != string::npos) {
                max_steps_ = line.substr (offset + Max_Steps_.length());		}

            if ((offset = line.find(Few__, 0)) != string::npos) {
                few__ = line.substr (offset + Few__.length());		}

            if ((offset = line.find(Random_Seed_Generator_, 0)) != string::npos) {
                random_seed_generator_ = line.substr (offset + Random_Seed_Generator_.length());		}

            if ((offset = line.find(No_Basis_To_Check_, 0)) != string::npos) {
                no_basis_to_check_ = line.substr (offset + No_Basis_To_Check_.length());		}

            if ((offset = line.find(File_Eval_Epectrum_, 0)) != string::npos) {
                file_eval_spectrum = line.substr (offset + File_Eval_Epectrum_.length());		}

            if ((offset = line.find(Seed_File_, 0)) != string::npos) {
                seed_file_in = line.substr (offset + Seed_File_.length());		}

            if ((offset = line.find(Loc_Obs_In_Basis_Of_States_, 0)) != string::npos) {
                file_Loc_obs_in_basis_of_states = line.substr (offset + Loc_Obs_In_Basis_Of_States_.length());		}

            if ((offset = line.find(Overlap_Out_File_, 0)) != string::npos) {
                overlap_out_file = line.substr (offset + Overlap_Out_File_.length());		}

            if ((offset = line.find(File_Dynamics_Out_, 0)) != string::npos) {
                file_dynamics_out = line.substr (offset + File_Dynamics_Out_.length());		}

            if ((offset = line.find(Save_All_Krylov_Space_Vecs_, 0)) != string::npos) {
                save_all_krylov_space_vecs_ = line.substr (offset + Save_All_Krylov_Space_Vecs_.length());		}

            if ((offset = line.find(Get_Overlap_Bw_Eigstates_And_K_Vecs_, 0)) != string::npos) {
                get_overlap_bw_eigstates_and_k_vecs_ = line.substr (offset + Get_Overlap_Bw_Eigstates_And_K_Vecs_.length());
            }

            if ((offset = line.find(Calculate_Local_Obs_For_States_To_Look_, 0)) != string::npos) {
                calculate_local_obs_for_states_to_look_ = line.substr (offset + Calculate_Local_Obs_For_States_To_Look_.length());
            }

            if ((offset = line.find(Get_Overlap_With_Basis_, 0)) != string::npos) {
                get_overlap_with_basis_ = line.substr (offset + Get_Overlap_With_Basis_.length()); }

            if ((offset = line.find(States_To_Look_, 0)) != string::npos) {
                states_to_look_ = line.substr (offset + States_To_Look_.length()); }

            if ((offset = line.find(Lanczos_Reorthogonalization_, 0)) != string::npos) {
                lanczos_reorthogonalization_ = line.substr (offset + Lanczos_Reorthogonalization_.length());		}

            if ((offset = line.find(Saving_Hamil_, 0)) != string::npos) {
                saving_hamil_ = line.substr (offset + Saving_Hamil_.length());		}

            if ((offset = line.find(Get_Full_Spectrum_, 0)) != string::npos) {
                get_full_spectrum_ = line.substr (offset + Get_Full_Spectrum_.length());		}

            if ((offset = line.find(Read_The_Seed_, 0)) != string::npos) {
                read_the_seed_ = line.substr (offset + Read_The_Seed_.length());		}//here

            if ((offset = line.find(Need_Few_Eig_Vecs_, 0)) != string::npos) {
                need_few_eig_vecs_ = line.substr (offset + Need_Few_Eig_Vecs_.length());		}

        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the Model class."<<endl;}


    eta=atof(eta_.c_str());
    check_omega=atof(check_omega_.c_str());
    omega_max=atof(omega_max_.c_str());
    omega_min=atof(omega_min_.c_str());
    d_omega=atof(d_omega_.c_str());


    Total_Random_States_for_FTLM=atoi(total_random_states_ftlm_.c_str());
    M_FTLM=atoi(m_ftlm_.c_str());
    Energy_Offset_FTLM=atof(energy_offset_ftlm_.c_str());
    Temperature_LTLM_Dynamics=atof(temperature_LTLM_dynamics_.c_str());
    Lanc_Error=atof(lanc_error_.c_str());
    max_steps=atoi(max_steps_.c_str());
    few_=atoi(few__.c_str());
    Random_seed_value=atoi(random_seed_generator_.c_str());
    no_basis_to_check=atoi(no_basis_to_check_.c_str());


    if(read_the_seed_=="true"){
        Read_the_seed=true;
    }
    else{
        Read_the_seed=false;
    }


    if(save_all_krylov_space_vecs_=="true"){
        save_all_Krylov_space_vecs=true;
    }
    else{
        save_all_Krylov_space_vecs=false;
    }


    if(calculate_local_obs_for_states_to_look_=="true"){
        calculate_local_obs_for_states_to_look=true;
    }
    else{
        calculate_local_obs_for_states_to_look=false;
    }


    if(get_overlap_with_basis_=="true"){
        get_overlap_with_basis=true;
    }
    else{
        get_overlap_with_basis=false;
    }

    stringstream states_to_look_ss;
    int n_states;
    states_to_look_ss<<states_to_look_;
    states_to_look_ss>>n_states;
    states_to_look.resize(n_states);

    for(int i=0;i<n_states;i++){
        states_to_look_ss>>states_to_look[i];
    }


    stringstream temperature_range_ftlm_ss;
    temperature_range_ftlm_ss << temperature_range_ftlm_;
    temperature_range_ftlm_ss >> delta_Temperature_FTLM;
    temperature_range_ftlm_ss >> Temprature_min_FTLM;
    temperature_range_ftlm_ss >> Temprature_max_FTLM;



    if(get_overlap_bw_eigstates_and_k_vecs_=="true"){
        get_overlap_bw_eigstates_and_K_vecs=true;
    }
    else{
        get_overlap_bw_eigstates_and_K_vecs=false;
    }

    if(lanczos_reorthogonalization_=="true"){
        lanczos_reorthogonalization=true;
    }
    else{
        lanczos_reorthogonalization=false;
    }

    if(saving_hamil_=="true"){
        Saving_Hamil=true;
    }
    else{
        Saving_Hamil=false;
    }


    if(need_few_eig_vecs_=="true"){
        need_few_eig_vecs=true;
        Get_Full_Spectrum=true;
    }
    else{
        need_few_eig_vecs=false;
    }

    if(get_full_spectrum_=="true"){
        Get_Full_Spectrum=true;
    }
    else{
        Get_Full_Spectrum=false;
    }

    //DO NOT CHANGE THESE, IT SHOULD BE CHANGED ONLY FROM FTLM_DYNAMICS OR (FTLM_STATIC)
    Get_insitu_FTLM_overlaps=false;
    Save_the_Seed=false;
    Get_SeedVec_from_another_routine=false;
    Eig_vecs_required=true;


}

template <typename Basis_type, typename Model_type>
void LANCZOS<Basis_type, Model_type>::Get_Dynamics_seed(Matrix_COO &Opr){

    Matrix_COO_vector_multiplication("not_U", Opr, Eig_vec, Dynamics_seed);

    double tmpnrm_type_double;
    tmpnrm_type_double=Norm(Dynamics_seed);
    Numerator_Dynamics=tmpnrm_type_double;
    double tmpnrm;
    tmpnrm=sqrt(tmpnrm_type_double);


    for(int j=0;j<Dynamics_seed.size();j++){
        Dynamics_seed[j] = (Dynamics_seed[j]/(tmpnrm));

    }
}

template <typename Basis_type, typename Model_type>
void LANCZOS<Basis_type, Model_type>::Get_Dynamics_seed(Mat_1_doub vec_){

    Dynamics_seed=vec_;
    double tmpnrm_type_double;
    tmpnrm_type_double=Norm(Dynamics_seed);
    Numerator_Dynamics=tmpnrm_type_double;
    double tmpnrm;
    tmpnrm=sqrt(tmpnrm_type_double);

    for(int j=0;j<vec_.size();j++){
        Dynamics_seed[j] = (Dynamics_seed[j]*(1.0/(tmpnrm)));

    }

}




//#endif
