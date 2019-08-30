#include "Lanczos_engine.h"
#include <stdlib.h>

#define PI 3.14159265
using namespace std;
//#define USE_COMPLEX
#ifdef USE_COMPLEX


void LANCZOS::Perform_LANCZOS(Matrix_COO &Hamil){


   // Print_Matrix_COO(Hamil);

    int seed_lanczos=Random_seed_value;
    int lanc_iter=0;
    double eps=Lanc_Error;
    int max_lanczos_states=min(max_steps,Hamil.nrows);
    double diff_E;
    double temp2, E0, E0_old, temp1;
    double temp2_type_double;
    double_type temp1_type_double,  temp3_type_double, temp4_type_double;

    int Target_state=0;
    Mat_1_real B2 ,Norms;
    Mat_1_doub A, red_eig_vec;
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




    srand(seed_lanczos);

    double tmpnrm_type_double,tmpnrm_type_double2;
    double tmpnrm,tmpnrm2;
    if(Dynamics_performed==false){

        if(Read_the_seed==false){
            for(int j=0;j<Hamil.nrows;j++){

                temp1_type_double.real((rand()%RAND_MAX));
                temp1_type_double.imag((rand()%RAND_MAX));

                temp1_type_double=(temp1_type_double)*(1.0/(RAND_MAX*1.0));

                Kvector_n.push_back(temp1_type_double);

            }

            tmpnrm_type_double=dot_product(Kvector_n,Kvector_n).real();
            tmpnrm=sqrt(tmpnrm_type_double);


            for(int j=0;j<Hamil.nrows;j++){

                Kvector_n[j] = (Kvector_n[j]/(tmpnrm));

            }
        }
        else{
            ifstream infile_seed(seed_file_in.c_str());

            for(int j=0;j<Hamil.nrows;j++){
                infile_seed>>tmpnrm;
                joker_double.real(tmpnrm);

                infile_seed>>tmpnrm;
                joker_double.imag(tmpnrm);

                Kvector_n.push_back(joker_double);
            }

        }


    }
    else{
        Kvector_n=Dynamics_seed;

        tmpnrm_type_double=dot_product(Kvector_n,Kvector_n).real();
        tmpnrm=sqrt(tmpnrm_type_double);


        for(int j=0;j<Hamil.nrows;j++){

            Kvector_n[j] = (Kvector_n[j]/(tmpnrm));

        }
        norm_dynamics=tmpnrm_type_double;
    }




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
            //cout<<"Overlarp b/w 0th and last Krylov vector : "<<dot_product(Krylov_space_vecs[0],Krylov_space_vecs[lanc_iter])<<endl;
        }

        Evals_Tri_all.resize(lanc_iter+1);
        clock_t Lanc_time = clock();

        temp1 = sqrt(dot_product(Kvector_n,Kvector_n).real());//*1.0e-10;

        Norms.push_back(temp1);
        if(lanc_iter==0){B2.push_back(0);}
        else{
            B2.push_back(tmpnrm*tmpnrm);
        }

        clock_t oprt_SB_time = clock();

        Matrix_COO_vector_multiplication("U", Hamil, Kvector_n, Kvector_np1); // saved in K_vector_np1




        cout<<"Time to operate Hamiltonian : "<<double( clock() - oprt_SB_time ) / (double)CLOCKS_PER_SEC<<endl;//cout<<"here"<<endl;

        temp3_type_double = dot_product(Kvector_n, Kvector_np1);
        temp2_type_double = dot_product(Kvector_n, Kvector_n).real();


        temp4_type_double = temp3_type_double/temp2_type_double;


        A.push_back(temp4_type_double);


        Subtract(Kvector_np1, A[lanc_iter], Kvector_n);	//

        if(lanc_iter!=0){

            Subtract(Kvector_np1, sqrt(B2[lanc_iter]), Kvector_nm1);

        }


        //Normalizaton of Knp1, added by myself, not included in std. Lanczos
        tmpnrm_type_double=dot_product(Kvector_np1,Kvector_np1).real();
        tmpnrm=sqrt(tmpnrm_type_double);
        for(int i=0;i<Kvector_np1.size();i++){
            Kvector_np1[i] = (Kvector_np1[i]/(tmpnrm));//*1.0e-10;
        }


        if(save_all_Krylov_space_vecs==true){
            cout<<"Overlarp b/w 0th and Kvector_np1 : "<<dot_product(Krylov_space_vecs[0],Kvector_np1)<<endl;
        }


        if(lanc_iter!=0){
            if(to_reorth==true){
                Mat_1_doub temp_dot;
                temp_dot.resize(Krylov_space_vecs.size());
                for(int i=0;i<Krylov_space_vecs.size();i++){
                    temp_dot[i]=dot_product(Kvector_np1,Krylov_space_vecs[i]);
                }

                for(int l=0;l<Kvector_n.size();l++){
                    val_temp=zero;
                    for(int i=0;i<Krylov_space_vecs.size();i++){
                        val_temp = val_temp + Krylov_space_vecs[i][l]*temp_dot[i];
                    }
                    Kvector_np1[l] = Kvector_np1[l] - val_temp;
                }
                temp_dot.clear();

                //Normalizaton of Knp1, added by myself, not included in std. Lanczos
                tmpnrm_type_double2=dot_product(Kvector_np1,Kvector_np1).real();
                tmpnrm2=sqrt(tmpnrm_type_double2);
               for(int i=0;i<Kvector_np1.size();i++){
                  Kvector_np1[i] = (Kvector_np1[i]/(tmpnrm2));//*1.0e-10;
                }

               cout<<"Reorthogonalization performed"<<endl;
            }
        }

        if(save_all_Krylov_space_vecs==true){
            cout<<"Overlarp b/w 0th and Kvector_np1 : "<<dot_product(Krylov_space_vecs[0],Kvector_np1)<<endl;
        }

        if(Get_Full_Spectrum==false){
            Diagonalize(A,B2,E0,red_eig_vec);
        }
        else{
            Evals_Tri_all[lanc_iter].resize(A.size());
            Diagonalize(A,B2,Evals_Tri_all,red_eig_vecs,lanc_iter,few_,states_to_look);
            red_eig_vec=red_eig_vecs[0];
            E0=Evals_Tri_all[lanc_iter][0];
        }


        diff_E=	fabs(E0-E0_old);
        if(!Dynamics_performed){
            cout<<"iter = "<<lanc_iter<<" diff_E = "<<diff_E<<" E0 = "<<E0<<" E0_old = "<<E0_old<<endl;
            GS_energy=E0;
        }
        if(lanc_iter==0){if(dot_product(Kvector_np1,Kvector_np1)==zero){diff_E = 0;}}
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


    cout<<"Perform_LANCZOS: "<<"LANCZOS(pass 2) STARTING FOR SUPERBLOCK Eigenvector, "<<", Size of Matrix(SB) = "<<Hamil.ncols<<endl<<endl;

    Lanc_iter_done=lanc_iter;

    Kvector_n.clear();	Kvector_nm1.clear(); Kvector_np1.clear();


    if(Dynamics_performed==false){


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
            srand(seed_lanczos);


            if(Read_the_seed==false){
                for(int j=0;j<Hamil.nrows;j++){

                    temp1_type_double.real((rand()%RAND_MAX));
                    temp1_type_double.imag((rand()%RAND_MAX));

                    temp1_type_double=(temp1_type_double)*(1.0/(RAND_MAX*1.0));

                    Kvector_n.push_back(temp1_type_double);

                }

                tmpnrm_type_double=dot_product(Kvector_n,Kvector_n).real();
                tmpnrm=sqrt(tmpnrm_type_double);


                for(int j=0;j<Hamil.nrows;j++){

                    Kvector_n[j] = (Kvector_n[j]/(tmpnrm));

                }
            }
            else{
                ifstream infile_seed(seed_file_in.c_str());

                for(int j=0;j<Hamil.nrows;j++){
                    infile_seed>>tmpnrm;
                    joker_double.real(tmpnrm);

                    infile_seed>>tmpnrm;
                    joker_double.imag(tmpnrm);

                    Kvector_n.push_back(joker_double);
                }

            }


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
                //cout<<"NOrm = "<<Norms[lanc_iter2]<<endl<<endl;

                // saved in K_vector_np1
                Matrix_COO_vector_multiplication("U", Hamil, Kvector_n, Kvector_np1);

                Subtract(Kvector_np1, A[lanc_iter2], Kvector_n);	//
                if(lanc_iter2!=0){
                    Subtract(Kvector_np1, sqrt(B2[lanc_iter2]), Kvector_nm1);
                }

                //Normalizaton of Knp1 , not included in std. Lanczos
                tmpnrm_type_double = dot_product(Kvector_np1,Kvector_np1).real(); //new
                tmpnrm=sqrt(tmpnrm_type_double);

                for(int i=0;i<Kvector_np1.size();i++){
                    Kvector_np1[i] = (Kvector_np1[i]/(tmpnrm));
                }


                Kvector_nm1=Kvector_n;
                Kvector_n=Kvector_np1;

            }
            }

            double norm_ev=dot_product(Eig_vecs[Ts], Eig_vecs[Ts]).real();
            for(int j=0;j<Hamil.nrows;j++){

                Eig_vecs[Ts][j] = (Eig_vecs[Ts][j]/(sqrt(norm_ev)));


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

                    Overlaps[Ts][bi].first=( Eig_vecs[Ts][bi]*conj(Eig_vecs[Ts][bi]) ).real() ;
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
                    outfile_overlap<<"  "<<Overlaps_bare[Ts][bi].first;
                }
                for(int Ts=0;Ts<Eig_vecs.size();Ts++){
                    outfile_overlap<<"  "<<Overlaps_bare[Ts][bi].second;
                }

                outfile_overlap<<endl;
            }




            if(Check_Ghosts==true){
                string out3 = "ghost.txt";
                ofstream file_out3(out3.c_str());

                for(int Ts1=0;Ts1<Eig_vecs.size();Ts1++){

                    for(int Ts2=0;Ts2<Eig_vecs.size();Ts2++){
                        file_out3<<dot_product(Eig_vecs[Ts1],Eig_vecs[Ts2])<<"  ";

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
    B2.clear();A.clear(); red_eig_vec.clear();Norms.clear();


}

void LANCZOS::Measure_macro_observables(Mat_1_string macro_obs, Hamiltonian_1_COO &Macro_oprts, int state_no){

    cout<<"Macro obs Measurement is started after completing Lanczos algorithm"<<endl;
    cout<<"Macro obs Measurement is started after completing Lanczos algorithm fo state_no = "<<states_to_look[state_no]<<endl;
    //Krylov_space_vecs



    double_type value;
    Mat_1_doub temp_vec;

    /*
    for(int opr_no=0;opr_no<one_point_obs.size();opr_no++){

        for(int site=0;site<T_sites;site++){

            value=0;
            for(int i=0;i<Krylov_space_vecs.size();i++){

                Matrix_COO_vector_multiplication("U", One_point_oprts[opr_no][site], Krylov_space_vecs[i], temp_vec);

                for(int j=0;j<Krylov_space_vecs.size();j++){
                    value = value + dot_product(Krylov_space_vecs[j],temp_vec)*red_eig_vecs[state][j]*red_eig_vecs[state][i];
                }

            }

            cout<<one_point_obs[opr_no]<<"["<<site<<"]  "<<value<<endl;


        }

    }
    */

    //Faster way:

    for(int opr_no=0;opr_no<macro_obs.size();opr_no++){


        Matrix_COO_vector_multiplication("U", Macro_oprts[opr_no], Eig_vecs[state_no], temp_vec);

        value = dot_product(Eig_vecs[state_no],temp_vec);

        cout<<macro_obs[opr_no]<<" =  "<<value<<endl;




    }


}




void LANCZOS::Measure_KE(Matrix_COO &Macro_oprt, int state_no){


    double_type value;
    Mat_1_doub temp_vec,vector_used;

    if(state_no==-121){
        vector_used=Evolving_State;
    }
    else{
        vector_used=Eig_vecs[state_no];
    }


    //Faster way:

        Matrix_COO_vector_multiplication("U", Macro_oprt, vector_used, temp_vec);

        value = dot_product(vector_used,temp_vec).real();

        cout<<value<<"#KE"<<endl;




}


void LANCZOS::Measure_Total_Energy(Matrix_COO &Macro_oprt, int state_no){


    double_type value;
    Mat_1_doub temp_vec,vector_used;

    if(state_no==-121){
        vector_used=Evolving_State;
    }
    else{
        vector_used=Eig_vecs[state_no];
    }


    //Faster way:

        Matrix_COO_vector_multiplication("U", Macro_oprt, vector_used, temp_vec);

        value = dot_product(vector_used,temp_vec).real();

        cout<<value<<"#TE"<<endl;

}




void LANCZOS::Time_evolution_type1(Matrix_COO &Hamil, double dt_, double &Energy_){

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
}

void LANCZOS::Measure_one_point_observables(Mat_1_string one_point_obs, Hamiltonian_2_COO &One_point_oprts, int T_sites, int state_no){

    cout<<"One-point Measurement is started after completing Lanczos algorithm"<<endl;
    //Krylov_space_vecs


    double_type value;
    Mat_1_doub temp_vec;
    Mat_1_doub vector_used;

    if(state_no==-121){
        vector_used=Evolving_State;
    }
    else{
        vector_used=Eig_vecs[state_no];
    }

    complex<double> sz;
    complex<double> jz_eff(0,0);
    complex<double> n_up(0,0);
    complex<double> n_dn(0,0);

    /*
    for(int opr_no=0;opr_no<one_point_obs.size();opr_no++){

        for(int site=0;site<T_sites;site++){

            value=0;
            for(int i=0;i<Krylov_space_vecs.size();i++){

                Matrix_COO_vector_multiplication("U", One_point_oprts[opr_no][site], Krylov_space_vecs[i], temp_vec);

                for(int j=0;j<Krylov_space_vecs.size();j++){
                    value = value + dot_product(Krylov_space_vecs[j],temp_vec)*red_eig_vecs[state][j]*red_eig_vecs[state][i];
                }

            }

            cout<<one_point_obs[opr_no]<<"["<<site<<"]  "<<value<<endl;


        }

    }
    */


    //Faster way:

    One_point_observables_values.resize(One_point_oprts.size());
    for(int opr_no=0;opr_no<One_point_oprts.size();opr_no++){
        One_point_observables_values[opr_no].resize(T_sites);

        for(int site=0;site<T_sites;site++){

            Matrix_COO_vector_multiplication("cx", One_point_oprts[opr_no][site], vector_used, temp_vec);

            value = dot_product(temp_vec,vector_used);

            cout<<one_point_obs[opr_no]<<"["<<site<<"]  "<<value<<endl;
            One_point_observables_values[opr_no][site]=value;

            if(one_point_obs[opr_no]=="n_3by2_3by2"){
                jz_eff += (3.0/2.0)*value;

            }
            if(one_point_obs[opr_no]=="n_3by2_m3by2"){
                jz_eff += (-3.0/2.0)*value;

            }
            if(one_point_obs[opr_no]=="n_3by2_1by2" || one_point_obs[opr_no]=="n_1by2_1by2"){
                jz_eff += (1.0/2.0)*value;

            }
            if(one_point_obs[opr_no]=="n_3by2_m1by2" || one_point_obs[opr_no]=="n_1by2_m1by2"){
                jz_eff += (-1.0/2.0)*value;

            }

            if((one_point_obs[opr_no]=="n_0_up" || one_point_obs[opr_no]=="n_1_up") || one_point_obs[opr_no]=="n_2_up"  ){
                n_up += (1.0)*value;

            }
            if((one_point_obs[opr_no]=="n_0_dn" || one_point_obs[opr_no]=="n_1_dn") || one_point_obs[opr_no]=="n_2_dn"  ){
                n_dn += (1.0)*value;

            }




        }

    }

    cout<<"Jz_eff_total = "<<jz_eff<<endl;
    cout<<"N_total = "<<n_up + n_dn <<endl;
    cout<<"Sz_total = "<<0.5*(n_up - n_dn) <<endl;

    vector_used.clear();

}


void LANCZOS::Measure_two_point_observables_smartly(Mat_1_string one_point_obs, Hamiltonian_2_COO &One_point_oprts, int T_sites, int state_no, string _model_){
    cout<<"Two-point Measurement is done smartly = "<<states_to_look[state_no]<<endl;
    //Krylov_space_vecs


    double_type value,value1,value2;
    Mat_1_doub temp_vec, temp_vec2;
    string mult_type;
    Mat_3_doub corr;
    Matrix_COO TEMP_COO, TEMP_COO2;
    Matrix_COO TEMP_COO3, TEMP_COO4;
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


    if(_model_=="3_orb_Hubbard_chain_GC"){

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





    double Total_excitons, Frenkel_excitons;

    cout<<"--------<n_{1/2,1/2}(site)> - <n_{1/2,1/2}(site)n_{3/2,1/2}(site2)>--------------------"<<endl;

    Frenkel_excitons=0;
    Total_excitons=0;
    //n_{3/2,1/2}=8, n_{1/2,1/2}=10
    for(int site=0;site<T_sites;site++){

        Matrix_COO_vector_multiplication("Full", One_point_oprts[10][site], Eig_vecs[state_no], temp_vec);
        value1 = dot_product(temp_vec, Eig_vecs[state_no]);

        for(int site2=0;site2<T_sites;site2++){

                Matrix_COO_vector_multiplication("Full", One_point_oprts[8][site2], Eig_vecs[state_no], temp_vec);
                Matrix_COO_vector_multiplication("Full", One_point_oprts[10][site], temp_vec, temp_vec2);

                value = dot_product(temp_vec2, Eig_vecs[state_no]);

                Total_excitons = Total_excitons + value1.real() - value.real();
                if(site == site2){
                  Frenkel_excitons = Frenkel_excitons + value1.real() - value.real();
                }
                cout<<value1 - value<<"  ";


        }
        cout<<endl;

    }

    cout<<"Frenkel_Exc/Total_Excitons = "<<Frenkel_excitons/Total_excitons<<"   " <<Frenkel_excitons<<"  "<<Total_excitons<<endl;

}

}


double LANCZOS::Measure_observable(Matrix_COO &OPR_, int state_no){}

void LANCZOS::Measure_four_point_observables(Hamiltonian_3_COO &Two_point_oprts, Mat_1_tetra_int sites_set, int state_no){

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

void LANCZOS::Measure_two_point_observables(Mat_1_string two_point_obs, Hamiltonian_3_COO &Two_point_oprts, int T_sites,  int state_no, bool PBC_check){

    cout<<"Two-point Measurement is started after completing Lanczos algorithm for state_no = "<<states_to_look[state_no]<<endl;
    //Krylov_space_vecs


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

void LANCZOS::Write_full_spectrum(){

    ofstream file_out(file_eval_spectrum.c_str());

    file_out<<"#----> Eval spectrum for each Lanczos iteration"<<endl;

    for(int i=0;i<Evals_Tri_all.size();i++){
        for(int j=0;j<Evals_Tri_all[i].size();j++){
            file_out<<scientific<<setprecision(10)<<i<<"   "<<Evals_Tri_all[i][j]<<endl;

        }
        file_out<<endl;

    }

}

void LANCZOS::Read_Lanczos_parameters(string filename){
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
    string need_few_eig_vecs_ ,Need_Few_Eig_Vecs_ = "need_few_eig_vecs = ";
    string get_overlap_with_basis_, Get_Overlap_With_Basis_ = "Get_overlap_with_basis = ";
    string states_to_look_, States_To_Look_ = "States_to_look = ";
    string Overlap_Out_File_ = "overlap_out_file = ";



    int offset;
    string line;
    ifstream inputfile(filename.c_str());


    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);


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

}

void LANCZOS::Get_Dynamics_seed(Matrix_COO &Opr){

    Matrix_COO_vector_multiplication("not_U", Opr, Eig_vec, Dynamics_seed);

    double tmpnrm_type_double=(dot_product(Dynamics_seed,Dynamics_seed)).real();
    Numerator_Dynamics=tmpnrm_type_double;
    double tmpnrm;
    tmpnrm=sqrt(tmpnrm_type_double);


    for(int j=0;j<Eig_vec.size();j++){
        Dynamics_seed[j] = (Dynamics_seed[j]/(tmpnrm));

    }
}

void LANCZOS::Get_Dynamics_seed(Mat_1_doub vec_){

    Dynamics_seed=vec_;
    double tmpnrm_type_double=(dot_product(Dynamics_seed,Dynamics_seed)).real();
    Numerator_Dynamics=tmpnrm_type_double;
    double tmpnrm;
    tmpnrm=sqrt(tmpnrm_type_double);


    for(int j=0;j<vec_.size();j++){
        Dynamics_seed[j] = (Dynamics_seed[j]*(1.0/(tmpnrm)));

    }

}




#endif
