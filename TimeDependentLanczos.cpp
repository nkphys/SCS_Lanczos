#include "TimeDependentLanczos.h"
#include <stdlib.h>

#define PI 3.14159265
using namespace std;
//#define USE_COMPLEX
//#ifdef USE_COMPLEX


#ifndef TimeDependentLanczos_functions
#define TimeDependentLanczos_functions


void TD_Lanczos::Create_Scheduler(){

    //assert(abs(Restart_Time)<0.0000000001);

    cout<<"No of Time slices = "<<No_TimeSlices<<endl;

    Time_.resize(No_TimeSlices);
    Gamma_.resize(No_TimeSlices);
    Js_.resize(No_TimeSlices);

    double time_normalized;
    double gamma_temp, js_temp;
    for(int time_ind=0;time_ind<No_TimeSlices;time_ind++){

        time_normalized=(time_ind*dt_)/(Time_max);

       // cout<<time_ind<<"  "<<time_normalized<<endl;
        //assert(time_normalized>=0 && time_normalized<=1.0);

        for(int time_ind2=0;time_ind2<Time_bare.size()-1;time_ind2++){
            if(time_normalized>=Time_bare[time_ind2] && time_normalized<=Time_bare[time_ind2+1] ){
                gamma_temp = (abs(Time_bare[time_ind2+1]-time_normalized)*Gamma_bare[time_ind2]
                        +abs(Time_bare[time_ind2]-time_normalized)*Gamma_bare[time_ind2+1])*
                        (1.0/(Time_bare[time_ind2+1]-Time_bare[time_ind2]));

                js_temp = (abs(Time_bare[time_ind2+1]-time_normalized)*Js_bare[time_ind2]
                        +abs(Time_bare[time_ind2]-time_normalized)*Js_bare[time_ind2+1])*
                        (1.0/(Time_bare[time_ind2+1]-Time_bare[time_ind2]));

                break;
            }
        }

        Time_[time_ind]=time_ind*dt_;
        Gamma_[time_ind]=gamma_temp;
        Js_[time_ind]=js_temp;
    }


    string created_schd_str = "Created_Scheduler.txt";
    ofstream created_schd_stream(created_schd_str.c_str());

    created_schd_stream<<"# time   Gamma    Js"<<endl;
    for(int time_ind=0;time_ind<Time_.size();time_ind++){
        created_schd_stream<<Time_[time_ind]<<"  "<<Gamma_[time_ind]<<"  "<<Js_[time_ind]<<endl;
    }

}



void TD_Lanczos::reading_input(){
    string filepath = inp_filename;


    ifstream inputfile(filepath.c_str());

    string evolution_type, Evolution_Type_ = "Evolution_Type = ";
    string timemax, TimeMax_ = "TimeMax = ";
    //string no_of_timeslices, No_of_TimeSlices_ = "No_of_TimeSlices = ";
    string time_otoc, Time_Otoc = "Time_otoc = ";

    string krylovstepsfortimeevo_, KrylovStepsForTimeEvo_ = "KrylovStepsForTimeEvo = ";
    string m_eigStates_, M_EigStates_ = "M_EigStates = ";

    string sitev_, SiteV_ = "SiteV = ";
    string sitew_, SiteW_ = "SiteW = ";

    string time_slice_width, Time_slice_width_ = "Time_slice_width = ";
    string modelname, ModelName_= "Model = ";


    string Scheduler_File_ = "Scheduler_File = ";
    string use_scheduler_, Use_Scheduler_ = "Use_Scheduler = ";


    int offset;
    string line;

    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);


            if ((offset = line.find(Evolution_Type_, 0)) != string::npos) {
                evolution_type = line.substr (offset+Evolution_Type_.length());	}

            if ((offset = line.find(TimeMax_, 0)) != string::npos) {
                timemax = line.substr (offset+TimeMax_.length());    }

            if ((offset = line.find(Time_Otoc, 0)) != string::npos) {
                time_otoc = line.substr (offset+Time_Otoc.length());    }

            if ((offset = line.find(KrylovStepsForTimeEvo_, 0)) != string::npos) {
                krylovstepsfortimeevo_ = line.substr (offset+KrylovStepsForTimeEvo_.length());    }

            if ((offset = line.find(M_EigStates_, 0)) != string::npos) {
                m_eigStates_ = line.substr (offset+M_EigStates_.length());    }

            if ((offset = line.find(SiteV_, 0)) != string::npos) {
                sitev_ = line.substr (offset+SiteV_.length());    }

            if ((offset = line.find(SiteW_, 0)) != string::npos) {
                sitew_ = line.substr (offset+SiteW_.length());    }

            if ((offset = line.find(Time_slice_width_, 0)) != string::npos) {
                time_slice_width = line.substr (offset+Time_slice_width_.length());    }

            if ((offset = line.find(ModelName_, 0)) != string::npos) {
                modelname = line.substr (offset+ModelName_.length());	}

            if ((offset = line.find(Scheduler_File_, 0)) != string::npos) {
                Scheduler_File = line.substr (offset + Scheduler_File_.length());		}

            if ((offset = line.find(Use_Scheduler_, 0)) != string::npos) {
                use_scheduler_ = line.substr (offset + Use_Scheduler_.length());		}
            
        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the Model class."<<endl;}

    EvolutionType=evolution_type;
    //    if(evolution_type != "WithConstModelType"){
    //        cout<<"At present only Evolution_Type = WithConstModelType is allowed"<<endl;
    //        assert(false);
    //    }

    Time_max = atof(timemax.c_str());
    Time_otoc =atof(time_otoc.c_str());
    //No_TimeSlices = atoi(no_of_timeslices.c_str());
    KrylovStepsForTimeEvo = atoi(krylovstepsfortimeevo_.c_str());
    M_EigStates = atoi(m_eigStates_.c_str());
    SiteV = atoi(sitev_.c_str());
    SiteW = atoi(sitew_.c_str());
    Time_slice_width= atof(time_slice_width.c_str());
    dt_ = Time_slice_width;
    No_TimeSlices = int(Time_max/Time_slice_width)+1;


    //ModelName_str=constmodelname;


    if(use_scheduler_ == "true"){
        Use_Scheduler=true;
        cout<<"Scheduler is used i.e. time dependent Hamiltonian"<<endl;
    }
    else{
        Use_Scheduler=false;
    }


    if(Use_Scheduler){
    string line2;
    double temp_t, temp_h, temp_J;
    ifstream scheduler_stream(Scheduler_File.c_str());
    while(getline(scheduler_stream,line2)){
    stringstream line_ss(line2);
    line_ss>>temp_t>>temp_h>>temp_J;
    Time_bare.push_back(temp_t);
    Gamma_bare.push_back(temp_h);
    Js_bare.push_back(temp_J);
    }
    }


   // cout<<"Here 1"<<endl;

}

void TD_Lanczos::Constructing_InitialState(){

    State_T0_construction_route="HamiltonianGS";

    if(State_T0_construction_route=="HamiltonianGS"){
        Constructing_State_T0_via_HamiltonianGS();
    }
    else if(State_T0_construction_route=="GivenAnsatz"){
        Constructing_State_T0_via_GivenAnsatz();
    }

}


void TD_Lanczos::Constructing_State_T0_via_HamiltonianGS(){

    MODEL_Spins _MODEL;
    BASIS_Spins _BASIS;


    _MODEL.Read_parameters(_BASIS, inp_filename);
    _BASIS.Construct_basis();

    //double Hx_factor,double Hz_factor, double Jpm_factor, double Jzz_factor
    _MODEL.Update_Hamiltonian_Params(_BASIS, Gamma_[0], 1.0, 1.0, Js_[0]);

   // _MODEL.Update_Hamiltonian_Params(_BASIS, 1.0, 1.0, 1.0, 1.0);

    _MODEL.Add_diagonal_terms(_BASIS);
    _MODEL.Add_non_diagonal_terms(_BASIS);
    _MODEL.Add_connections(_BASIS);

    LANCZOS<BASIS_Spins, MODEL_Spins> _LANCZOS_GS(_BASIS, _MODEL);
    _LANCZOS_GS.Dynamics_performed=false;
    _LANCZOS_GS.Read_Lanczos_parameters(inp_filename);
    _LANCZOS_GS.TimeEvoPerformed=false;
    _LANCZOS_GS.Perform_LANCZOS(_MODEL.Hamil);
    Print_vector_in_file(_LANCZOS_GS.Eig_vec,"seed_GS.txt");

    cout<<"--HERE 0------"<<endl;

    Psi0 = _LANCZOS_GS.Eig_vec;

    cout<<"--GS observables--------"<<endl;
    _MODEL.MeasureLocalOprs(_BASIS, _LANCZOS_GS.Eig_vec);
    _MODEL.MeasureTwoPointOprs(_BASIS, _LANCZOS_GS.Eig_vec);

}

void TD_Lanczos::Constructing_State_T0_via_GivenAnsatz(){

}

//Act_Operator(BASIS_Spins &basis, Mat_1_doub &Vec_in, Mat_1_doub &Vec_out, string opr_str, int opr_site)

void TD_Lanczos::Calculate_OTOC(string opr_type, int siteV, int siteW){

    double_type WtV, Otoc;

    MODEL_Spins _MODEL;
    BASIS_Spins _BASIS;
    _MODEL.Read_parameters(_BASIS, inp_filename);
    _BASIS.Construct_basis();

    Mat_1_doub Psi2, Psi3, Psi1;

   int time_slice_otoc = int(Time_otoc/Time_slice_width);

   cout<<"time_slice_otoc = "<<time_slice_otoc<<endl;

    //Act V Opr at site=siteV
    _MODEL.Act_Operator(_BASIS, Psi0, Psi2, opr_type, siteV);

   //t=0---->t_otoc
   Perform_TimeEvolution(Psi2, Psi3, 1, time_slice_otoc);

   //Act W Opr at site=siteW
   _MODEL.Act_Operator(_BASIS, Psi3, Psi2, opr_type, siteW);

   //t_otoc----->t=0
   Perform_TimeEvolution(Psi2, Psi3, time_slice_otoc-1, 0);

    WtV=dot_product(Psi3,Psi0); //<Psi0|Psi3>

    cout<<"WtV = "<<WtV<<endl;



    //t=0---->t_otoc
    Perform_TimeEvolution(Psi0, Psi2, 1, time_slice_otoc);

    //Act W Opr at site=siteW
    _MODEL.Act_Operator(_BASIS, Psi2, Psi1, opr_type, siteW);

    //t_otoc----->t=0
    Perform_TimeEvolution(Psi1, Psi2, time_slice_otoc-1, 0);

    //Act V Opr at site=siteV
    _MODEL.Act_Operator(_BASIS, Psi2, Psi1, opr_type, siteV);


     Otoc=dot_product(Psi3,Psi1); //<Psi1|Psi3>

     cout<<"OTOC = "<<Otoc<<endl;

}



void TD_Lanczos::Perform_TimeEvolution(Mat_1_doub &Psi_initial, Mat_1_doub &Psi_final, int time_slice_init,
                                       int time_slice_final){


    int time_arrow = (time_slice_final - time_slice_init)/(abs(time_slice_final - time_slice_init));

    assert(abs(time_arrow)==1);

    int time_value=0;
    double_type Energy;
    cout<<"Starting time evolution------------"<<endl;

    Mat_1_doub Psi_old = Psi_initial;

    int time_slice = time_slice_init;


    //for(int time_slice=time_slice_init;time_slice<=No_TimeSlices;time_slice++){
    while(
          ((time_arrow>0) && (time_slice<=time_slice_final) && (time_slice>=time_slice_init))
          ||
          ((time_arrow<0) && (time_slice>=time_slice_final) && (time_slice<=time_slice_init))

          ){


        MODEL_Spins _MODEL;
        BASIS_Spins _BASIS;

        //_MODEL.Extenstion_to_FilePaths="Timeslice"+to_string(time_slice)+".txt";
        _MODEL.Read_parameters(_BASIS, inp_filename);
        _BASIS.Construct_basis();

        _MODEL.Update_Hamiltonian_Params(_BASIS, Gamma_[time_slice], 1.0, 1.0, Js_[time_slice]);

        _MODEL.Add_diagonal_terms(_BASIS);
        _MODEL.Add_non_diagonal_terms(_BASIS);
        _MODEL.Add_connections(_BASIS);

        LANCZOS<BASIS_Spins, MODEL_Spins> _LANCZOS(_BASIS, _MODEL);
        _LANCZOS.Dynamics_performed=false;
        _LANCZOS.Read_Lanczos_parameters(inp_filename);

        cout<<"Here 1"<<endl;
        M_=min(M_EigStates,_BASIS.basis_size);
        _LANCZOS.TimeEvoPerformed=true;
        _LANCZOS.M_TimeEvo= M_;
        _LANCZOS.dt_TimeEvo = dt_*PI*time_arrow;
        _LANCZOS.Get_SeedVec_from_another_routine=true;
        _LANCZOS.Seed_used = Psi_old;
        _LANCZOS.save_all_Krylov_space_vecs=true; // increases RAM, but speed up by 2
        _LANCZOS.lanczos_reorthogonalization=true;
        _LANCZOS.Eig_vecs_required=true;
        _LANCZOS.Get_Full_Spectrum=true;
        _LANCZOS.Lanc_Error=-10.0;
        _LANCZOS.max_steps=min(KrylovStepsForTimeEvo,_BASIS.basis_size);
        _LANCZOS.states_to_look.resize(M_);
        for(int l=0;l<M_;l++){
         _LANCZOS.states_to_look[l]=l;
        }

        _LANCZOS.Perform_LANCZOS(_MODEL.Hamil);


        assert(_LANCZOS.Evals_Tri_all.size()==_LANCZOS.Lanc_iter_done);

        Psi_old = _LANCZOS.Vec_new_TimeEvo;

        cout<<"Obs. for time = "<<time_slice*dt_<<endl;
        if(time_slice%1==0){
            _MODEL.MeasureLocalOprs(_BASIS, Psi_old);
            _MODEL.MeasureTwoPointOprs(_BASIS, Psi_old);
            _MODEL.MeasureEnergy(_BASIS, Psi_old, Energy);
            cout<<time_slice<<"  "<<time_slice*dt_<<" <Psi(t)|H|Psi(t)> : "<<Energy<<endl;
        }

        time_value += Time_slice_width;
        time_slice +=time_arrow;
    }

Psi_final=Psi_old;


}

void TD_Lanczos::Perform_TD_Lanczos(){

        int time_value=0;
        double_type Energy;
        cout<<"Starting time evolution------------"<<endl;

        Mat_1_doub Psi_old = Psi0;
        for(int time_slice=1;time_slice<=No_TimeSlices;time_slice++){


            MODEL_Spins _MODEL;
            BASIS_Spins _BASIS;


            //_MODEL.Extenstion_to_FilePaths="Timeslice"+to_string(time_slice)+".txt";
            _MODEL.Read_parameters(_BASIS, inp_filename);
            _BASIS.Construct_basis();

            _MODEL.Update_Hamiltonian_Params(_BASIS, Gamma_[time_slice], 1.0, 1.0, Js_[time_slice]);

            _MODEL.Add_diagonal_terms(_BASIS);
            _MODEL.Add_non_diagonal_terms(_BASIS);
            _MODEL.Add_connections(_BASIS);

            LANCZOS<BASIS_Spins, MODEL_Spins> _LANCZOS(_BASIS, _MODEL);
            _LANCZOS.Dynamics_performed=false;
            _LANCZOS.Read_Lanczos_parameters(inp_filename);

            cout<<"Here 1"<<endl;
            M_=min(100,_BASIS.basis_size);
            _LANCZOS.TimeEvoPerformed=true;
            _LANCZOS.M_TimeEvo= M_;
            _LANCZOS.dt_TimeEvo = dt_*PI;
            _LANCZOS.Get_SeedVec_from_another_routine=true;
            _LANCZOS.Seed_used = Psi_old;
            _LANCZOS.save_all_Krylov_space_vecs=true; // increases RAM, but speed up by 2
            _LANCZOS.lanczos_reorthogonalization=true;
            _LANCZOS.Eig_vecs_required=true;
            _LANCZOS.Get_Full_Spectrum=true;
            _LANCZOS.Lanc_Error=-10.0;
            _LANCZOS.max_steps=min(100,_BASIS.basis_size);
            _LANCZOS.states_to_look.resize(M_);
            for(int l=0;l<M_;l++){
             _LANCZOS.states_to_look[l]=l;
            }

            _LANCZOS.Perform_LANCZOS(_MODEL.Hamil);


            assert(_LANCZOS.Evals_Tri_all.size()==_LANCZOS.Lanc_iter_done);

            Psi_old = _LANCZOS.Vec_new_TimeEvo;

            cout<<"Obs. for time = "<<time_slice*dt_<<endl;
            if(time_slice%1==0){
                _MODEL.MeasureLocalOprs(_BASIS, Psi_old);
                _MODEL.MeasureTwoPointOprs(_BASIS, Psi_old);
                _MODEL.MeasureEnergy(_BASIS, Psi_old, Energy);
                cout<<time_slice<<"  "<<time_slice*dt_<<" <Psi(t)|H|Psi(t)> : "<<Energy<<endl;
            }

            time_value += Time_slice_width;
        }



}

#endif
