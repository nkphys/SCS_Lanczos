#include "TimeDependentLanczos.h"
#include <stdlib.h>

#define PI 3.14159265
using namespace std;
//#define USE_COMPLEX
//#ifdef USE_COMPLEX


#ifndef TimeDependentLanczos_functions
#define TimeDependentLanczos_functions

void TD_Lanczos::reading_input(){
string filepath = inp_filename;

// string model_name = "Model = ";
// string restricted_basis_, Restricted_Basis_ = "Restricted_Basis = ";
// string processors_, Processors_ = "Processors = ";
// int offset;
// string line;

ifstream inputfile(filepath.c_str());

    string evolution_type, Evolution_Type_ = "Evolution_Type = ";
    string timemax, TimeMax_ = "TimeMax = ";
    string no_of_timeslices, No_of_TimeSlices_ = "No_of_TimeSlices = ";
    string constmodelname, ConstModelName_= "Model = ";


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

             if ((offset = line.find(No_of_TimeSlices_, 0)) != string::npos) {
                no_of_timeslices = line.substr (offset+No_of_TimeSlices_.length());    }

            

        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the Model class."<<endl;}

    EvolutionType=evolution_type;
    if(evolution_type != "WithConstModelType"){
        cout<<"At present only Evolution_Type = WithConstModelType is allowed"<<endl;
        assert(false);
    }

    Time_max = atof(timemax.c_str());
    No_TimeSlices = atoi(no_of_timeslices.c_str());
    Time_slice_width= Time_max/(No_TimeSlices);



    if(evolution_type == "WithConstModelType"){

        if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);


            if ((offset = line.find(ConstModelName_, 0)) != string::npos) {
                constmodelname = line.substr (offset+ConstModelName_.length());	}
            
        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the Model class."<<endl;}

    ConstModelName_str=constmodelname;
    }


}

void TD_Lanczos::Constructing_State_T0_route(){

if(State_T0_construction_route=="ViaHamiltonian"){
Constructing_State_T0_via_Hamiltonian();
}
else if(State_T0_construction_route=="GivenAnsatz"){
Constructing_State_T0_via_GivenAnsatz();
}
   
}


void TD_Lanczos::Constructing_State_T0_via_Hamiltonian(){

}

void TD_Lanczos::Constructing_State_T0_via_GivenAnsatz(){

}


void TD_Lanczos::Perform_TD_Lanczos(){

if(ConstModelName_str=="SpinOnly"){

    int time_value=0;
    for(int time_slice=0;time_slice<=No_TimeSlices;time_slice++){

        
        MODEL_Spins _MODEL;
        BASIS_Spins _BASIS;

        _MODEL.Extenstion_to_FilePaths="Timeslice"+to_string(time_slice)+".txt";
        _MODEL.Read_parameters(_BASIS, inp_filename);
        _BASIS.Construct_basis();


        _MODEL.Add_diagonal_terms(_BASIS);
        _MODEL.Add_non_diagonal_terms(_BASIS);
        _MODEL.Add_connections(_BASIS);


        

        time_value += Time_slice_width;
    }


}

}

#endif
