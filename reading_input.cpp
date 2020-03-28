#include "reading_input.h"

void reading_input(string filename, string & Model_Name, bool &Do_Dynamics, bool &Static_Finite_Temp, bool &Dynamics_Finite_Temp, bool &Perform_RIXS, bool &Restricted_Basis){


    string filepath = filename;
    string model_name = "Model = ";
    string do_dynamics_,Do_Dynamics_ = "Perform_Dynamics = ";
    string do_static_finite_temp_, Static_Finite_Temp_ = "Static_Finite_Temperature = ";
    string do_dynamics_finite_temp_, Dynamics_Finite_Temp_ = "Dynamics_Finite_Temperature = ";
    string perform_rixs_, Perform_RIXS_ = "Perform_RIXS = ";
    string restricted_basis_, Restricted_Basis_ = "Restricted_Basis = ";
    int offset;
    string line;


    ifstream inputfile(filepath.c_str());


    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);



            if ((offset = line.find(model_name, 0)) != string::npos) {
                Model_Name = line.substr (offset+model_name.length());				}

            if ((offset = line.find(Do_Dynamics_, 0)) != string::npos) {
                do_dynamics_ = line.substr (offset+Do_Dynamics_.length());				}

            if ((offset = line.find(Restricted_Basis_, 0)) != string::npos) {
                restricted_basis_ = line.substr (offset+Restricted_Basis_.length());	}

            if ((offset = line.find(Perform_RIXS_ , 0)) != string::npos) {
                perform_rixs_ = line.substr (offset+Perform_RIXS_ .length());				}

            if ((offset = line.find(Static_Finite_Temp_ , 0)) != string::npos) {
                do_static_finite_temp_ = line.substr (offset+Static_Finite_Temp_ .length());	}

            if ((offset = line.find(Dynamics_Finite_Temp_ , 0)) != string::npos) {
                do_dynamics_finite_temp_ = line.substr (offset+Dynamics_Finite_Temp_ .length());	}

        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file 1st time."<<endl;}

if(do_dynamics_=="true"){
    Do_Dynamics=true;
}
else{
    Do_Dynamics=false;
}

if(perform_rixs_=="true"){
    Perform_RIXS=true;
}
else{
    Perform_RIXS=false;
}

if(do_static_finite_temp_=="true"){
    Static_Finite_Temp=true;
}
else{
    Static_Finite_Temp=false;
}

if(do_dynamics_finite_temp_=="true"){
    Dynamics_Finite_Temp=true;
}
else{
    Dynamics_Finite_Temp=false;
}


if(restricted_basis_=="true"){
    Restricted_Basis=true;
}
else{
    Restricted_Basis=false;
}


}





