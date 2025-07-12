#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>

#include "../include/structs.h"
#include "../include/matrix.h" 

#include "load_function.cpp"
#include "fem_funcs.cpp"
#include "gauss.cpp"

void solve(std::string file_name, bool print_conf = 1, bool print_H = 0, bool print_C = 0, bool print_P=0){
    std::string file = "../" + file_name;

    std::vector<Fem::Node> nodes = load_nodes(file);
    std::vector<Fem::Element> elements = load_elements(file);
    Fem::GlobalData configuration = load_configuration(file);
    load_bc(file, nodes);

    if(print_conf){
        print_config(configuration);
    }

    Fem::Ref_triangle ref_el;

    Fem::Matrix H_glob(configuration.node_number, configuration.node_number);
    Fem::Matrix C_glob(configuration.node_number, configuration.node_number);
    Fem::Matrix P_glob(configuration.node_number, 1);
    
    for(int i =0; i<configuration.elem_number;i++){ //pętla po elementach
        elements[i].H_local = calc_local_H(elements[i], ref_el, nodes, configuration.conductivity); // macierz H lokalna
        elements[i].H_bc = calc_local_Hbc(elements[i], nodes, configuration.alfa); //lokalna macierz Hbc
        elements[i].H_local = elements[i].H_local+elements[i].H_bc; //sumowanie lokalnej H i Hbc
        elements[i].P = calc_p_vec(elements[i], nodes, configuration.alfa, configuration.tot);
        elements[i].C = calc_c(elements[i], ref_el, nodes, configuration.density, configuration.specific_heat);

        //aggregate(H_glob, elements[i], elements[i].H_local);
        //aggregate(C_glob, elements[i], elements[i].C);
        //aggregate_p_vec(P_glob, elements[i], elements[i].P);
    }
    std::cout<<"Local matrices and vectors constructed\n";
    
    for(int i =0; i<configuration.elem_number;i++){
        aggregate(H_glob, elements[i], elements[i].H_local);
        aggregate(C_glob, elements[i], elements[i].C);
        aggregate_p_vec(P_glob, elements[i], elements[i].P);
    }

    std::cout<<"Agregation finished\n";

    if(print_H){
        std::cout<<"Global H:\n"<<H_glob<<"\n";
    }

    if(print_C){
        std::cout<<"Global C:\n"<<C_glob<<"\n";
    }

    if(print_P){
        std::cout<<"Global P:\n"<<P_glob<<"\n";
    }

    Fem::Matrix Global(configuration.node_number, configuration.node_number);

    std::vector<double> t0(configuration.node_number);
    std::vector<double> t(configuration.node_number);

    for(int i=0; i<configuration.node_number;i++){
        t0[i] = configuration.init_temperature;
    }

    for(int i=0; i<configuration.node_number; i++){
        for(int j=0; j<configuration.node_number; j++){
            Global[i][j] = H_glob[i][j] + (C_glob[i][j] / configuration.time_step);
        }
    }
    
    //std::cout << "MIN: " << *std::min_element(t0.begin(), t0.end()) << " MAX: " << *std::max_element(t0.begin(), t0.end()) << std::endl;

    for(float i=configuration.time_step; i<configuration.total_time; i+=configuration.time_step){
        for(int j=0; j<configuration.node_number; j++){
            double rhs = 0;
            for(int k = 0; k<configuration.node_number; k++){
                rhs += (C_glob[j][k] / configuration.time_step) * t0[k];
            }
            rhs += P_glob[j][0];
            t[j] = rhs;
        }

        t = Gauss(Global, t);
        std::cout<<"Temperature at time " << i << "s:\n";
        
        /*
        for(size_t a=0; a<t.size(); a++){
            std::cout<<"Temperature in node "<<a<<" = "<<t[a]<<", ";
        }*/

        std::cout << "\nMIN: " << *std::min_element(t.begin(), t.end()) << " MAX: " << *std::max_element(t.begin(), t.end()) << std::endl;

        t0=t;
    }
}
