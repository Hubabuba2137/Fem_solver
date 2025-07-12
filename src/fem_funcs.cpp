#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

#include "../include/structs.h"
#include "../include/matrix.h" 
#include <algorithm>

std::vector<double> pc_xi = {0.1667, 0.6667, 0.1667};
std::vector<double> pc_eta = {0.1667, 0.1667, 0.6667};
std::vector<double> weights_pc = {0.1667, 0.1667, 0.1667};

std::vector<double> n1_bc_val = {0.8873, 0.5000, 0.1127, 0, 0, 0, 0.1127, 0.5000, 0.8873 };
std::vector<double> n2_bc_val = {0.1127, 0.5000, 0.8873, 0.8873, 0.5000, 0.1127, 0, 0, 0};
std::vector<double> n3_bc_val = {0, 0, 0, 0.1127, 0.5000, 0.8873, 0.8873, 0.5000, 0.1127,};


std::vector<double> weights_bc = {0.2777, 0.4444, 0.2777};

float dist(Fem::Node a, Fem::Node b){
    return sqrt(pow( a.x - b.x,2)+pow(a.y - b.y,2));
}


Fem::Matrix calc_jacobian_mat(Fem::Element &element, std::vector<Fem::Node> &nodes){
    Fem::Matrix jacobian(2,2);

    Fem::Node n1 = nodes[element.node_ids[0]-1];
    Fem::Node n2 = nodes[element.node_ids[1]-1];
    Fem::Node n3 = nodes[element.node_ids[2]-1];

    //      |x2-x1  x3-x1|
    // J =  |            |
    //      |y2-y1  y3-y1|

    jacobian[0][0] = n2.x-n1.x;
    jacobian[0][1] = n3.x-n1.x;
    jacobian[1][0] = n2.y-n1.y;
    jacobian[1][1] = n3.y-n1.y;

    return jacobian;
}

float calc_jacobian(Fem::Matrix &J){
    return (J[0][0]*J[1][1])-(J[1][0]*J[0][1]);
}

Fem::Matrix inverse_jacobian_matrix(Fem::Matrix &J){
    float inv_det_J = 1/calc_jacobian(J);

    Fem::Matrix inv_jac(2,2);

    inv_jac[0][0] = inv_det_J*J[1][1];
    inv_jac[1][1] = inv_det_J*J[0][0];
    inv_jac[1][0] = -inv_det_J*J[0][1];
    inv_jac[0][1] = -inv_det_J*J[1][0];

    return inv_jac;
}

Fem::Matrix calc_local_H(Fem::Element &local_el, Fem::Ref_triangle &ref_el, std::vector<Fem::Node> &nodes, float conductivity)
{
    Fem::Matrix J = calc_jacobian_mat(local_el, nodes);
    double detJ = calc_jacobian(J);
    Fem::Matrix invJ = inverse_jacobian_matrix(J);

    Fem::Matrix dNdx(3,1), dNdy(3,1);
    for(int i=0; i<3; ++i){
        dNdx[i][0] = invJ[0][0]*ref_el.dNdxi[i] + invJ[0][1]*ref_el.dNdeta[i];
        dNdy[i][0] = invJ[1][0]*ref_el.dNdxi[i] + invJ[1][1]*ref_el.dNdeta[i];
    }

    Fem::Matrix H_local = (dNdx*dNdx.transpose()) + (dNdy*dNdy.transpose());

    double area = 0.5 * std::fabs(detJ);
    double scale = conductivity * area;

    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            H_local[i][j] *= scale;

    return H_local;
}

Fem::Matrix calc_local_Hbc(Fem::Element &local_el, std::vector<Fem::Node> &nodes, float alfa){
    Fem::Matrix hbc(3,3);

    for(int i=0; i<3; i++){//pętla po bokach
        //n1 - i
        //n2 - (i+1)%3

        //odejmujemy tutaj 1 od id żeby mieć przeniesienia
        //z indeksowania o bazie 1, które jest w elemencie
        //na indeksowanie z bazą 0, które jest w vector nodes
        int id1 = local_el.node_ids[i]-1;
        int id2 = local_el.node_ids[(i+1)%3]-1;
        
        if(!nodes[id1].bc && !nodes[id2].bc){
            continue;
        }

        //std::cout<<"\t"<<id1<<", "<<id2<<"\n"; //sprawdzenie jakie węzły są brane pod uwagę w czasie całkowania
        //sprawdzenie czy np. nody bez warunków brzegowych nie są brane pod uwagę
        float det_J = 0.5*dist(nodes[id1], nodes[id2]);
        Fem::Matrix hbc_i(3,3);

        for(int j=0; j<3; j++){//pętla po pc
            int k = i*3;
            //std::cout<<"k = "<<k<<"\n";
            //std::cout<<"k+j = "<<k+j<<"\n";
            Fem::Matrix h_pc(3,1);
            h_pc[0][0] = n1_bc_val[j+k];
            h_pc[1][0] = n2_bc_val[j+k];
            h_pc[2][0] = n3_bc_val[j+k];

            //std::cout<<h_pc[0][0]<<", "<<h_pc[1][0]<<", "<<h_pc[2][0]<<"\n";
            Fem::Matrix h_pc_temp(3,3);

            h_pc_temp = h_pc*h_pc.transpose();

            //std::cout<<h_pc_temp;
            
            for(int a=0; a<3; a++){
                for(int b=0; b<3; b++){
                    h_pc_temp[a][b] = h_pc_temp[a][b]*weights_bc[j];
                }
            }
            //std::cout<<"h_pc *"<<weights_bc[j]<<"=\n"<<h_pc_temp<<"\n";

            hbc_i = hbc_i + h_pc_temp;
            //std::cout<<hbc_i<<"\n";
        }
        
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                hbc_i[i][j] = hbc_i[i][j]*det_J*alfa;
            }
        }

        //std::cout<<hbc_i<<"\n";
        hbc = hbc+hbc_i;
    }

    return hbc;
}

Fem::Matrix calc_p_vec(Fem::Element &local_el, std::vector<Fem::Node> &nodes, float alfa, float T_ext){
    Fem::Matrix p_vec(3,1);

    for(int i =0; i<3; i++){
        int id1 = local_el.node_ids[i]-1;
        int id2 = local_el.node_ids[(i+1)%3]-1;
        
        if(!nodes[id1].bc && !nodes[id2].bc){
            continue;
        }

        float det_J = 0.5*dist(nodes[id1], nodes[id2]);
        Fem::Matrix p_i(3,1);

        for(int j=0; j<3;j++){//pętla po pc na jednym boku
            int k = i*3;//upewniamy się że przechodzimy przez wszystkie 9 boków

            Fem::Matrix p_temp(3,1);
            p_temp[0][0] = n1_bc_val[j+k];
            p_temp[1][0] = n2_bc_val[j+k];
            p_temp[2][0] = n3_bc_val[j+k];

            for(int i=0; i<3; i++){
                p_temp[i][0] = p_temp[i][0]*T_ext*weights_bc[j];
            }
            //std::cout<<p_temp<<"\n";

            p_i = p_i + p_temp;
        }
        
        for(int i=0; i<3; i++){
            p_i[i][0] = p_i[i][0]*alfa*det_J;
        }
        //std::cout<<p_i<<"\n";
        
        p_vec = p_vec + p_i;
    }

    return p_vec;
}

Fem::Matrix calc_c(Fem::Element &local_el, Fem::Ref_triangle ref_el, std::vector<Fem::Node> &nodes, float density, float conductivity){
    Fem::Matrix c_mat(3,3);
    Fem::Matrix el_jac = calc_jacobian_mat(local_el, nodes);
    float det_J = calc_jacobian(el_jac);
    //std::cout<<det_J<<"\n";

    for(int i=0; i<3; i++){
        Fem::Matrix c_i(3,1);
        c_i[0][0] = ref_el.N1(pc_xi[i], pc_eta[i]);
        c_i[1][0] = ref_el.N2(pc_xi[i], pc_eta[i]);
        c_i[2][0] = ref_el.N3(pc_xi[i], pc_eta[i]);

        //std::cout<<c_i<<"\n";
        Fem::Matrix c_temp = c_i*c_i.transpose();
        //std::cout<<"Before\n"<<c_temp<<"\n";

        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){
                c_temp[j][k] *= weights_pc[i]*std::fabs(det_J)*density*conductivity;
            }
        }
        //std::cout<<"After\n"<<c_temp<<"\n";

        c_mat = c_mat + c_temp;
        //std::cout<<c_mat<<"\n";
    }

    return c_mat;
}

void aggregate(Fem::Matrix &Global, Fem::Element element, Fem::Matrix &Local){
    for(int i=0; i<3;i++){
        for(int j=0; j<3;j++){
            int glob_i = element.node_ids[i]-1;
            int glob_j = element.node_ids[j]-1;
            Global[glob_i][glob_j] += Local[i][j]; 
        }
    }
}

void aggregate_p_vec(Fem::Matrix &P_vec, Fem::Element element, Fem::Matrix &Local){
    for(int i=0; i<3; i++){
        int glob_i = element.node_ids[i]-1;
        P_vec[glob_i][0] += Local[i][0];

    }
}