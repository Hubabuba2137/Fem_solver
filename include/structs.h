#pragma once

#include "../src/matrix.cpp"

namespace Fem{

    struct Robin_bc{
        float alfa;
        float t_ext;
    };
    struct Neumann_bc{
        float q_n;
    };

    struct Bc_type{
        Robin_bc robin;
        Neumann_bc neumann;
    };

    struct Node{
        float x;
        float y;

        bool bc;

        Node(){
        this-> x = 0; 
        this-> y = 0; 
        }

        Node(float x, float y){
            this->x = x;
            this->y = y;
            this->bc = false;
        }
    };

    struct Element{
        int id;
        int node_ids[3];

        Matrix H_local;
        Matrix H_bc;
        Matrix P;
        Matrix C;

        Element(int n1, int n2, int n3) : H_local(3, 3), H_bc(3,3), P(3,1), C(3,3) {
            this->node_ids[0] = n1;
            this->node_ids[1] = n2;
            this->node_ids[2] = n3;
        }
    };

    struct Ref_triangle{
        std::vector<float> dNdxi = {-1, 1, 0};
        std::vector<float> dNdeta = {-1, 0, 1};

        float N1(float xi, float eta){
            return 1-xi-eta;
        }
        float N2(float xi, float eta){
            return xi;
        }
        float N3(float xi, float eta){
            return eta;
        }
    };

    struct GlobalData{
        float total_time;
        float time_step;
        float conductivity;
        float alfa;
        float tot;
        float init_temperature;
        float density;
        float specific_heat;
        int node_number;
        int elem_number;

        GlobalData(){
            this->total_time=0;
            this->time_step=0;
            this->conductivity=0;
            this->alfa=0;
            this->tot=0;
            this->init_temperature=0;
            this->density=0;
            this->specific_heat=0;
            this->node_number=0;
            this->elem_number=0;
        }
    };
}