#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <set>

struct Node {
    int id;
    double x, y;
};

struct Triangle {
    int id;
    int node1, node2, node3;
};

struct Bar {
    int id;
    int node1, node2;
};

void nas_to_mesh(std::string file_in, std::string file_out) {
    std::ifstream input_file(file_in);
    std::ofstream output_file(file_out);
    
    if (!input_file.is_open()) {
        std::cerr << "Error: Cannot open input file " << file_in << std::endl;
        return;
    }
    
    if (!output_file.is_open()) {
        std::cerr << "Error: Cannot open output file " << file_out << std::endl;
        return;
    }
    
    std::vector<Node> nodes;
    std::vector<Triangle> triangles;
    std::vector<Bar> bars;
    std::set<int> boundary_nodes;
    
    std::string line;
    
    // Parse NAS file
    while (std::getline(input_file, line)) {
        std::istringstream iss(line);
        std::string keyword;
        iss >> keyword;
        
        if (keyword == "GRID") {
            Node node;
            
            // NASTRAN uses fixed-width fields
            // GRID format: GRID, ID, CP, X1, X2, X3
            if (line.length() >= 48) {
                std::string id_str = line.substr(8, 8);
                std::string cp_str = line.substr(16, 8);
                std::string x_str = line.substr(24, 8);
                std::string y_str = line.substr(32, 8);
                
                node.id = std::stoi(id_str);
                node.x = std::stod(x_str);
                node.y = std::stod(y_str);
            } else {
                // Fallback for free format
                std::istringstream iss_fallback(line);
                std::string temp;
                int coord_sys;
                iss_fallback >> temp >> node.id >> coord_sys >> node.x >> node.y;
            }
            
            nodes.push_back(node);
        }
        else if (keyword == "CTRIA3") {
            Triangle triangle;
            int prop_id;
            iss >> triangle.id >> prop_id >> triangle.node1 >> triangle.node2 >> triangle.node3;
            triangles.push_back(triangle);
        }
        else if (keyword == "CBAR") {
            Bar bar;
            int prop_id;
            double orient1, orient2, orient3;
            iss >> bar.id >> prop_id >> bar.node1 >> bar.node2 >> orient1 >> orient2 >> orient3;
            bars.push_back(bar);
            
            // Add nodes from bars to boundary nodes set
            boundary_nodes.insert(bar.node1);
            boundary_nodes.insert(bar.node2);
        }
        else if (keyword == "ENDDATA") {
            break;
        }
    }
    
    input_file.close();
    
    // Write output file
    
    // Write nodes section
    output_file << "*Nodes" << std::endl;
    for (const auto& node : nodes) {
        output_file << node.id << ", " << node.x << ", " << node.y << std::endl;
    }
    
    // Write elements section (renumber starting from 1)
    output_file << "*Elements" << std::endl;
    int element_counter = 1;
    for (const auto& triangle : triangles) {
        output_file << element_counter << ", " << triangle.node1 << ", " 
                    << triangle.node2 << ", " << triangle.node3 << std::endl;
        element_counter++;
    }
    
    // Write boundary conditions section
    output_file << "*BC" << std::endl;
    for (const auto& node_id : boundary_nodes) {
        output_file << node_id << std::endl;
    }
    
    output_file.close();
}