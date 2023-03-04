#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <sstream>


class Solution{
public:
    int price;
    int n_vertices;
    std::map<int, int> X_vertices;

    Solution(int price, int n_vertices, std::map<int, int> X_vertices){
        this->price = price;
        this->n_vertices = n_vertices;
        this->X_vertices = X_vertices;
    }

    Solution(Solution* old_solution){
        price = old_solution->price;
        n_vertices = old_solution->n_vertices;
        X_vertices = old_solution->X_vertices;
    }

    Solution(Solution* old_solution, int new_vertex){
        price = old_solution->price;
        n_vertices = old_solution->n_vertices;
        X_vertices = old_solution->X_vertices;
        this->add_vertex(new_vertex);
    }

    int get_price(){
        return price;
    }

    int get_n_vertices(){
        return n_vertices;
    }

    bool check_cond(int a){
        if (a == n_vertices){
            return true;
        }
        return false;
    }

    std::map<int, int> get_X_vertices(){
        return X_vertices;
    }
private:
    void add_vertex(int vertex){
        X_vertices[vertex] = 1;
    }
};

class Graph{
public:
    int n;
    int k;
    int a;
    std::vector<std::vector<int>> edges_w;
    Solution* bestSol;

    Graph(int n, int k, int a, std::vector<std::vector<int>> edges_w){
        this->n = n;
        this->k = k;
        this->a = a;

        this->edges_w = edges_w;
    }

    int get_n(){
        return n;
    }

    int get_a(){
        return a;
    }

    int get_k(){
        return k;
    }

    int get_edge_w(int u_i, int v_i){
        return edges_w[u_i][v_i];
    }

    void solve(Solution* solution){
        // Průběžně počítáme váhu řezu pro již přiřazené uzly,
        // pokud je tato váha větší než dosud nalezené minimum,
        // tuto větev ukončíme.
        if (this->bigger_than_min(solution)){
            return;
        }

        if (this->smaller_than_min(solution)){
            bestSol = solution;
        }

        int 

    }


private:

    int get_weights_sum(std::map<int, int> X_vertices){
        int sum = 0;

        for (int row = 0; row < n; row++){
            for (int column = 0; column < n; column++){
                if (X_vertices[row] == 1 and X_vertices[column] == 0){
                    sum += edges_w[row][column];
                }
            }
        }
    }

    bool bigger_than_min(Solution* solution){
        if (solution->get_price() > bestSol->get_price()){
            return true;
        }
        return false;
    }

    bool smaller_than_min(Solution* solution){
        if (solution->get_price() < bestSol->get_price()){
            return true;
        }
        return false;
    }


};

Graph* load_data(std::string & dir_path, std::string & file_name, int a){
    if (a < 5)
        return nullptr;

    std::stringstream ss(file_name.substr(0, file_name.size() - 4));
    std::string val;
    std::vector<std::string> tokens;
    while(std::getline(ss, val, '_')){
        tokens.push_back(val);
    }

    int n = std::stoi(tokens[1]);
    int k = std::stoi(tokens[2]);

    if (n >= 150 or n < 10)
        return nullptr;
    if (k >= n * (3.0/4) or k < 5)
        return nullptr;
    if (a > (n/2))
        return nullptr;

    std::vector<std::vector<int>> edges_w;
    std::cout << dir_path + file_name << std::endl;
    std::ifstream data(dir_path + file_name, std::ios::in);
    if (!data.is_open()) {
        return nullptr;
    }

    int tmp;
    data >> tmp;
    if (tmp != n)
        return nullptr;

    for (int row = 0; row < n; row++){
        std::vector<int> edges_w_row;
        for (int column = 0; column < n; column++){
            int edge_w = 0;
            data >> edge_w;
            edges_w_row.push_back(edge_w);
        }
        edges_w.push_back(edges_w_row);
    }

    Graph * graph = new Graph(n, k, a, edges_w);
    return graph;
}



int main() {

    std::string dir_path = "/home/anna/skola/PDP/ni_pdp/graf_mro/";
    std::string file_name = "graf_10_5.txt";
    Graph* graph = load_data(dir_path, file_name, 5);
    std::cout << graph->get_edge_w(9, 2);

    return 0;
}
