#include <iostream>
#include <vector>
#include <string>
#include <map>
#include<set>
#include <fstream>
#include <sstream>
#include <limits>

class Graph{
public:
    // celkem pocet vrcholu
    int n;
    // prumerny stupen vrcholu
    int k;
    // pocet vrcholu mnoziny X
    int a;

    std::vector<std::vector<int>> edges_w;
    std::vector<int> best_sol;
    int best_price;

    int calls;

    Graph(int n, int k, int a, std::vector<std::vector<int>> edges_w){
        this->n = n;
        this->k = k;
        this->a = a;

        this->edges_w = edges_w;

        best_price = std::numeric_limits<int>::max();
        calls = 0;
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

    std::vector<int> get_best_sol(){
        return best_sol;
    }

    int get_best_price(){
        return best_price;
    }

    int get_calls(){
        return calls;
    }

    void make_min_edge_cut(){
        std::vector<int> X_vertices = init_X_vertices();
        solve(0, X_vertices, 0, 0);
    }

private:
    std::vector<int> init_X_vertices(){
        std::vector<int> X_vertices;
        for (int i = 0; i < n; i++){
            X_vertices.push_back(-1);
        }
        return X_vertices;
    }

    int get_weights_sum(std::vector<int>& X_vertices){
        int sum = 0;

        for (int row = 0; row < n; row++){
            for (int column = row + 1; column < n; column++){
                if ((X_vertices[row] == 1 and X_vertices[column] == 0)
                or (X_vertices[row] == 0 and X_vertices[column] == 1)){
                    // Průběžně počítáme váhu řezu pro již přiřazené uzly,
                    // pokud je tato váha větší než dosud nalezené minimum,
                    // tuto větev ukončíme
                    sum += edges_w[row][column];
                    if (sum > best_price){
                        return -1;
                    }
                }
            }
        }

        return sum;
    }

    bool check_lower_cost_estimation(int vertex, std::vector<int>& solution, int price){
        int sum = price;
        //std::vector<int> solution_i = solution;
        for (int i = vertex; i < n; i++){
            solution[i] = 0;
            int new_price_0 = get_weights_sum(solution);

            solution[i] = 1;
            int new_price_1 = get_weights_sum(solution);

            sum += std::min(new_price_0 - price, new_price_1 - price);

            solution[i] = -1;
        }

        if (sum > best_price){
            return false;
        }
        return true;
    }

    // BB-DFS pro hledani minimalniho hranoveho rezu
    void solve(int vertex, std::vector<int>& X_vertices, int ones, int price){
        calls += 1;

        if (ones > a){
            return;
        }

        if (price == -1){
            return;
        }

        if (!check_lower_cost_estimation(vertex, X_vertices, price)){
            return;
        }

        if (vertex == n){
            if (ones == a && price < best_price){
                best_sol = X_vertices;
                best_price = price;
            }
            return;
        }

        // volani rekurze

        // left son
        std::vector<int> X_vertices_left = X_vertices;
        X_vertices_left[vertex] = 0;
        int price_left = get_weights_sum(X_vertices_left);

        solve(vertex + 1, X_vertices_left, ones, price_left);

        //int vertex_right = vertex + 1;
        std::vector<int> X_vertices_right = X_vertices;
        X_vertices_right[vertex] = 1;
        int price_right = get_weights_sum(X_vertices_right);

        solve(vertex + 1, X_vertices_right, ones + 1, price_right);
    }

};

Graph* load_data(std::string& dir_path, std::string& file_name, int a){
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


std::ostream & operator<< (std::ostream& os, std::vector<int> solution)
{
    int n = solution.size();
    for (int i = 0; i < n; i++){
        os << solution[i] << " ";
    }

    return os;
}



int main() {

    std::string dir_path = "/home/anna/skola/PDP/ni_pdp/graf_mro/";


    std::string file_name = "graf_30_10.txt";
    Graph* graph = load_data(dir_path, file_name, 10);

    graph->make_min_edge_cut();
    std::cout << "solution " << graph->get_best_sol() << std::endl;
    std::cout << "best price " << graph->get_best_price() << std::endl;
    std::cout << "total calls " << graph->get_calls() << std::endl;


    /*
    std::string file_name_1 = "graf_20_7.txt";
    Graph* graph_1 = load_data(dir_path, file_name_1, 10);
    Solution first_solution_1 = Solution(graph_1->get_n());
    graph_1->solve(&first_solution_1);
    //std::cout << graph->get_best_sol() << std::endl;
    std::cout << "best price " << graph_1->get_best_price() << std::endl;
    std::cout << "total calls " << graph_1->get_calls() << std::endl;
    return 0;
     */
}
