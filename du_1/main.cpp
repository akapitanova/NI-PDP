#include <iostream>
#include <vector>
#include <string>
#include <map>
#include<set>
#include <fstream>
#include <sstream>
#include <limits>


class Solution{
private:
    int price;
    int n_vertices;
    // TODO prepsat na vector
    std::map<int, int> X_vertices;

public:
    Solution(int price, int n_vertices, std::map<int, int> X_vertices){
        this->price = price;
        this->n_vertices = n_vertices;
        this->X_vertices = X_vertices;
    }

    Solution(int n){
        price = 0;
        n_vertices = 0;
        for (int i = 0; i < n; i++){
            X_vertices.insert({i, 0});
        }
    }

    int get_price(){
        return price;
    }

    int get_n_vertices(){
        return n_vertices;
    }

    std::map<int, int> get_X_vertices(){
        return X_vertices;
    }

};

class Graph{
public:
    // celkem pocet vrcholu
    int n;
    // prumerny stupen vrcholu
    int k;
    // pocet vrcholu mnoziny X
    int a;

    std::vector<std::vector<int>> edges_w;
    Solution* best_sol;
    int best_price;

    int calls;

    Graph(int n, int k, int a, std::vector<std::vector<int>> edges_w){
        this->n = n;
        this->k = k;
        this->a = a;

        this->edges_w = edges_w;

        best_sol = nullptr;
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

    int get_edge_w(int u_i, int v_i){
        return edges_w[u_i][v_i];
    }

    Solution* get_best_sol(){
        return this->best_sol;
    }

    int get_best_price(){
        return best_price;
    }

    int get_calls(){
        return calls;
    }

    // BB-DFS pro hledani minimalniho hranoveho rezu
    void solve(Solution* solution){
        calls += 1;
        int n_vertices = solution->get_n_vertices();

        // ukončíme větev, pokud velikost menší podmnožiny uzlů překročí a
        if (n_vertices > a){
            return;
        }

        if (n_vertices == a){
            // nalezeno prvni validni reseni
            if(!best_sol){
                best_sol = solution;
                best_price = solution->get_price();
            }

            // cena menší než cena dosavadního nejlepšího řešení
            // provede se aktualizace nejlepšího řešení
            else if (this->smaller_than_min(solution)){
                best_sol = solution;
                best_price = solution->get_price();
            }
            return;
        }

        // volani rekurze
        // TODO vector - jako pointer
        std::map<int, int> X_vertices = solution->get_X_vertices();
        std::set<int> free_vertices = get_free_vertices(X_vertices);

        // TODO iterovat pres vector, continue pokud je obsah 1
        for (auto it = free_vertices.begin(); it != free_vertices.end(); ++it){
            // TODO podminka if obsah 1 -> continue
            // TODO az zde delat deep copy
            std::map<int, int> X_vertices_i = X_vertices;
            X_vertices_i[*it] = 1;

            int price_i = get_weights_sum(X_vertices_i);
            if (price_i == -1){
                continue;
            }

            int n_vertices_i = n_vertices + 1;

            Solution solution_i = Solution(price_i, n_vertices_i, X_vertices_i);

            // Větev ukončíme, i pokud v daném mezistavu pomocí dolního odhadu váhy zbývajícího řezu zjistíme,
            // že z něj nelze žádným způsobem vytvořit přípustný koncový stav s menší vahou řezu,
            // než je nejlepší dosud nalezené řešení
            //if (!check_lower_cost_estimation(&solution_i)){
            //    continue;
            //}
            solve(&solution_i);
        }
    }


private:

    // return vertices in graph not already in X set (value in map is = 0)
    std::set<int> get_free_vertices(std::map<int, int> X_vertices){
        std::set<int> free_vertices;
        for (int i = 0; i < this->n; i++){
            if (X_vertices[i] == 0){
                free_vertices.insert(i);
            }
        }
        return free_vertices;
    }

    //TODO reference
    int get_weights_sum(std::map<int, int> X_vertices){
        int sum = 0;

        for (int row = 0; row < n; row++){
            for (int column = 0; column < n; column++){
                if (X_vertices[row] == 1 and X_vertices[column] == 0){
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

    bool check_lower_cost_estimation(Solution* solution){
        if(!best_sol) return true;

        int price = solution->get_price();
        std::map<int, int> X_vertices = solution->get_X_vertices();
        std::set<int> free_vertices = get_free_vertices(X_vertices);
        int sum = price;
        //TODO zase zbytecny set
        for (auto it = free_vertices.begin(); it != free_vertices.end(); ++it){
            std::map<int, int> X_vertices_i = X_vertices;
            X_vertices_i[*it] = 1;
            int new_price = get_weights_sum(X_vertices_i);

            sum += new_price - price;
        }

        if (sum > best_price){
            return false;
        }
        return true;
    }


    bool smaller_than_min(Solution* solution){
        if (solution->get_price() < best_price){
            return true;
        }
        return false;
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


std::ostream & operator<< (std::ostream& os, Solution* solution)
{
    std::map<int, int> X_vertices = solution->get_X_vertices();
    for (auto i : X_vertices){
        if (i.second == 1){
            os << i.first << ", ";
        }
    }
    os << std::endl;

    return os;
}

int main() {

    std::string dir_path = "/home/anna/skola/PDP/ni_pdp/graf_mro/";


    std::string file_name = "graf_20_7.txt";
    Graph* graph = load_data(dir_path, file_name, 7);
    Solution first_solution = Solution(graph->get_n());
    graph->solve(&first_solution);
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
