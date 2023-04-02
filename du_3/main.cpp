#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <queue>
#include <fstream>
#include <sstream>
#include <limits>
#include <omp.h>



std::ostream & operator<< (std::ostream& os, std::vector<int> solution)
{
    int n = solution.size();
    for (int i = 0; i < n; i++){
        os << solution[i] << " ";
    }
    os << std::endl;

    os << "X vertices: ";
    for (int i = 0; i < n; i++){
        if (solution[i]){
            os << i << " ";
        }
    }
    os << std::endl;

    os << "Y vertices: ";
    for (int i = 0; i < n; i++){
        if (!solution[i]){
            os << i << " ";
        }
    }
    return os;
}


class Graph{
    // celkem pocet vrcholu
    int n;
    // prumerny stupen vrcholu
    int k;
    // pocet vrcholu mnoziny X
    int a;
    std::vector<std::vector<int>> edges_w;

public:
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

    std::vector<std::vector<int>> get_edges_w(){
        return edges_w;
    }
};

class State{
    int vertex;
    int price;
    int ones;
    int zeros;
    std::vector<int> solution;

public:
    State(int vertex, int price, int ones, int zeros, std::vector<int> solution){
        this->vertex = vertex;
        this->price = price;
        this->ones = ones;
        this->zeros = zeros;
        this->solution = solution;
    }

    int get_vertex(){
        return vertex;
    }

    int get_price(){
        return price;
    }

    int get_ones(){
        return ones;
    }

    int get_zeros(){
        return zeros;
    }

    std::vector<int> get_solution(){
        return solution;
    }

    void change_solution(int index, int value, Graph graph){
        solution[index] = value;
        price = get_weights_sum(graph);
    }

    int get_weights_sum(Graph graph){
        int sum = 0;
        std::vector<std::vector<int>> edges_w = graph.get_edges_w();

        for (int row = 0; row < graph.get_n(); row++){
            for (int column = row + 1; column < graph.get_n(); column++){
                if ((solution[row] == 1 and solution[column] == 0)
                    or (solution[row] == 0 and solution[column] == 1)){
                    sum += edges_w[row][column];
                }
            }
        }
        return sum;
    }



};

void print_final_weights(Graph graph, std::vector<int> solution){
    std::vector<std::vector<int>> edges_w = graph.get_edges_w();

    for (int row = 0; row < graph.get_n(); row++){
        for (int column = row + 1; column < graph.get_n(); column++){
            if ((solution[row] == 1 and solution[column] == 0)
                or (solution[row] == 0 and solution[column] == 1)){
                if (edges_w[row][column] != 0) {
                    std::cout << "(" << row << ", " << column << ") " << edges_w[row][column] << std::endl;
                }
            }
        }
    }
}

class MinEdgeCut{
    std::vector<int> best_solution;
    int min_price;
    int calls;

public:
    MinEdgeCut(){
        this->min_price = std::numeric_limits<int>::max();
        this->calls = 0;
    }

    std::vector<int> get_best_solution(){
        return best_solution;
    }

    int get_min_price(){
        return min_price;
    }

    int get_calls(){
        return calls;
    }

    bool check_lower_cost_estimation(Graph graph, State state, int& min_price){
        int price = state.get_price();
        int sum = price;
        int vertex = state.get_vertex();
        std::vector<int> solution = state.get_solution();
        int n = graph.get_n();

        for (int i = vertex; i < n; i++){
            solution[i] = 0;
            int new_price_0 = state.get_weights_sum(graph);

            solution[i] = 1;
            int new_price_1 = state.get_weights_sum(graph);

            sum += std::min(new_price_0 - price, new_price_1 - price);

            solution[i] = -1;
        }

        if (sum > min_price){
            return false;
        }
        return true;
    }


    void bfs(Graph& graph, std::queue<State> * q){

        State state = q->front();
        q->pop();

        int vertex = state.get_vertex();
        int ones = state.get_ones();
        int zeros = state.get_zeros();
        int n = graph.get_n();
        int a = graph.get_a();
        int price = state.get_price();
        std::vector<int> solution = state.get_solution();

        // kontrola, ze se jedna o legitimni stavy
        if (zeros > n - a){
            return;
        }
        if (ones > a){
            return;
        }

        State state_0 = State(vertex + 1, price, ones, zeros + 1, solution);
        state_0.change_solution(vertex, 0, graph);

        State state_1 = State(vertex + 1, price, ones + 1, zeros, solution);
        state_1.change_solution(vertex, 1, graph);

        // pridani stavu do fronty
        q->push(state_0);
        q->push(state_1);
    }


    // BB-DFS pro hledani minimalniho hranoveho rezu
    void solve_seq(Graph& graph, State& current, int& calls, int& min_price, std::vector<int>& best_solution){

        #pragma omp atomic
        calls += 1;

        int vertex = current.get_vertex();
        int n = graph.get_n();
        int a = graph.get_a();

        if (current.get_zeros() > n - a){
            return;
        }
        if (current.get_ones() > a){
            return;
        }

        if (current.get_price() > min_price){
            return;
        }

        if (vertex == n){
            if (current.get_ones() == a && current.get_price() < min_price){
                #pragma omp critical
                {
                    if (current.get_ones() == a && current.get_price() < min_price){
                        min_price = current.get_price();
                        best_solution = current.get_solution();
                    }
                }
            }
            return;
        }

        State state_0 = State(vertex + 1, current.get_price(), current.get_ones(), current.get_zeros() + 1, current.get_solution());
        state_0.change_solution(vertex, 0, graph);

        State state_1 = State(vertex + 1, current.get_price(), current.get_ones() + 1, current.get_zeros(), current.get_solution());
        state_1.change_solution(vertex, 1, graph);

        if (check_lower_cost_estimation(graph, state_0, min_price)){
            solve_seq(graph, state_0, calls, min_price, best_solution);
        }

        if (check_lower_cost_estimation(graph, state_1, min_price)){
            solve_seq(graph, state_1, calls, min_price, best_solution);
        }

    }

    std::vector<int> init_solution(int n){
        std::vector<int> solution;
        for (int i = 0; i < n; i++){
            solution.push_back(-1);
        }
        return solution;
    }

    void make_min_edge_cut(Graph& graph, int threads, int init_states_number){
        std::vector<int> solution = init_solution(graph.get_n());
        State root = State(0, 0, 0, 0, solution);

        std::queue<State> * q = new std::queue<State>();
        q->push(root);

        // generovani pocatecnich stavu pomoci bfs
        while (true) {
            if (int(q->size()) >= init_states_number) break;

            if (q->front().get_vertex() == graph.get_n() - 1) break;

            bfs(graph, q);
        }

        std::vector<State> states;
        while (!q->empty()){
            states.push_back(q->front());
            q->pop();
        }

        delete q;

        omp_set_num_threads(threads);

        #pragma omp parallel for schedule(dynamic) shared (states, min_price, best_solution, calls, graph)
            for (auto & state: states){
                solve_seq(graph, state, calls, min_price, best_solution);
            }

    }


    void print_results(Graph& graph, std::string file_name){
        std::cout << "file: " << file_name << std::endl;
        std::cout << "best price: " << get_min_price() << std::endl;
        std::cout << "total calls: " << get_calls() << std::endl;
        std::cout << "solution: " << get_best_solution() << std::endl;
        print_final_weights(graph, best_solution);
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




int main(int argc, char * argv[]) {

    if (argc < 3){
        return 1;
    }

    std::string dir_path = "../graf_mro/";
    std::string file_name(argv[1]);
    int a = std::stoi(argv[2]);
    int threads = std::stoi(argv[3]);
    int init_states_number = std::stoi(argv[4]);

    Graph* graph = load_data(dir_path, file_name, a);
    if (!graph){
        return 1;
    }

    MinEdgeCut min_edge_cut = MinEdgeCut();

    min_edge_cut.make_min_edge_cut(*graph, threads, init_states_number);

    min_edge_cut.print_results(*graph, file_name);

}
