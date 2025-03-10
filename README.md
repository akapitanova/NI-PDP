# Parallel Algorithm for Minimum Edge Cut  

## Author  
**Anna Kapitánová**  
Magisterské studium, FIT ČVUT  

## Project Overview  
This project implements parallel algorithms to solve the **Minimum Edge Cut** problem in a weighted, undirected graph.  
The solution explores **OpenMP Functional Parallelism**, **OpenMP Data Parallelism**, and **MPI** approaches, comparing their efficiency.  

📄 **For detailed analysis and results, see the full report:** [PDP.pdf](PDP.pdf)  

## Problem Definition  
Given a graph **G(V, E)** with **n** nodes and edge weights in **[80,120]**,  
find a partition of nodes into **X** and **Y** such that:  
- **|X| = a**, **|Y| = n - a**  
- The sum of edge weights between X and Y is **minimized**  

## Implementations  

### 1️⃣ **Sequential Algorithm**  
- **Branch and Bound (B&B) with Depth-First Search (DFS)**  
- Prunes non-optimal branches using **lower bound estimation**  

### 2️⃣ **OpenMP Parallelization**  
#### **Functional Parallelism**  
- Uses **omp parallel shared** to create tasks  
- Work is distributed via **omp task**  

#### **Data Parallelism**  
- **BFS** generates initial states, stored in a queue  
- **omp parallel for** distributes work dynamically  
- Uses **#pragma omp critical** to update the best solution safely  

### 3️⃣ **MPI Parallelization**  
- **Master-Slave model**  
- **Master process** distributes tasks  
- **Slave processes** compute solutions using OpenMP  
- **Communication:** `MPI_Send` and `MPI_Recv`  
- **Global minimum updates shared between processes**  

## How to Run  
### Sequential Execution  
```sh
./program graph_data.txt 5
