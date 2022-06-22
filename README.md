## About The Project

This is the tool for automatic dominance breaking for MiniZinc modelling language. 

You can use MiniZinc to model constraint satisfaction and optimisation problems
in a high-level, solver-independent way. The tool will analyze the model and generate 
dominannce breaking constraints if there are any. 

<!-- GETTING STARTED -->
## Getting Started

To get a MiniZinc up and running follow these simple steps.

### Prerequisites

- [CMake](https://cmake.org/) (>=3.4)
- A recent C++ compiler - Compilation is tested with recent versions of Clang,
  GCC, and Microsoft Visual C++.

Use the following command to clone the repo and the submodules
```
git clone --recurse-submodules https://github.com/AllenZzw/auto-dom.git
```
Then, enter the folder ```auto-dom``` and use the following command to compile the tool. 
```
bash build.sh
```

The package has been tested on Linux and Mac environment. 

### Usage

Once the tool is built on your machine, you can find a binary named "auto-dom" in the "build" directory. 
You can start using the tool to analyze your discrete optimisation models.  
The following code segment shows a MiniZinc model for the well known Knapsack problem.

```minizinc
int: n = 10; % number of items
int: W = 40; % weight limit
array [1..n] of int: w = [5, 7, 8, 3, 12, 9, 5, 6, 7, 4]; % weight of each item
array [1..n] of int: v = [3, 11, 5, 9, 11, 8, 13, 11, 2, 15];  % value of each item

array [1..n] of var 0..1: x;   % take item or not

constraint sum (i in 1..n) (w[i]*x[i]) <= W;

solve maximize sum (i in 1..n) (v[i]*x[i]);
```

Using the following command 
```
build/auto-dom knapsack.mzn 
```

The output will be the dominance breaking constraints for the knapsack problem like 
```
constraint x[1] != 1 \/ x[4] != 0;
constraint x[1] != 1 \/ x[7] != 0;
constraint x[1] != 0 \/ x[9] != 1;
constraint x[1] != 1 \/ x[10] != 0;
constraint x[2] != 0 \/ x[3] != 1;
constraint x[2] != 0 \/ x[6] != 1;
constraint x[2] != 1 \/ x[7] != 0;
constraint x[2] != 1 \/ x[8] != 0;
constraint x[2] != 0 \/ x[9] != 1;
constraint x[2] != 1 \/ x[10] != 0;
...
```

More details can be found in the follwing paper:
```
@inproceedings{jlee2022exploiting,
	title={Exploiting Functional Constraints in Generating Dominance Breaking Nogoods for Constraint Optimization},
	author={Jimmy H.M. Lee and Allen Z. Zhong},
	booktitle={Proceedings of the 28th International Conference on Principles and Practice of Constraint Programming},
	year={2022},
  	organization={Schloss Dagstuhl-Leibniz-Zentrum f{\"u}r Informatik}
}
```
