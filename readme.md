# Links

[Proposal](proposal.md)

[Milestone Report](milestone_report.md)

[Final Report](final_report.pdf)

# Running the Code
Adjusting the number of threads should be done within the source code with `omp_set_num_threads`.
### Contractions:
- Compile the code on the GHC machines:
  
    `cd contractions`

    `g++ -std=c++17 -fopenmp -ltbb -DTBB_INTERFACE_NEW=1 contraction_hierarchies.cpp -o contractions`

- Run the code:

    `./contractions`

### Delta Stepping:
- Compile the code on the GHC machines:

    `cd deltaStepping`
    
    `make`

- Run the code:

    `./planner`
