# FaMoS – Fast Model Learning for Hybrid Cyber-Physical Systems using Decision Trees

This repository provides the code for the paper

_Swantje Plambeck, Aaron Bracht, Nemanja Hranisavljevic, and Goerschwin Fey. Famos - fast model learning for hybrid cyber-physical systems using decision trees. In International Conference on Hybrid Systems: Computation and Control (HSCC), 2024_

Bibtex entry
```
@inproceedings{HSCC2024,
  author          = {Swantje Plambeck and Aaron Bracht and Nemanja Hranisavljevic and Goerschwin Fey},
  title           = {FaMoS - Fast Model Learning for Hybrid Cyber-Physical Systems using Decision Trees},
  booktitle       = {International Conference on Hybrid Systems: Computation and Control (HSCC)},
  year            = {2024},
}
```

## Prerequisites

The following MATLAB toolboxes (an of course MATLAB R2023b) are required to use the framework (as of October 2023):

- Signal Processing Toolbox (various use-cases)
- Robust Control Toolbox (mainly used for LMIs)
- Statistics and Machine Learning Toolbox (mainly used for Decision Tree Learning)
- Symbolic Math Toolbox (mainly used for Automata Learning)
- Control System Toolbox (mainly used to convert discrete state space models to continous ones)

## Verify Results

1. Set this directory (in which this README resides) as your working directory in MATLAB
2. Expand the folder "ExampleSystems" by pressing the "+"-icon next to it on the left hand side of MATLAB. (Do not double click as this would change the working directory)
3. Expand an example of interest in the folder "ExampleSystems" as done in step 2. The following examples were used for the evaluation - other examples might not run properly with the current version:
    | Name in Paper                | Folder                                  |
    |------------------------------|-----------------------------------------|
    | Simple Heating System        | ExampleSystems/SimpleHeatingSystem      |
    | Variable Heating System      | ExampleSystems/InputHeatingSystem       |
    | Complex Tank System          | ExampleSystems/ComplexTankSystem        |
    | Two State Hybrid Automaton   | ExampleSystems/TwoStateHybridAutomaton  |
    | Buck Converter               | ExampleSystems/BuckConverter            |
4. In the respective folder double click on the file "short_main_trace.m" which is needed to execute FaMoS or HAutLearn for this specific example. Again, other executable files in this folder may not properly run with the current version. One exception is "createTrace.m" which was used to generate the training data (training1.mat - training10.mat), in which you can take a look to better understand the underlying system.
5. All "short_main_trace.m" files are configured such that FaMoS will be used when you press "Run", i.e. clustering using DTW with possibly an optional LMI refinement and DTL for model learning. You can reconfigure the file in the following way:
- IMPORTANT NOTE: When executing the file for the first time, MATLAB will ask you if the responding folder should be added to the path. Please do so. Additionally, if you want to change to another example, remove the previous example folder from the path to avoid shadowing. To do so right click on the folder of the example system and select "Remove from path" in the menu. There are no subfolders in the example system folders, so it does not matter what you select in the next step.
- "methodCluster" (roughly line 10) changes the clustering method, values are provided in the comment.
- "methodTraining" (roughly line 11) changes the training method, values are provided in the comment. Note: Change training set size to provided value in paper to achieve the same results, as explained in the following bullet point.
- "variedMetricSteps" (roughly lines 40 - 46) allows you to try out many values for a parameter to find a good parameterization. In the standard configuration the varied parameter is training set size. To reconstruct our results for HAutLearn's training method you have to change line 45 to reflect the provided training set size. Feel free to vary other parameters, "tolLi" is probably the most interesting. You will find the corresponding value for "variedMetric" in the adjacent comment. Only the training paras can be varied.
- Here is an overview of how variable names are called in "short_main_trace.m" compared to the names in the paper if you want to dive in really deep:
    | Parameter Paper                  | Description                      | Variable Name                                                                             |
    |----------------------------------|----------------------------------|-------------------------------------------------------------------------------------------|
    | $T$                              | Sampling period                  | Ts                                                                                        |
    | $l$                              | Number of output variables       | num_var                                                                                   |
    | $k$                              | Number of input variables        | num_ud                                                                                    |
    | $w$                              | Window width                     | windowSize                                                                                |
    | $w_{size}$                       | Minimum distance of changepoints | windowSize                                                                                |
    | $d_{max}$                        | Maximum order of derivative      | max_deriv                                                                                 |
    | $d_{c}$                          | Order of derivative              | offsetCluster                                                                             |
    | $\alpha, \rho_{min}, \rho_{max}$ | Similarity thresholds            | facThres, thresClusterMin, thresClusterMax                                                |
    | $\sigma$                         | LMI tolerance                    | sigma                                                                                     |
    | $k_c$                            | Order of difference Equations    | offsetCluster                                                                             |
    |                                  | Use time for transition learning | useTime (DTL, false for all examples) Time (HAL, false for all examples, triggers errors) |
    | $\eta, \gamma, \lambda$          | PTA paras HAutLearn              | eta, gamma, lambda                                                                        |
    | $\psi$                           | Relexation para for LIs          | tolLi                                                                                     |
    |                                  | Legacy (unused) paras            | fixedIntervalLenght, precisionDTL                                                         |
- The conformance degree is saved in the variable "ConfDeg" after executing the program. "ClusterCorrect" and "ClusterFalse" are the results of the clustering evaluation. "t_cluster" and "t_train" are the measured times of the clustering and training time, respectively. Note that we did three back-to-back runs and selected the minimal duration oberserved for our time performance measurements provided in the table. "correct" and "false" are leftovers from an earlier evaluation format.
- Last note: If the clustering is too off, the traning will produce errors. Set a break point in "traceMain.m" (folder InterOperability, again open the folder by clicking on the "+" to not change the working directory) line 71 to see the clustering evaluation results.
6. If you want to see the prediction against the ground truth trace, copy the following in the command window after running the program: "hold on; plot(pred_trace.x(:,1)); plot(trace(evalData(1)).x(:,1));". Rerun the program with a different training method and rerun this command to see the other prediction plotted against the ground truth. Finally, replace the "1" right before each semicolon with an "2" to see the second output variable (if the system has at least two) plotted in this scenario.
7. **Happy Replication!** 
8. If you want to dive into the source code, read on further (not only this bullet point). A good starting point is "traceMain.m" (folder InterOperability, again open the folder by clicking on the "+" to not change the working directory). If you want to see a function definition, just put your cursor on it and press CTRL + D (you need to run an example first so that all folders are added to properly the path). Many of the specific functions of these two frameworks start with "Fn...". Set breakpoints along the execution path, run an example and view the workspace variables to better understand what is going on. I guess that is all. Simple, right? xD

## Repository Structure

In the following subsections the general layout of the files in this repository is explained.

### HAutLearn

Provides the Matlab code to reconstruct the underlying Hybrid Automaton of captured/simulated traces using the HAutLearn Frame developed by Yang, Beg, Kenigsberg and Johnson. The reconstructed HA is outputted in multiple formats, the XML is the easiest to use. System Mode IDs and Cluster IDs are not interchangeable!

### Proposed Algorithm

Provides the Matlab code to reconstruct the underlying Hybrid Automaton of captured/simulated traces using the proposed algorithm which is based on HAutLearn but uses DTW for clustering and DTL for transition learning. The transitions are implicitly reconstructed with a DT, thus no fully reconstructed HA is outputted.

### Example Systems

All examplary systems used to evaluate the accuracy and performance of the proposed algorithm are provided in this folder. Each subfolder is dedicated to one example and includes a script to generate traces with different initial conditions and varying initial system modes. Multiple main programs are included which are used to carry out isolated test, comparisons with HAutLearn, and so on.

The following unique features are associated to each exemplary system:

| System                                  | Unique Feature                        |
|-----------------------------------------|---------------------------------------|
| Buck Converter                          | Based on Real Physical Model          |
| Complex Hybrid Automaton                | Large Scale Example                   |
| Input Heating Sytsem                    | Input Variables used                  |
| Multiple Room Heating System            | Simplified Example of Benchmark       |
| Multiple Variables Heating System       | Two Output Variables                  |
| Simple Heating System                   | Simple Reference System               |
| Smart Factory                           | Real Data (very much WiP)             |
| Three State Hybrid Automaton            | Different ODE Degrees in System Modes |
| Three State Timed Hybrid Automaton      | Timed Transitions                     |
| Three State Stochastic Hybrid Automaton | Stochastic Transitions                |
| Two State Hybrid Automaton              | ODEs of Degree Two                    |


