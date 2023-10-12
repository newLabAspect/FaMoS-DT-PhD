# Hybrid System Learning

To use all MATLAB scripts the parent directory (e.g. ".../hybrid-system-learning") needs to be selected as the working directory in MATLAB!

## HAutLearn

Provides the Matlab code to reconstruct the underlying Hybrid Automaton of captured/simulated traces using the HAutLearn Frame developed by Yang, Beg, Kenigsberg and Johnson. The reconstructed HA is outputted in multiple formats, the XML is the easiest to use. System Mode IDs and Cluster IDs are not interchangeable!

## Proposed Algorithm

Provides the Matlab code to reconstruct the underlying Hybrid Automaton of captured/simulated traces using the proposed algorithm which is based on HAutLearn but uses DTW for clustering and DTL for transition learning. The transitions are implicitly reconstructed with a DT, thus no fully reconstructed HA is outputted. This directory also includes some examples used for development purposes.

## Example Systems

All examplary systems used to evaluate the accuracy and performance of the proposed algorithm are provided in this folder. Each subfolder is dedicated to one example and includes a script to generate traces with different initial conditions and varying initial system modes. Multiple main programs are included which are used to carry out isolated test, comparisons with HAutLearn, and so on.

The following unique features are associated to each exemplary system:

| System                                  | Unique Feature                        |
|-----------------------------------------|---------------------------------------|
| Simple Heating System                   | Simple Reference System               |
| Multiple Variables Heating System       | Two Output Variables                  |
| Two State Hybrid Automaton              | ODEs of Degree Two                    |
| Three State Hybrid Automaton            | Different ODE Degrees in System Modes |
| Three State Timed Hybrid Automaton      | Timed Transitions                     |
| Three State Stochastic Hybrid Automaton | Stochastic Transitions                |
| Complex Hybrid Automaton                | Large Scale Example                   |

## Prerequisites

The following MATLAB toolboxes are required to use the framework (as of October 2023):

- Signal Processing Toolbox (various use-cases)
- Robust Control Toolbox (mainly used for LMIs)
- Statistics and Machine Learning Toolbox (mainly used for Decision Tree Learning)
- Symbolic Math Toolbox (mainly used for Automata Learning)
- Control System Toolbox (mainly used to convert discrete state space models to continous ones)

