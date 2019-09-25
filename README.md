# Fault_Diagnosis
Identification of causal models for Fault Propagation Analysis on chemical processes.

Because the variables in industrial processes are closely related to each other, the 
occurrence of a disturbance usually causes its spread through the plant. Causal models 
are a useful tool for fault propagation analysis, as well as to locate the original fault
that caused the disturbance. To be able to form this model it is necessary to find the 
cause-effect relationships between the process variables. Variables with greater causal 
impact on the others have the highest probability of being considered as the root cause. 
In this work, several data management methods are proposed and analyzed to extract the 
causality from historical data of the process. As a case of test a continuously stirred 
heated tank was used.

## Getting Started

This project was completely developed in C# and implement 3 mathematical methods used
for detect causality between time series. Also implement the clases to create a 
qualitative model of the founded relationships using Signed Directed Graph (SDG).

### Prerequisites

Visual Studio Comunity 2015.

### Installing

Run the "DiagnosticoDeFallas.sln" file.

## Authors

* **Gabriel Pérez García** - *Initial work* - [ggarciaperez9](https://github.com/ggarciaperez9)
