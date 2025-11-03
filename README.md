# Gene Editor Kinetics Calculator: Cas_Solver+

Programmable endonuclease-based editing strategies are rapidly being produced for gene therapy development, GMO production, biological research, and more. Studying the kinetics of an editor's DNA cleavage and subsequent cell repair activites provides essential information for the development and optimisation of gene editing strategies.

Cas_solver+ provides robust quantification of the underlying dynamics of gene editor activity across cell types and genomic targets. The script analyzes time-course editing data in parallel with and without potent DNA repair inhibitors (e.g., AZD7648 and ART558), allowing accurate estimation of fixed cutting/repair rate coefficients, recurrent cleavage, and precise repair frequencies.


## Publications

- [Unveiling the Cut-and-Repair Cycle of Designer Nucleases in Human Stem and T Cells via CLEAR-time dPCR (Nature Communications)](https://www.researchsquare.com/article/rs-4577114/v1)

## Applications
* Measure precise DNA repair and recurrent Cas9 cleavage
* Quantify the rates, probabilities, and frequencies of any experimentally measured DSB repair products
* Optimize the speed and safety of novel gene editing platforms
* Elucidate deeper biological mechanisms of DSB repair

## Schematic

<img width="600" align="center" alt="Overview" src="https://github.com/A-Chalk/DSB-Kinetics-Calculator/assets/167306438/c1fd58b2-a67c-4b29-aa92-d06c2bfdc6e4">

**a.** Programmable endonucleases like Cas9 recognise and cleave target sequences, producing a number of intended and unintended DSB repair products. **b.** If the cell manages to repair an induced DSB, it can either generate a reaction terminating mutation (no longer recognized by the endonuclease), or it can precisely restore the wildtype sequence, which can then be reidentified by Cas9 and cleaved again. This recurrent cleavage is what ultimately drives the generation of mutation products. **c.** By sampling from an edited cell population over a timeseries and processing the DNA with CLEAR-Time ddPCR safety screening, a model can be generated describing the activity of DSB generation and cellular repair.

## Pipeline Installation and Usage
* To use Cas_solver+ (1.7), download the latest script version shared in this repository and Cas9_Solver_Input.xlsx (includes dummy test data). Functional on R version 4.4.1 (_but older/newer versions may also function_).
* The following R packages require installation:
  - deSolve
  - tidyverse
  - readxl
  - foreach
  - openxlsx
  - doParallel
  
    
* Upon opening Cas_solver+ (1.7) in a program like Rstudio:
  - **Specify the location of the downloaded Cas9_Solver_Input.xlsx on your system (line 32)**
  - Modify the number of cores you would like to allocate to computation (line 19, default = 3 less than total number of detected cores)
  - Specify the number of bootstrap replicates you want to generate (line 24, default = num_bootstrap_reps <- 10, _recommend 1000 for final analysis_)
  - Whether or not to export results as a csv (export_results, line 25; default = TRUE)
  - Number of sample replicates you would like analysed (line 33, default = 3)

## Model Description

**Cas_solver+ (1.7)**: Solves for rate coefficient $k$ values based on experimental data and fits curves (_requires Cas9 Only and Cas9 + Repair inhibitor_ samples)

This pipeline employs a developed version of the three-state model of DSB repair in which Cas9 cleavage, precise DNA repair, and mutatagenic repair are each constrained by fixed rate coefficients $k_{dsb}(t)$, $k_{pr}$, and $k_{m}$ respectively.

$$Wildtype \xrightleftharpoons[k_{pr}]{k_{dsb}(t)}DSB\xrightarrow{k_m}Mutation$$

There are many possible reaction terminating products ($k_m$) detected by CLEAR-Time ddPCR, including NHEJ indels, large deletions, and targeted gene integration, each with distinct rate coefficients $k_{in}$, $k_{ld}$, $k_{ti}$. 

This model can be described continuously across the timeseries by a set of ordinary differential equations:

* Wildtype rate | $\frac{dWildtype}{dt} = (k_{pr} * DSB)-(k_{dsb} * D(t) * WT)$<br/>
* DSB rate | $\frac{dDSB}{dt} = (k_{dsb} * D(t) * WT) - (k_{pr} + k_{in} + k_{ld}) * DSB$<br/>
* Precise repair rate | $\frac{dPR}{dt} = (k_{pr} * DSB)-(k_{dsb} * D(t) * PR)$<br/>
* Indel rate | $\frac{dIndel}{dt} = k_{in} * DSB$<br/>
* Large deletion rate | $\frac{dLD}{dt} = k_{ld} * DSB$<br/>
* Targeted integration rate | $\frac{dTI}{dt} = k_{ti} * DSB$<br/>
* Cas9 trafficking delay | $D(t) = 1 - 2^{-(\frac{t}{delay})}$

The pipeline uses these ODEs to estimate the $k$ values based on observed mutation allele population frequencies. Notice that $k_{dsb}$ is multiplied by a delay function $D(t)$, which takes into account a brief initial lag in Cas9 cleavage activity resulting from gene editor trafficking into the nucleus post-electroporation (this variable can adjust for slower/faster Cas9 delivery strategies and quantify Cas9 searching times).

## Experimental Setup

Experimental data must be generated prior to analysis with a comprehensive mutation detection strategy like [CLEAR-Time ddPCR safety profiling]([https://github.com/A-Chalk/MEGA](https://www.nature.com/articles/s41467-025-65182-4)) applied in a timeseries on gene edited cells. Using such a technique, it is possible measure shifts in all mutation frequencies over time. Analysing kinetics timeseries data with Cas_solver reveals underlying dynamics of DSB generation and subsequent repair.

Accurate kinetics information can be generated by editing your cell population alongside repair inhibitors that block cell repair. Our approach involves nucleofecting cells with complexed Cas9 ribonucleoprotein (RNP), splitting the cells, and adding AZD7648 + ART558 in DMSO to half. Additionally, if you are interested in studying donor template integration, you can split each group again and transfect half with an AAV (etc).

**This will produce five experimental groups** to track over the timecourse:
1. Mock untreated cells  
2. RNP only  
3. RNP + inhibitors  
4. RNP + donor  
5. RNP + inhibitors + donor  


## Data Generation
Information on how quantify the absolute frequencies of each mutation within a gene-edited cell population can be found in the [CLEAR-Time ddPCR paper](https://www.researchsquare.com/article/rs-4577114/v1). This approach provides frequencies of indels (grouped together), DSBs, large deletions, targeted integration, and can also be modified to add translocations, MMEJ indels, and other chromosomal aberrations. 


## Authors
* Alexander Chalk
* Nathan White
* Giandomenico Turchiano

## License

This project is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html).
Commercial use is restricted under this license; please contact the author for inquiries.

## References
1. Brinkman, E. K. *et al.* (2018). Kinetics and Fidelity of the Repair of Cas9-Induced Double-Strand DNA Breaks. *Mol. Cell* **70**, 801–813.e6.
2. Ben-Tov, D. *et al.* (2023). Uncovering the Dynamics of Precise Repair at CRISPR/Cas9-induced DSBs. *bioRxiv*, doi:10.1101/2023.01.10.523377
3. Cucinotta, F. A. *et al.* (2008). Biochemical Kinetics Model of DSB Repair and Induction of γ-H2AX Foci. *Radiat. Res.* **169**, 214–222.
4. Rose, J. C. *et al.* (2017). Rapidly inducible Cas9 and DSB-ddPCR to probe editing kinetics. *Nat. Methods* **14**, 891–896.


