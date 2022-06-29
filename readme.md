# Mathematical Models of Protein Complex Formation Systems (PMID: 35729235)

Hsu, KL., Yen, HC.S. & Yeang, CH. Cooperative stability renders protein complex formation more robust and controllable. Sci Rep 12, 10490 (2022).

Protein complexes play crucial roles in a wide range of biological processes. Although numerous prior studies have identified several mechanisms involved in complex formation, few studies have directly assessed how these mechanisms benefit complex formation in terms of robustness and controllability. Understanding those properties would provide insights into how natural protein complexes are regulated and inform circuit design principles for artificial protein complexes. Here, we focus on one such mechanism, i.e., cooperative stability, whereby subunits are stabilized upon assembly into protein complexes. Through simulations of simple models of protein complex formation, we show that a protein complex subjected to cooperative stability is simultaneously robust against parameter configurations and sensitive to changes in synthesis of an individual subunit. Moreover, our *in silico* experiment allowed us to evaluate comprehensively the performance of four complex formation circuits in terms of eight system characteristics. Our work reveals the advantages and unique role of cooperative stability and opens up a new avenue for designing and controlling protein complex formation. 


## Installation

We recommend that you use Python 3.7 and Jupyter Notebook to run the simulations.  Please also install libraries listed in the requirements.txt 

- Install Python and Jupyter Notebook at Anaconda (https://www.anaconda.com/ )

- Using pip or conda to install/upgrade the required libraries

  ```bash
  $ python -m pip install --upgrade pip
  $ pip install -r requirements.txt
  ```

## Implementation

### Jupyter Notebook

- Once you installed the Anaconda, you can launch the Jupyter Notebook via Anaconda-Navigator or Anaconda prompt:

  ```bash
  $ Jupyter Notebook
  ```

- Open the Main.ipynb and start to play with the models interactively.

### Python Interactive mode

- Alternatively, you can directly use the Python interactive mode without Jupyter Notebook. Use command line and change directory to the ./mathematical_modeling_of_complex_formation/

- Enter the Python interactive mode:

  ```bash
  $ python 
  ```

- Import the functions from src/ and explore the models interactively.

  ```bash
  >>> from src.Robustness import draw_dimer_explore
  >>> draw_dimer_explore(lam1=0.027*5, lam3=0.027, lam12_equal=True)
  Normalized dimer (nM):	 2509.4782457253486
  Upper limit: 5.83469387755102
  Lower limit: 0.19421768707482992
  Tolerance score: 4.9089104699894115
  >>>
  ```

## Documentation 

We provide functions to explore the robustness and controllability of protein complex formation systems and evaluate/compare the performance of given circuits in eight system's characteristics. All the source codes are available in src/.

Model parameters: 

- lam1, lam2 and lam3 are the degradation rate constants ($hour^{-1}$) of p1 monomer, p2 monomer and dimer.
- K is the association constant ($nM^{-1}$) of p1 and p2.
- c1 and c2 are the synthesis rates ($nM/hour$)  of p1 and p2 monomers.

### Robustness

- draw_dimer_explore

  ```python
  from src.Robustness import draw_dimer_explore
  draw_dimer_explore(lam1=0.135, lam2=0.027, lam3=0.027, K=0.05, c2=98, fold=2, lam12_equal=True, normalize=True)
  ```

- tolerance_score_heatmap

  ```python
  from src.Robustness import tolerance_score_heatmap
  tolerance_score_heatmap(lam1=0.27, lam3=0.027,c_min=1,c_max=1000,k_min=0.001,k_max=1,resolution1 = 500, resolution2=500)
  ```

- tolerance_score_heatmap3D

  ```python
  from src.Robustness import tolerance_score_heatmap3D
  tolerance_score_heatmap3D(lam_min=0.027, lam_max=0.7, k_min=0.001, k_max=1, c2=10, offset=0, resolution=100)
  ```

### Controllability

- draw_dimer_response

  ```python
  from src.Controllability import draw_dimer_response
  draw_dimer_response(lam1=0.135, lam2=0.027, lam3=0.027, K=0.01, c2=98, fold=2, lam12_equal=True, normalize=True, compare=True)
  ```

- draw_trimer_response

  ```python
  from src.Controllability import scenario, draw_trimer_response
  draw_trimer_response([scenario.stable, scenario.three], xmiddle=98, who='c3')
  ```

- controllability_heatmap 

  ```python
  from src.Controllability import controllability_heatmap
  controllability_heatmap(lam1=0.27, lam3=0.027,c_min=1,c_max=1000,k_min=0.001,k_max=1,resolution1 = 500, resolution2=500, up=True)
  ```

### Performance evaluation 

- import the models and solvers

  ```python
  from src.Models_and_solvers import open_circuit, closed_circuit
  ```

- performance

  ```python
  from src.Performance_evaluation import performance
  
  open_circuit_performance = performance(open_circuit.dynamics)
  closed_circuit_performance = performance(closed_circuit.dynamics)
  ```

- evaluation_ss

  ```python
  from src.Performance_evaluation import evaluation_ss
  
  l = '_efficiency_n'
  o_wo, o_with, c_wo, c_with = evaluation_ss(open_circuit_performance.efficiency, closed_circuit_performance.efficiency, 
                 n=10, fc=1, l=l)
  ```

- evaluation_d

  ```python
  from src.Performance_evaluation import evaluation_d
  
  l = '_recover_n'
  o_wo, o_with, c_wo, c_with = evaluation_d(open_circuit_performance.recovery_time, closed_circuit_performance.recovery_time, 
                 n=10, shot=5, l=l)
  ```

- save_data

  ```python
  from src.Performance_evaluation import save_data
  
  result_ls = [o_wo, o_with, c_wo, c_with]
  for i in range(4):
      save_data('data/'+variables.circuits[i]+l, result_ls[i])
  ```

- read_data

  ```python
  from src.Performance_evaluation import read_data
  
  result_dict = {}
  for i in variables.circuits:
      for j in variables.suffix_n:
          key = i + '_' + j
          filename = 'data/' + key + '.txt'
          result = read_data(filename)
          result_dict[key] = result
  ```

### Performace comparison 

- variables 

  ```python
  from src.Multiple_comparison import variables
  from itertools import permutations as perm 
  
  circuit_name = ['owo', 'owith', 'cwo','cwith']
  z = perm(circuit_name, 4)
  
  ranking_ls = []
  for i in z:
      ranking_ls.append(i)
      
  global_variable = variables(result_dict, ranking_ls)
  ```

- scores_for_radar

  ```python
  from src.Multiple_comparison import scores_for_radar
  
  for i in global_variable.suffix_n:
      globals()[i+'ls']=scores_for_radar(i, global_variable, normalize=True)
  ```

  







