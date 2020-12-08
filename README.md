A Variable Strategy Reinforced Lin-Kernighan-Helsgaun Algorithm (VSR-LKH) 
----
This repository contains code to use our algorithm to solve the TSP and reproduce results from the paper: <br> <br>
Combining Reinforcement Learning with Lin-Kernighan-Helsgaun Algorithm for the Traveling Salesman Problem (AAAI 2021) <br>
Jiongzhi Zheng, Kun He, Jianrong Zhou, Yan Jin, Chu-min Li <br> <br>

Installation
----
On a Unix/Linux machine execute the following commands: <br> <br>

cd VSR-LKH <br>
make <br> <br>

An executable file called LKH will now be available in the directory VSR-LKH. <br>
Then enter the .par file name corresponding to the TSP instance, such as [u574.par](./u574.par), to run the program. <br> <br>

File Description
----
VSR-LKH was achieved on top of the famous TSP heuristic, Lin-Kernighan-Helsgaun (LKH) algorithm. You can learn LKH from its open source website, http://akira.ruc.dk/~keld/research/LKH/, and understand the effect of each file and the meaning of the parameters from the PDF files in the directory [DOC](./DOC). The description of files that differ between VSR-LKH and LKH is as follows: <br> <br>

* The statement of reinforcement learning parameters and Q-value: [LKH.h](./SRC/INCLUDE/LKH,h) <br>
* The initialization of Q-values. The candidate sets is sorted according to Q-values: [GenerateCandidates.c](./SRC/GenerateCandidates.c), [AdjustCandidateSet.c](./SRC/AdjustCandidateSet.c), [ResetCandidateSet.c](./SRC/ResetCandidateSet.c). <br>
* The reinforcement learning process in k-opt: [BestKOptMove.c](./SRC/BestKOptMove.c) (for parameters PATHING_A = 2, PATHING_C = 3), [Best5OptMove.c](./SRC/Best5OptMove.c) (for parameters PATHING_A = 1, PATHING_C = 0). <br>
* The initialization of parameters: [ReadParameters.c](./SRC/ReadParameters.c) <br>
* The reinforcement learning parameters adaptation process: [FindTour.c](./SRC/FindTour.c) <br> <br>

Data
----
We give three TSP instances (u574, u1060, rl11849) as examples. All the TSPLIB Format instance (.tsp) can be calculated by VSR-LKH. The all 111 TSPLIB instances can be found in the directory [TSPLIB_DATA](./TSPLIB_DATA) or http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/. <br> <br>

Contact
----
Questions and suggestions can be sent to jzzheng@hust.edu.cn. <br> <br>

Citation
----
If you find this code and data useful, please consider citing the original work by authors: <br>
```
@article{zheng2021RL-TSP,
  title={Combining Reinforcement Learning with Lin-Kernighan-Helsgaun Algorithm for the Traveling Salesman Problem},
  author={Jiongzhi Zheng and Kun He and Jianrong Zhou and Yan Jin and Chu-min Li},
  journal={AAAI Conference on Artificial Intelligence},
  year={2021}
}
```
