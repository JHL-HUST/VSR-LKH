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
Then enter the .par file name corresponding to the TSP instance, such as ./u574.par, to run the program. <br> <br>

File Description
----
VSR-LKH was achieved on top of the famous TSP heuristic, Lin-Kernighan-Helsgaun (LKH) algorithm. You can learn LKH from its open source website, http://akira.ruc.dk/~keld/research/LKH/, and understand the effect of each file and the meaning of the parameters from the PDF files in the directory ./DOC/. The description of files that differ between VSR-LKH and LKH is as follows: <br> <br>

* The initialization of Q-values. The candidate sets is sorted according to Q-values: ./SRC/GenerateCandidates.c, ./SRC/AdjustCandidateSet.c, ./SRC/ResetCandidateSet.c. <br>
* The reinforcement learning process in k-opt: ./SRC/BestKOptMove.c (for parameters PATHING_A = 2, PATHING_C = 3), ./SRC/Best5OptMove.c (for parameters PATHING_A = 1, PATHING_C = 0). <br>
* The initialization of reinforcement learning parameters: ./SRC/LKHmain.c <br>
* The reinforcement learning parameters adaptation process: ./SRC/FindTour.c <br> <br>

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
