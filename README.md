# HERE-release

Implementation of the paper:

[Efficient and Robust Point Cloud Registration via Heuristics-guided Parameter Search] (T-PAMI 2024).

To apply this algorithm to your data, you may need to adjust the noise_bound (Line-42 in ./demo.cpp) and BRANCH_ACCURACY (Line-22 in ./include/utils.h). For example, smaller BRANCH_ACCURACY leads to higher accuracy but requires more time cost.

## Requirements
- CMake: 3.22.1
- Eigen: 3.4.0
- Boost: 1.65.1
- PCL: 1.9.1
- OpenMP: 201511

## Citation

```bibtex
@article{Huang_2024_TPAMI,
    author       = {Huang, Tianyu and Li, Haoang and Peng, Liangzu and Liu, Yinlong and Liu, Yun-Hui},
    title        = {Efficient and Robust Point Cloud Registration via Heuristics-guided Parameter Search},
    journal      = {IEEE Transactions on Pattern Analysis and Machine Intelligence},
    year         = {2024}
    publisher={IEEE}
}
```

## Our another work regarding 3D registration
[Scalable 3D Registration via Truncated Entry-wise Absolute Residuals](https://github.com/tyhuang98/TEAR-release) (CVPR 2024)


## Acknowledgements
The Branch-and-Bound algorithm part of our code is built on [GORE](https://cs.adelaide.edu.au/~aparra/project/gore/). We thank the authors for their excellent work!
