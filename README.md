Design of high-performance parallel Graph interface supporting efficient Dynamic batch updates.

Research in graph-structured data has grown rapidly due to graphs' ability to represent complex real-world information and capture intricate relationships, particularly as many real-world graphs evolve dynamically through edge/vertex insertions and deletions. This has spurred interest in programming frameworks for managing, maintaining, and processing such dynamic graphs. In our report, we evaluate the performance of [PetGraph (Rust)], [Stanford Network Analysis Platform (SNAP)], [SuiteSparse:GraphBLAS], [cuGraph], [Aspen], and [our custom implementation] in tasks including loading graphs from disk to memory, cloning loaded graphs, applying in-place edge deletions/insertions, and performing a simple iterative graph traversal algorithm. Our implementation demonstrates significant performance improvements: it outperforms PetGraph, SNAP, SuiteSparse:GraphBLAS, cuGraph, and Aspen by factors of `177x`, `106x`, `76x`, `17x`, and `3.3x` in graph loading; `20x`, `235x`, `0.24x`, `1.3x`, and `0x` in graph cloning; `141x`/`45x`, `44x`/`25x`, `13x`/`11x`, `28x`/`34x`, and `3.5x`/`2.2x` in edge deletions/insertions; and `67x`/`63x`, `86x`/`86x`, `2.5x`/`2.6x`, `0.25x`/`0.24x`, and `1.3x`/`1.3x` in traversal on updated graphs with deletions/insertions.

Below, we plot the runtime (in seconds, logarithmic scale) for loading a graph from file into memory with PetGraph, SNAP, SuiteSparse:GraphBLAS, cuGraph, Aspen, and **Our DiGraph** for each graph in the dataset.

![Image](https://github.com/user-attachments/assets/3ef09c02-3fa1-4558-9cef-45531093a47a)

Next, we plot the runtime (in milliseconds, logarithmic scale) of deleting a batch of `10^−7|𝐸|` to `0.1|𝐸|` randomly generated edges into a graph, **in-place**, in multiples of `10`. Here, we evaluate PetGraph, SNAP, SuiteSparse:GraphBLAS, cuGraph, Aspen, and Our DiGraph on each graph in the dataset. The left subfigure presents overall runtimes using the geometric mean for consistent scaling, while the right subfigure shows runtimes for individual graphs.

![Image](https://github.com/user-attachments/assets/7252a1fc-abb2-4f38-afd5-0479e08c6468)


Below, we plot the runtime of inserting a batch of edges into a graph, **in-place**, using PetGraph, SNAP, SuiteSparse:GraphBLAS, cuGraph, Aspen, and **Our DiGraph**.

![Image](https://github.com/user-attachments/assets/a5edbd2e-012e-4d52-9edc-fc5f83bf2d21)

Finally, we plot the runtime of traversing a graph using a simple iterative algorithm (`42`-step reverse walks from each vertex in a graph) on graphs with edge deletions. We evaluate PetGraph, SNAP, SuiteSparse:GraphBLAS, cuGraph, Aspen, and **Our DiGraph** on each graph in the dataset.

![Image](https://github.com/user-attachments/assets/60ec8b22-e4d0-4ee7-8752-21bda26add01)

Refer to our technical report for more details: \
[Performance Comparison of Graph Representations Which Support Dynamic Graph Updates][report].

<br>

> [!NOTE]
> You can just copy `main.sh` to your system and run it. \
> For the code, refer to `main.cxx`.

<br>
<br>


### Code structure

The code structure is as follows:

```bash
- inc/_*.hxx: Utility functions
- inc/**.hxx: Common graph functions
- inc/_memory.hxx: The memory allocators (FAA, AA, CP2AA)
- inc/csr.hxx: Compressed Sparse Row (CSR) data structure functions
- inc/Graph.hxx: Graph data structure
- inc/io.hxx: Graph file reader (MTX format)
- inc/main.hxx: Main header
- main.cxx: Experimentation code
- main.sh: Experimentation script
- process.js: Node.js script for processing output logs
```

Note that each branch in this repository contains code for a specific experiment. The `main` branch contains code for the final experiment. If the intention of a branch in unclear, or if you have comments on our technical report, feel free to open an issue.

<br>
<br>


## References

- [Algorithm 1037: SuiteSparse:GraphBLAS: Parallel Graph Algorithms in the Language of Sparse Linear Algebra; Timothy A. Davis et al. (2023)](https://dl.acm.org/doi/full/10.1145/3577195)
- [Low-latency graph streaming using compressed purely-functional trees; Laxman Dhulipala et al. (2019)](https://dl.acm.org/doi/abs/10.1145/3314221.3314598)
- [cuGraph C++ primitives: vertex/edge-centric building blocks for parallel graph computing; Seunghwa Kang et al. (2023)](https://ieeexplore.ieee.org/abstract/document/10196665)
- [SNAP: A General-Purpose Network Analysis and Graph-Mining Library; Jure Leskovec et al. (2016)](https://dl.acm.org/doi/abs/10.1145/2898361)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)
- [How can I convert a std::string to int?](https://stackoverflow.com/a/7664227/1413259)
- [Fastest way to read numerical values from text file in C++ (double in this case)](https://stackoverflow.com/a/5678975/1413259)
- [What's the difference between istringstream, ostringstream and stringstream? / Why not use stringstream in every case?](https://stackoverflow.com/a/3292168/1413259)
- [c++ stringstream is too slow, how to speed up?](https://stackoverflow.com/a/5830907/1413259)
- [Best Approach to read huge files utilizing multithreading; Stephan van Hulst :: Coderanch](https://coderanch.com/t/699934/java/Approach-read-huge-files-utilizing)
- [How to get current time and date in C++?](https://stackoverflow.com/a/997988/1413259)
- [Signed variant of size_t in standard C++ library](https://stackoverflow.com/q/65496071/1413259)
- [Is 'signed size_t' different from 'ssize_t'?](https://stackoverflow.com/q/20744349/1413259)
- [How to create a temporary directory?](https://stackoverflow.com/a/4632032/1413259)
- [How to amend a commit without changing commit message (reusing the previous one)?](https://stackoverflow.com/a/10365442/1413259)
- [Syntax for a single-line while loop in Bash](https://stackoverflow.com/a/1289029/1413259)
- [How can I save username and password in Git?](https://stackoverflow.com/a/35942890/1413259)
- [How do I tell git to use fewer cores/threads when compressing?](https://superuser.com/a/539478/305990)
- [Containers library :: cppreference](https://en.cppreference.com/w/cpp/container)
- [Date and time utilities :: cppreference](https://en.cppreference.com/w/cpp/chrono)
- [Standard library header &lt;string&gt; :: cppreference](https://en.cppreference.com/w/cpp/header/string)
- [Standard library header &lt;algorithm&gt; :: cppreference](https://en.cppreference.com/w/cpp/header/algorithm)

<br>
<br>


[![](https://img.youtube.com/vi/yqO7wVBTuLw/maxresdefault.jpg)](https://www.youtube.com/watch?v=yqO7wVBTuLw)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)


[PetGraph (Rust)]: https://github.com/petgraph/petgraph
[Stanford Network Analysis Platform (SNAP)]: https://github.com/snap-stanford/snap
[SuiteSparse:GraphBLAS]: https://github.com/GraphBLAS/LAGraph
[cuGraph]: https://github.com/rapidsai/cugraph
[Aspen]: https://github.com/ldhulipala/aspen
[our custom implementation]: https://github.com/puzzlef/graph-openmp
[sheets-o1]: https://docs.google.com/spreadsheets/d/102WZCbN0cGFns8VlCoY_b-_dh-5C9JhPKgUg2d32WU0/edit?usp=sharing
[report]: https://arxiv.org/abs/2502.13862
