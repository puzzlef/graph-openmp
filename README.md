[OpenMP]-based parallel graph implementation.

I have been trying to parallelize the graph data structure from bottom up.
- Each vertex has a **list of edges**, which is a **sorted vector of pairs** of *edge id* and *weight*.
- A sorted list has better locality, lookup, and update time.
- To amortize cost of edge deletions and insertions, they are *simply queued up* and lazily updated.
- *Either* a batch of deletions *or* insertions must be done at a time to use lazy behaviour (not a mix).
- The lazy update is done manually by calling `update()` with a temporary per-thread buffer.
- The buffer is needed to do an *in-place* `set_union_last_unique()` / `set_difference_unique()`.
- Different threads must work on *different* vertices.

<!-- [![](https://i.imgur.com/Jp1UDS5.png)][sheetp] -->
<!-- [![](https://i.imgur.com/7lA6tWb.png)][sheetp] -->
<!-- [![](https://i.imgur.com/170NBzh.png)][sheetp] -->
<!-- [![](https://i.imgur.com/rdyR5Uo.png)][sheetp] -->

<br>

The last time i tried parallelizing graph reading with 12 threads, performance
as 1/2 that of sequential, with converting each line to numbers taking the most
time. It seems that the default stream operator (`>>`) uses locks for some
reason to support locales, and is not recommended to use in parallel code.

[![](https://i.imgur.com/WgPl7nr.png)][sheetp]

<br>

I switched to using `strtoull()` and `strtod()` for the conversion, switched to
`istringstream` from `stringstream`, and parallelized edge addition into the
graph (along with parallel update). On `web-Google` graph we are now getting a
speedup of **~6x** for graph loading. I now log system date and time for
detecting any more bottlenecks. The parallel edge update is as follows. Only one
of the 12 threads will add the edge to the graph if the source vertex belongs to
its chunk.

<!-- [![](https://i.imgur.com/ONf0uPi.png)][sheetp] -->
[![](https://i.imgur.com/EBCyi5u.png)][sheetp]
<!-- [![](https://i.imgur.com/lB8xouh.png)][sheetp] -->
<!-- [![](https://i.imgur.com/pMRwmIa.png)][sheetp] -->
<!-- [![](https://i.imgur.com/VJMyves.png)][sheetp] -->

<br>

Graph functions such as `symmetricize()` are now simply parallelized as below.

[![](https://i.imgur.com/tJNoNYO.png)][sheetp]

<br>

We are able to load graphs with more that *4 billion edges* in under **10**
**minutes**. Most of the time spent is loading the graph is the time needed to
*read the file*, and *add edges to the graph data structure*. We get a net
loading rate of about **10 million edges per second**. *Updating the graph*,
which involves sorting the edges of each vertex while picking only the last
unique entry, is about *50 times faster*.

[![](https://i.imgur.com/9FIflvU.png)][sheetp]
[![](https://i.imgur.com/PSALt0j.png)][sheetp]
[![](https://i.imgur.com/bkqzHLa.png)][sheetp]

<br>

All outputs are saved in a [gist] and a small part of the output is listed here.
Some [charts] are also included below, generated from [sheets]. The input data
used for this experiment is available from the [SuiteSparse Matrix Collection].
This experiment was done with guidance from [Prof. Kishore Kothapalli] and
[Prof. Dip Sankar Banerjee].


[OpenMP]: https://www.openmp.org

<br>

```bash
$ g++ -std=c++17 -O3 main.cxx
$ ./a.out ~/Data/GAP-road.mtx
$ ./a.out ~/Data/GAP-twitter.mtx
$ ...

# 2022-12-17 21:04:36 Loading graph /home/subhajit/Data/GAP-road.mtx ...
# 2022-12-17 21:04:36 OMP_NUM_THREADS=24
# 2022-12-17 21:04:42 readMtxOmpW(): vertices=269.8ms, read=3173.7ms, parse=218.2ms, edges=2290.6ms, update=193.3ms
# 2022-12-17 21:04:42 order: 23947347 size: 57708624 [directed] {}
# 2022-12-17 21:04:42 [06160.472 ms] readMtxOmpW
#
# 2022-12-17 21:04:43 Loading graph /home/subhajit/Data/GAP-twitter.mtx ...
# 2022-12-17 21:04:43 OMP_NUM_THREADS=24
# 2022-12-17 21:08:07 readMtxOmpW(): vertices=711.1ms, read=138848.7ms, parse=10906.0ms, edges=50429.3ms, update=3494.1ms
# 2022-12-17 21:08:07 order: 61578415 size: 1468364884 [directed] {}
# 2022-12-17 21:08:07 [204464.219 ms] readMtxOmpW
#
# ...
```

<br>
<br>


## References

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
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)

<br>
<br>

[![](https://i.imgur.com/LWJP5Hy.jpg)](https://www.youtube.com/watch?v=iHsqqgvwUxk)<br>
[![ORG](https://img.shields.io/badge/org-ionicf-green?logo=Org)](https://ionicf.github.io)
[![DOI](https://zenodo.org/badge/576609870.svg)](https://zenodo.org/badge/latestdoi/576609870)


[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[gist]: https://gist.github.com/wolfram77/1b2c4f07a1dc46a330b1dc14afa2b4ab
[charts]: https://imgur.com/a/94omat4
[sheets]: https://docs.google.com/spreadsheets/d/16M2A3ucmqjSr1JL-WWnT-_h3Edggsk_-MvT2MALfHPs/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vToyQvdBJF-sc9QZ8X2cL6udirwEhWmQnusLT6HgtxYkdRnwrcoJoMDpwC0RMh1Dzh5a4cdrmkDMlRg/pubhtml
