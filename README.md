# HistChecker
This repository contains a filter program for checking properties involving Homeomorphically Irreducible Spanning Trees (HISTs). It is used for the article "J. Goedgebeur, K. Noguchi, J. Renders and C.T. Zamfirescu, HIST-Critical Graphs and Malkevitch’s Conjecture, manuscript".

The latest version of this program can be obtained from <https://github.com/JarneRenders/HistChecker>.

This program can be used the verify whether or not a graph is HIST-free, HIST-critical or simply to count how many HISTs are present.

### Installation

This requires a working shell and `make`. On Windows an easy way to simulate this is by using Windows Subsystem for Linux (WSL).

- Compile using: 
  * `make` to create a binary for the 64 bit version
  * `make 128bit` to create a binary for the 128 bit version
  * `make 128bitarray` to create a binary for the 128 bit version which implements bitsets using arrays
  * `make all` to create all of the above

The 64 bit version is much faster than the 128 bit version, hence it is recommended to use this one when running the program on graphs with 64 vertices or fewer.
There is a difference in implementation between both of the 128 bit versions, but it is not entirely clear which one is fastest in which situations. Both can handle input graphs up to order 128.
Use `make clean` to remove all binaries created using `make`.

### Usage of HistChecker

This helptext can be found by executing `./histChecker -h`.

Usage: `./histChecker [-1|-f] [-c] [-hv]`

Program for checking HISTs in a graph.

Graph are read from stdin in graph6 format. Graphs are sent to stdout in graph6 format. For more information on the format, see <http://users.cecs.anu.edu.au/~bdm/data/formats.txt>.

Without any parameters the program outputs graphs which are HIST-free to stdout.
```
  -1  : outputs graphs which are HIST-critical; cannot be used with -f
  -f  : count the number of HISTs in each graph and give information on the
        smallest counts encountered; combine with -o# to output graphs with a
        specific number of HISTs
  -o# : use only with -f; outputs the graphs with # HISTs
  -c  : outputs the complement of graphs which would normally be output
  -h  : print this helptext
  -v  : output extra information to stderr
```

### Examples
We assume file.g6 is a file containing graphs in graph6 format in the current working directory. All this commands should be executed from the directory which contains the histChecker binary.

`./histChecker < file.g6`
Outputs all graphs from file.g6 which are HIST-free to stdout. 
(Append `1> newfile.g6` to the command to write the contents of stdout to a file called newfile.g6.)

`./histChecker -c < file.g6`
Outputs all graphs from file.g6 which contain a HIST to stdout.

`./histChecker -1 < file.g6`
Outputs all graphs from file.g6 which are HIST-critical to stdout.

`./histChecker -1 -c < file.g6`
Outputs all graphs from file.g6 which are not HIST-critical to stdout.

`./histChecker -f < file.g6`
Counts the number of HISTs in each of the graphs of file.g6 and keeps track of the three fewest numbers encountered and how many times these numbers were attained as well as giving an example graph for each number of HISTs. Nothing is sent to stdout.

`./histChecker -f -o10 < file.g6`
Does the same as above, but sends all graphs from file.g6 with exactly 10 HISTs to stdout.
