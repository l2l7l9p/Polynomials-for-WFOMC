# Weak Connectedness Polynomial and Strong Connectedness Polynomial

This tool is for computing Weak Connectedness Polynomial (WCP), Strong Connectedness Polynomial (SCP), Extended WCP and Tutte Polynomial from a given formula.

## Input format

1. First-order sentence with at most two logic variables, see [fol_grammar.py](sampling_fo2/parser/fol_grammar.py) for details, e.g.,
   * `\forall X: (\forall Y: (R(X, Y) <-> Z(X, Y)))`
   * `\forall X: (\exists Y: (R(X, Y)))`
   * `\exists X: (F(X) -> \forall Y: (R(X, Y)))`
   * ..., even more complex sentence...
2. Domain: 
   * `domain=3` or
   * `domain={p1, p2, p3}`
3. Weighting (optional): `positve_weight negative_weight predicate`
4. Cardinality constraint (optional): 
   * `|P| = k`
   * `|P| > k`
   * `|P| >= k`
   * `|P| < k`
   * `|P| <= k`
   * ...

### Example input file

2 colored graphs:

```
\forall X: (\forall Y: ((E(X,Y) -> E(Y,X)) &
                        (R(X) | B(X)) &
                        (~R(X) | ~B(X)) &
                        (E(X,Y) -> ~(R(X) & R(Y)) & ~(B(X) & B(Y)))))

V = 10
```

2 regular graphs:

```
\forall X: (~E(X,X)) &
\forall X: (\forall Y: ((E(X,Y) -> E(Y,X)) &
                        (E(X,Y) <-> (F1(X,Y) | F2(X,Y))) &
                        (~F1(X, Y) | ~F2(X,Y)))) &
\forall X: (\exists Y: (F1(X,Y))) & 
\forall X: (\exists Y: (F2(X,Y)))

V = 6
|E| = 12
```

> **Note: You can also directly input the SC2 sentence**

2 regular graphs (sc2):

```
\forall X: (~E(X,X)) &
\forall X: (\forall Y: (E(X,Y) -> E(Y,X))) &
\forall X: (\exists_{=2} Y: (E(X,Y)))

V = 6
```

Sampling possible worlds from `friends-smokes` MLN:

```
\forall X: (~fr(X,X)) &
\forall X: (\forall Y: (fr(X,Y) -> fr(Y,X))) &
\forall X: (\forall Y: (aux(X,Y) <-> (fr(X,Y) & sm(X) -> sm(Y)))) &
\forall X: (\exists Y: (fr(X,Y)))

person = 10
2.7 1 aux
```

> **Note: You can also directly input the MLN in the form defined in [mln_grammar.py](sampling_fo2/parser/mln_grammar.py)**

```
~friends(X,X).
friends(X,Y) -> friends(Y,X).
2.7 friends(X,Y) & smokes(X) -> smokes(Y)
\forall X: (\existes Y: (fr(X,Y))).
# or 
\exists Y: (fr(X,Y)).

person = 10
```

More examples are in [models](models/)

### Installation

Install the package:

```
$ pip install -e .
```

### How to use

see `$ python sampling_fo2/main.py -h`. The usage is as follows.

```
usage: main.py [-h] --input INPUT [--output_dir OUTPUT_DIR] [--func {wcp,scp,ewcp,tutte}] --pred PRED [--debug]

Computing Polynomials from a C2 sentence

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        mln file
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
  --func {wcp,scp,ewcp,tutte}, -f {wcp,scp,ewcp,tutte}
                        the function wanted
  --pred PRED, -p PRED  the special binary predicate
  --debug
```

E.g., to get the WCP for the predicate `E` in the sentence encoding an arbitrary undirected graph, the command is:

```
$ python sampling_fo2/main.py -i models/arbitrary-undirected-graph.wfomcs -f wcp -p E
```

## References