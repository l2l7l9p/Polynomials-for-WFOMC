# Weak Connectedness Polynomial and Strong Connectedness Polynomials

This tool is for computing Weak Connectedness Polynomial (WCP), Strict Strong Connectedness Polynomial (SSCP), Non-strict Strong Connectedness Polynomial (NSCP), Extended WCP and Tutte Polynomial from a given sentence.

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
  --func {wcp,nscp,sscp,ewcp,tutte}, -f {wcp,nscp,sscp,ewcp,tutte}
                        the function wanted
  --pred PRED, -p PRED  the special binary predicate
  --debug
```

E.g., to get the WCP for the predicate `E` in the sentence encoding an arbitrary undirected graph, the command is:

```
$ python sampling_fo2/main.py -i models/arbitrary-undirected-graph.wfomcs -f wcp -p E
```

## References

```
@article{DBLP:journals/ai/WangPWK24,
  author       = {Yuanhong Wang and
                  Juhua Pu and
                  Yuyi Wang and
                  Ondrej Kuzelka},
  title        = {Lifted algorithms for symmetric weighted first-order model sampling},
  journal      = {Artif. Intell.},
  volume       = {331},
  pages        = {104114},
  year         = {2024},
  url          = {https://doi.org/10.1016/j.artint.2024.104114},
  doi          = {10.1016/J.ARTINT.2024.104114},
  timestamp    = {Fri, 31 May 2024 21:06:28 +0200},
  biburl       = {https://dblp.org/rec/journals/ai/WangPWK24.bib},
  bibsource    = {dblp computer science bibliography, https://dblp.org}
}
```

```
@inproceedings{DBLP:conf/lics/WangP0K23,
  author       = {Yuanhong Wang and
                  Juhua Pu and
                  Yuyi Wang and
                  Ondrej Kuzelka},
  title        = {On Exact Sampling in the Two-Variable Fragment of First-Order Logic},
  booktitle    = {{LICS}},
  pages        = {1--13},
  year         = {2023},
  url          = {https://doi.org/10.1109/LICS56636.2023.10175742},
  doi          = {10.1109/LICS56636.2023.10175742},
  timestamp    = {Thu, 20 Jul 2023 11:32:59 +0200},
  biburl       = {https://dblp.org/rec/conf/lics/WangP0K23.bib},
  bibsource    = {dblp computer science bibliography, https://dblp.org}
}
```
If you use the WFOMC code, please cite
```
@inproceedings{DBLP:conf/uai/BremenK21,
  author       = {Timothy van Bremen and
                  Ondrej Kuzelka},
  editor       = {Cassio P. de Campos and
                  Marloes H. Maathuis and
                  Erik Quaeghebeur},
  title        = {Faster lifting for two-variable logic using cell graphs},
  booktitle    = {Proceedings of the Thirty-Seventh Conference on Uncertainty in Artificial
                  Intelligence, {UAI} 2021, Virtual Event, 27-30 July 2021},
  series       = {Proceedings of Machine Learning Research},
  volume       = {161},
  pages        = {1393--1402},
  publisher    = {{AUAI} Press},
  year         = {2021},
  url          = {https://proceedings.mlr.press/v161/bremen21a.html},
  timestamp    = {Fri, 17 Dec 2021 17:06:27 +0100},
  biburl       = {https://dblp.org/rec/conf/uai/BremenK21.bib},
  bibsource    = {dblp computer science bibliography, https://dblp.org}
}
```