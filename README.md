# LARCC 

Linear Algebraic Representation to Compute with (Co)Chains 

"*With increased complexity of geometric data, topological models play an increasingly important role beyond boundary representations, assemblies, finite elements, image processing, and other traditional modeling applications. While many graph- and index- based data structures have been proposed, no standard representation has emerged as of now. Furthermore, such representations typically do not deal with representations of mappings and functions and do not scale to support parallel processing, open source, and client-based architectures. We advocate that a proper mathematical model for all topological structures is a (co)chain complex: a sequence of (co)chain spaces and (co)boundary mappings. This in turn implies all topological structures may be represented by a collection of sparse matrices. We propose a Linear Algebraic Representation (LAR) scheme for mod 2 (co)chain complexes using compressed sparse matrices and show that it supports variety of topological computations using standard matrix algebra, without any overhead in space or running time.*" (from [Linear algebraic representation for topological structures](http://www.sciencedirect.com/science/article/pii/S001044851300184X), CAD 2014)

<!-- 
## API

see [docs](https://github.com/cvdlab/plasm.js/blob/master/docs/Readme.md)
 -->

## Installation

### pre-requisite

install [pyplasm](https://github.com/plasm-language/pyplasm)

### install

type in a terminal, from your workspace directory

```
git clone git@github.com:cvdlab/lar-cc
cd larlib
python setup.py install 
```

### update

type in a terminal, from your `lar-cc` directory

```
git pull origin master
cd larlib
python setup.py install 
```

## Authors

- [Alberto Paoluzzi](http://paoluzzi.dia.uniroma3.it)
- [Giorgio Scorzelli](http://www.dia.uniroma3.it/~scorzell)

## Contributors

- Alberto Paoluzzi
- Giorgio Scorzelli

## License

(The MIT License)

Copyright (c) 2013 Alberto Paoluzzi for CVD Lab ([http://dia.uniroma3.it/~cvdlab](http://dia.uniroma3.it/~cvdlab))

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
'Software'), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

