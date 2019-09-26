```
                                                 .                 
                            ____     _          ":"                
  _ __   ___ _ __ ___   ___|___ \ __| |       ___:____     |"\/"|  
 | '_ \ / _ \ '_ ` _ \ / _ \ __) / _` |     ,'        `.    \  /   
 | | | |  __/ | | | | | (_) / __/ (_| |     |  O        \___/  |   
 |_| |_|\___|_| |_| |_|\___/_____\__,_|   ~^~^~^~^~^~^~^~^~^~^~^~^~
                                            ~^~^ ~^~ ^~^~^ ~^~^    
```

# nemo2d (beta v0.9)#

It is a lightweight, easy to understand 2D hydro code which illustrates the usage
of the 'Discontinuous Galerkin Spectral Element Method' (DGSEM) scheme.

For reference and quick insights into the datastructures a simple 'Finite-Volume'
scheme (FV) is provided as well.

The midterm goal to grow this code into an 'easy-to-comprehend',
'copy-and-paste-ready' snippet library for integration into bigger projects or as
starting point for small scale projects like teaching/workshops/thesis-writing, etc.

The manageable code size invites for easy fiddling and experimentation.

Happy coding!

## Installation ##

```bash
SHELL-PROMPT> git clone https://github.com/jmark/nemo2d.git
SHELL-PROMPT> cd nemo2d # ... and modify Makefile to your needs
SHELL-PROMPT> make
```
## Usage ##

```bash
SHELL-PROMPT> ./build/nemo2d
```

## Parallization ##

Parallization is provided via OpenMP.

```bash
SHELL-PROMPT> export OMP_NUM_THREADS=${NUMBER_OF_THREADS}
SHELL-PROMPT> ./build/nemo2d
```

## Lizence ##
```
Copyright (c) 2019 Johannes Markert <johannes.markert@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
