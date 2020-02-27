# Black Holes in Quadratic Gravity

Hyun Lim and Aaron Held

We study some possible dynamical black holes in 
Quadratic gravity. Out interests foucus on stability
problem and possible gravitational wave signals in this
theory. We consider 3+1 decomposition to enable dynamical
studies. Also, we analyze structure of the equations of
motion from the theory and possible effective field theory
apporoach to have linearized version of system. Numerical
simulations of this system would be performed by `Dendro`
toolkit.


## Get this project with the code

In this project, we are using `Dendro` framework. 
To get all the code with other materials, you need to
clone it via
```{engine=sh}
git clone --recursive https://github.com/hlim88/quadGrav.git
```

Note that the code exists in separate repo. I updated it
regularly to maintaint latest version but you also want to
check yourself to contain latest figure using below command
```{engine=sh}
git submodule foreach git pull origin master
```
This will call latest commit from the repo

### Build the code 

Detailed instruction to build our code is described in
[here](https://github.com/hlim88/Dendro-5.01/blob/master/INSTALL_QG.md)

### For developer

If you are interested in updating/modifying code, 
please contact Hyun Lim (hylim1988@gmail.com) to learn how
you can use submodule accordingly. Some part of development is 
not open-source yet.

## Directory description

`notes` directory contains theory/formulation for this project.
Type `make` to complie the note

`scripts` directory contains `Mathematica` file that peforms 3+1
decomposition of the equations. It is saved in `Wolfram Language`
format. Users need a Wolfram compatible program such as`Mathematica`
or `WolframOne` etc.

Several references explain about the project. You can find
references from [here](https://github.com/hlim88/quadGrav/wiki/References)

## Contacts

If you have any questions on this project, please contact Hyun Lim 
(hylim1988@gmail.com) or Aaron Held (held.aaron@gmail.com)
