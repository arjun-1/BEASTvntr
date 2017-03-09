[![Build Status](https://travis-ci.org/arjun-1/BEASTvntr.svg?branch=master)](https://travis-ci.org/arjun-1/BEASTvntr) [![GitHub version](https://badge.fury.io/gh/arjun-1%2FBEASTvntr.svg)](https://badge.fury.io/gh/arjun-1%2FBEASTvntr) [![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://raw.githubusercontent.com/arjun-1/BEASTvntr/master/COPYING.LESSER)
# BEASTvntr

A package for [BEAST2](https://github.com/CompEvol/beast2) 2.4.3 or higher, which can infer phylogeny for VNTR (=Variable Number of Tandem Repeat) data.

## Installation 
To install BEASTvntr from the BEAST2 Package Manager in Beauti, go to **File > Manage Packages**, select *BEASTvntr* and click **Install/Upgrade**.

To install BEASTvntr manually, download the [latest release](https://github.com/arjun-1/BEASTvntr/releases/download/v0.1.1/BEASTvntr.addon.v0.1.1.zip) and extract the contents into its own folder in `~/.beast/2.4/`.
## Building from Source

These instructions will get you a copy of the BEASTvntr project up and running for development and testing purposes.

### Prerequisities

Make sure Apache Ant is installed.

### Installing

Download and build the BEAST2 project

```bash
git clone https://github.com/CompEvol/beast2.git
cd beast2
ant
```
Download and build the package
```bash
cd ../
git clone https://github.com/arjun-1/BEASTvntr.git
cd BEASTvntr
ant addon
```
Install the package
```bash
cp -r release/add-on ~/.beast/2.4/BEASTvntr
```
## Example
These instructions will show how to infer phylogeny for VNTR data of a set of taxa provided in a [paper](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0007815) by Comas.
### Setting up the XML file
First download [comas2009_VNTR.csv](https://raw.githubusercontent.com/arjun-1/BEASTvntr/master/examples/csv/comas2009_VNTR.csv) which contains the repeats in CSV format. Start Beauti, either via its shortcut or by running `java -cp build/dist/launcher.jar beast.app.beauti.BeautiLauncher` in `beast2/`. In the Beauti window, click **File > Import Alignment** and select *comas2009_VNTR.csv*. In the window that appears, we can either select repeats (single partition), repeats (multiple partitions) or nucleotides to import. Select *Repeats (single partition)* and click **OK**.

After selecting repeats, we must specify the minimum and maximum repeat which will bound our state space. For *Minimum repeat* specify **1** and for *Maximum repeat* specify **15**, and click **OK**.

We can now specify the parameters of the substitution model. Click on the **Site Model** tab and set *Gamma Category Count* to **6** to allow rate variations among the different loci. Then, click on **estimate** for the *Shape* parameter to estimate the gamma shape parameter of the distribution of rates among sites.

Now, click on the **Priors** tab. For the coalescent model, select **Coalescent Constant Population**

To save the configuration in an XML file, click **File > save** and save as **comas2009_VNTR.xml**.
### Running BEAST2
Start BEAST2, either via its shortcut or by running `java -cp build/dist/launcher.jar beast.app.beastapp.BeastLauncher` in `beast2/`, and select comas2009_VNTR.xml. If everything went right, you should see the MCMC run starting:
```text
Start likelihood: -6033.8021475437245 
Writing file comas2009_VNTR.log
Writing file comas2009_VNTR.trees
         Sample      posterior ESS(posterior)     likelihood          prior
              0     -6024.6925              N     -5830.8518      -193.8406 --
           1000     -3573.9178         2.0        -3260.2564      -313.6613 --
           2000     -3320.7594         3.0        -3007.2831      -313.4763 --
           3000     -3155.1852         4.0        -2849.9808      -305.2044 --
```

## Background
To infer phylogeny, BEASTvntr uses a model described in a [paper](http://www.genetics.org/content/168/1/383.long) by Sainudiin. Of this model, several implementations are available.

The standard implementation is called *Sainudiin Vanilla* and can model mutational bias, mutation rate proportionality, and any multi-step mutations. *Sainudiin Computed Frequencies Vanilla* is a variant on the above model, where the frequencies of the states in the root node are given by the stationary distribution, which is calculated from the other model parameters. 

All these implementations use a modified for expression for the mutational bias `beta`, which is described in a [paper](http://www.genetics.org/content/188/1/151.long) by Wu. In this expression, the bias `beta` for expansion given a mutation event, depends on the parameters `b0, b1`. The implementations *Sainudiin* and *Sainudiin Computed Frequencies* are an adaptation of the standard *Vanilla* models, which use a transformation of these parameters, given by:
```
b0 =  biasMagnitude / sqrt(1 + 1 / (focalPoint - minimum repeat)^2)
b1 = -biasMagnitude / sqrt(1 + (focalPoint - minimum repeat)^2)
```
This was done so that the equation
```
beta(b0, b1, focalPoint) = 1 - beta(b0, b1, focalPoint)
```
is always satisfied, i.e. for that `focalPoint` the bias for expansion is equal to that of contraction. This means that `focalPoint` can intuitively be interpreted as the focal point of the mutational bias.

In addition, the proportionality of the mutation rate to the number of repeats, `a1`, has been transformed into:
```
oneOnA1 = 1 / a1
```
For the cases where the mutation rate of the minimum repeat is significantly lower than that of the other repeats, `a1` blows up, whilst `oneOnA1` remains constrained.

##Known Issues
During a MCMC run, it is possible that the likelihood makes a sudden unrealistic increase, and that the trace of estimated parameters becomes a flat line:

<img src="https://cloud.githubusercontent.com/assets/8102654/16612531/bd0c3032-4367-11e6-8b60-1873ff80aef8.png" alt="alt text" width="680" height="472">

The cause of this issue might be that too many parameters are being estimated in the model. If you encounter such an issue, doing any of the following might resolve it:  
* When using the beagle library, pass `-beagle_single` as an option to beast, or try other scaling options.  
* Use *Sainudiin Frequencies Computed* instead of *Sainudiin* as substitution model. The *Sainudiin Frequencies Computed*   does not estimate the frequencies of the repeats, thus this greatly reduces any over-parametrization.  
* Use more restrictive priors on any of the parameters `biasMagnitude, focalPoint, g, oneOnA1` of the model.  Bounding `oneOnA1` from above seems to be most helpful.  
* Remove any duplicate VNTR sequence in the imported alignment.

See [this](https://groups.google.com/forum/#!topic/beast-users/ScG6PEZTADE) forum post for more information on a similar problem.

## Author

* **Arjun Dhawan** - *Initial work* - [arjun-1](https://github.com/arjun-1)

## License

This project is licensed under version 3 of the GNU Lesser General Public License - see the [COPYING.LESSER](COPYING.LESSER) file for details.
