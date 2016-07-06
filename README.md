# BEASTvntr

A package for BEAST2 which can infer phylogeny for VNTR (=Variable Number of Tandem Repeat) data.

## Pre-built Binaries
...
## Building from Source

These instructions will get you a copy of the BEASTvntr project up and running on for development and testing purposes.

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
cp release/add-on ~/.beast/2.4/BEASTvntr
```
## Example
These instructions will show how to infer phylogeny for VNTR data of a set of taxa provided in a paper of [Comas](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0007815)
### Setting up the XML file
First download [comas2009_VNTR.csv](examples/comas2009_VNTR.csv) which contains the repeats in CSV format. Start Beauti. In the Beauti window, click **File > Import Alignment** and select *comas_VNTR.csv*. In the window that appears, we can either select repeats or nucleotides to import. Select *Repeats* and click **OK**.

After selecting repeats, we must specify the minimum and maximum repeat which will bound our state space. For *Minimum repeat* specify **1** and for *Maximum repeat* specify **15**, and click **OK**.

We can now specify the parameters of the substitution model. Click on the **Site Model** tab and set *Gamma Category Count* to **6** to allow rate variations among the different loci. Then, click on **estimate** for the *Shape* parameter to estimate the gamma shape parameter of the distribution of rates among sites.

Now, click on the **Priors** tab. For the coalescent model, select **Coalescent Constant Population**

To save the configuration in an XML file, click **File > save** and save as **comas2009_VNTR.xml**.
### Running BEAST2
Start BEAST2 and select comas2009_VNTR.xml. If everything went right, you should see the MCMC run starting:
```text
Start likelihood: -6033.8021475437245 
Writing file comas_VNTR.log
Writing file comas_VNTR.trees
         Sample      posterior ESS(posterior)     likelihood          prior
              0     -6024.6925              N     -5830.8518      -193.8406 --
           1000     -3573.9178         2.0        -3260.2564      -313.6613 --
           2000     -3320.7594         3.0        -3007.2831      -313.4763 --
           3000     -3155.1852         4.0        -2849.9808      -305.2044 --
```

## Background
To infer phylogeny, BEASTvntr uses two implementations of a model explained in a paper of [Sainudiin](http://www.genetics.org/content/168/1/383.long). The first implementation is called *Sainudiin* and can model a mutational bias, mutation rate proportionality, and any multi-step mutations. The second implementation called *SainudiinStepWise*, is the same as *Sainudiin*, except that it only allows single step mutations, and that the frequencies of the repeats of the root node are calculated from other model parameters. 

These implementations use a modified for expression for the mutational bias beta however, which is described in a paper of [Wu](http://www.genetics.org/content/188/1/151.long). In this expression, the bias `beta` for expansion given a mutation event, depends on the parameters `b_0, b_1`. However, BEASTvntr uses a reparametrization of these parameters:
```
b_0 =  r_b / sqrt( 1 + 1 / i_eq^2 )
b_1 = -r_b / sqrt( 1 + i_eq^2 )
```
This was done so that the equation
```
beta(b_0, b_1, i_eq) = 1 - beta(b_0, b_1, i_eq)
```
is always satisfied, i.e. for a particular `i_eq` the bias for expansion is equal to that of contraction. This means that `i_eq` can intuitively be interpreted as the focal point of the mutational bias.
##Known Issues
During a MCMC run, it is possible that the likelihood makes a sudden unrealistic increase, and that the trace of estimated parameters becomes a flat line:

<img src="https://cloud.githubusercontent.com/assets/8102654/16612531/bd0c3032-4367-11e6-8b60-1873ff80aef8.png" alt="alt text" width="680" height="472">

The cause of this issue might be that too many parameters are being estimated in the model. If you encounter such an issue, doing any of the following might resolve it:  
1. Pass `-beagle_scaling none` as an option to beast.  
2. Use *SainudiinStepWise* instead of *Sainudiin* as substitution model. The *SainudiinStepWise*   does not estimate the frequencies of the repeats, thus this greatly reduces any over-parametrization.  
3. Use more restrictive priors on any of the parameters `r_b, i_eq, g, a_1` of the model.

See [this](https://groups.google.com/forum/#!topic/beast-users/ScG6PEZTADE) forum post for more information on a similar problem.

## Author

* **Arjun Dhawan** - *Initial work* - [arjun-1](https://github.com/arjun-1)

## License

This project is licensed under version 3 of the GNU Lesser General Public License - see the [COPYING.LESSER](COPYING.LESSER) file for details