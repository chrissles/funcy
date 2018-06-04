## Test environments

### rhub: 
* Ubuntu Linux 16.04 LTS, R-release, GCC
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit

### local: 
* R version 3.5.0 (2018-04-23), platform: x86_64-w64-mingw32 (64-bit), Win 7
* R version 3.4.4 (2018-03-15), platform: x86_64-pc-linux-gnu (64-bit), Debian

## R CMD check results

0 errors | 0 warnings | 3 notes


    N  checking CRAN incoming feasibility
       Maintainer: 'Christina Yassouridis <chris.yassou@gmail.com>'
   
       New maintainer:
         Christina Yassouridis <chris.yassou@gmail.com>
       Old maintainer(s):
         Christina Yassouridis <c_ya02@yahoo.de>
     
This change of email address was announced by c_ya02@yahoo.de to CRAN@R-project.org on 2018-06-01.      

    * checking compiled code ... NOTE
    File ‘funcy/libs/funcy.so’:
      Found no calls to: ‘R_registerRoutines’, ‘R_useDynamicSymbols’

    It is good practice to register native routines and to disable symbol
    search.

    See ‘Writing portable packages’ in the ‘Writing R Extensions’ manual.

I will take care of this note in the next submission.


    * checking examples ...
    ** running examples for arch 'i386' ... NOTE
    Examples with CPU or elapsed time > 5s
           user system elapsed
    funcit 5.15   0.08    5.24
    plot   5.10   0.04    5.14
    ** running examples for arch 'x64' ... NOTE
    Examples with CPU or elapsed time > 5s
           user system elapsed
    funcit 5.60   0.08    5.69
    plot   5.02   0.06    5.08
     
    * checking examples ... NOTE
    Examples with CPU or elapsed time > 5s
                user system elapsed
    rIMethods 18.252  0.128  18.352
    funcit     9.036  0.344   9.921

Execution times for both platforms are only slightly above the threshold of 5s, the linux machine is just slow.


## Reverse dependencies

This package has no reverse dependencies. 
