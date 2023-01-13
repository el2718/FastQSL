# FastQSL

Since the server of http://staff.ustc.edu.cn/~rliu/qfactor.html will stop, those codes will be updated on this site from July 2022.

A GPU version is also provided https://github.com/peijin94/FastQSL.

-----------------------------
## Dependencies
### For *.pro
* IDL https://www.l3harrisgeospatial.com/Software-Technology/IDL or
* GDL https://gnudatalanguage.github.io
  * note: setting an environmental variable of GDL_PATH is necessary for write_png. 
If you used one of the binary packages available for Linux then it depends on the distribution:  
-Ubuntu & Fedora:  /usr/share/gnudatalanguage/lib   
-ArchLinux: /usr/lib/gdl  
-Gentoo: /usr/local/share/gdl   
**Please append such line to '~/.bashrc'** (for Ubuntu)    
export GDL_PATH=/usr/share/gnudatalanguage/lib 

### For *.f90
* ifort https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran or
* gfortran https://gcc.gnu.org/wiki/GFortran

## Usage
Please see the beginning of qfactor.pro  
* Demos:  
  IDL> .run demo_charge4.pro

## Cite as

* [Zhang, P., Chen, J.*, Liu, R. and Wang, C., 2022, FastQSL: A Fast Computation Method for Quasi-separatrix Layers. The Astrophysical Journal, 937, 26](https://iopscience.iop.org/article/10.3847/1538-4357/ac8d61)

```bibtex
@ARTICLE{2022ApJ...937...26Z,
       author = {{Zhang}, PeiJin and {Chen}, Jun and {Liu}, Rui and {Wang}, ChuanBing},
        title = "{FastQSL: A Fast Computation Method for Quasi-separatrix Layers}",
      journal = {\apj},
     keywords = {Solar magnetic fields, GPU computing, 1503, 1969},
         year = 2022,
        month = sep,
       volume = {937},
       number = {1},
          eid = {26},
        pages = {26},
          doi = {10.3847/1538-4357/ac8d61},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022ApJ...937...26Z},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
