## aimd

This is designed to help researchers to perform diffusion analysis from *ab inito* molecular dynamic 
(AIMD) simulations. It contains some useful python classes which can be integrated into your python
analysis code, and one easy-to-use executable python script which can perform diffusion analysis on
the terminal without further modification.


The executable python script (script/analyze_aimd.py) has two main functions:

1. Analyze diffusion properties from AIMD calculations at one single temperature. 

2. Analyze diffusion properties from AIMD calculations at multiple temperatures. 

    More information can be found by the behind code after you install the library

    ```bash
    analyze_aimd -h
    analyze_aimd diffusivity -h
    analyze_aimd arrhenius -h
    ```

## Citing

If you use this library, please cite the following papers:

    He, Xingfeng, Yizhou Zhu, Alexander Epstein, and Yifei Mo. "Statistical variances 
    of diffusional properties from ab initio molecular dynamics simulations." 
    npj Computational Materials 4, no. 1 (2018): 18.
    
http://dx.doi.org/10.1038/s41524-018-0074-y

    He, Xingfeng, Yizhou Zhu, and Yifei Mo. "Origin of fast ion diffusion in super-ionic 
    conductors." Nature communications 8 (2017): 15893.
    
http://dx.doi.org/10.1038/ncomms15893

    Mo, Yifei, Shyue Ping Ong, and Gerbrand Ceder. "First principles study of the Li10GeP2S12 
    lithium super ionic conductor material." Chemistry of Materials 24, no. 1 (2011): 15-17.
    
http://dx.doi.org/10.1021/cm203303y


    Ong, Shyue Ping, William Davidson Richards, Anubhav Jain, Geoffroy Hautier, Michael Kocher, 
    Shreyas Cholia, Dan Gunter, Vincent L. Chevrier, Kristin A. Persson, and Gerbrand Ceder. 
    "Python Materials Genomics (pymatgen): A robust, open-source python library for materials 
    analysis." Computational Materials Science 68 (2013): 314-319. 

https://doi.org/10.1016/j.commatsci.2012.10.028

## Install and testing steps

1. Install all dependency and this library
    
    ```bash
    python setup.py install
    ``` 
    
    if you have no root access, you may need to use
    
    ```bash
    python setup.py install --user
    ```
    
2. Try to import python classes in your python console

    ```python
    from aimd.diffusion import DiffusivityAnalyzer, ErrorAnalysisFromDiffusivityAnalyzer, \
    ArreheniusAnalyzer
    ```

3. The setup.py will automatically create an executable file analyze_aimd into your PATH. Try to
call it from terminal and read the documentations:

    ```bash
    analyze_aimd -h
    analyze_aimd diffusivity -h
    analyze_aimd arrhenius -h
    ```
    
4. Use the provided test files to perform diffusion analysis.
    
    a. go to the folder aimd/tests/tests_files/latp_md
    
    b. run in terminal
     ```bash
     analyze_aimd diffusivity Li+ RUN_ 10 29 3.2
     ``` 
    
    c. You can ignore the warning msgs. You will get the diffusion results. 
    The diffusivity is ~7.6e-5 cm^2/s, conductivity is ~583 mS/cm
    
    d. go to folder aimd/tests/tests_files/arrhenius 
    
    e. run in terminal
     ```bash
     analyze_aimd arrhenius D_T.csv -p POSCAR -T 300 -s Li+
     analyze_aimd arrhenius D_T.csv -p POSCAR -T 300 -s Li+ --plot
     ``` 
    
    f. You will get the arrhenius relationship of LATP. The conductivity at 300K is predicted to be 
    ~1.08 mS/cm, Ea is ~0.258 eV +- 0.017 eV. If your x11 windows settings are correct or you are at
    local computer, you will have a plot window pop up to show the arrhenius relationship for cmd with
    -plot option.

## License


Python library aimd is released under the MIT License. The terms of the license are as
follows:

    The MIT License (MIT) Copyright (c) 2018 UMD 
     
    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
    documentation files (the "Software"), to deal in the Software without restriction, including without limitation 
    the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and 
    to permit persons to whom the Software is furnished to do so, subject to the following conditions:
     
    The above copyright notice and this permission notice shall be included in all copies or substantial portions of 
    the Software.
     
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO 
    THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
    SOFTWARE.
