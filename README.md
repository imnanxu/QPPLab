# A generalized toolbox for QPP detection, analyzsis, and visulization for rodent or human brains
This is a generally applicable MATLAB toolbox, which detects, analyzes, and visualizes up to 5 QPPs and QPP regressed functional connectivity maps from fMRI timeseries of the brain across rodents and humans. This toolbox is a significant modification and extension of [QPP_Scripts_v0620](https://github.com/GT-EmoryMINDlab/QPP_Scripts_v0620), and can analyze QPPs for resting brain and for task-evoked/stimulated brains across species.

# Table of Contents
* 1 - [Prerequisite & Resources](#section-1)
    * 1.1 [Prerequisite](#section-1-1)    
    * 1.2 [Input file (./Inputs/)](#section-1-2)
        * 1.2.1 [Parameters to include (*.mat)](#section-1-2-1)
        * 1.2.2 [(optional) st0_ROIreOrg.m](#section-1-2-2)
    * 1.3 [Resources (./resources/)](#section-1-2)	
* 2 - [Main Scripts](#section-2)
    * 2.1 [(Step 1) Run 'st1_ParamsSet.m'](#section-2-1)
    * 2.2 [(Step 2) Run 'st2_QPPanalysis.m'](#section-2-2)
        * 2.2.1 [Parameters to be prespecified](#section-2-3-1)
        * 2.2.2 [Automated QPP analysis](#section-2-3-2)
    * 2.3 [(Step 3) Run 'st3_QPPFCvisual.m'](#section-2-3)
        * 2.3.1 [EPI registration estimation & wm/csf mask generation](#section-2-3-1)
        * 2.3.2 [Non-brain tissue noise estimation by PCA](#section-2-3-2)
* 3 - [Input File](#section-3)
* 4 - [Output Files](#section-4)
* 4 - [References](#section-4)


<a name="section-1"></a>
## 1. Dependencies
1. [FSL5.0](https://web.mit.edu/fsl_v5.0.10/fsl/doc/wiki/FslInstallation(2f)Linux.html) (Jenkinson et al., 2012), [AFNI](https://afni.nimh.nih.gov/download) (Cox, 1996; Cox & Hyde, 1997), and [ANTs](https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS) (Avants et al., 2022): A Linux (i.e., Ubuntu, Fedora, and Red Hat) or MacOS system with the above 3 software packages installed is required. For Windows systems, it's possible to install the three packages through a Windows Subsystem for Linux (WSL, see "SoftwareInstallation_fsl_afni_ants.txt" for more details).
2. [Matlab](https://www.mathworks.com/) (The Mathworks Inc., Natick, MA, USA, R2018a or a later version) and [NIfTI and ANALYZE toolbox](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image) (Chen, 2022) are required for calling PCNN3D (Chou et al., 2011), which is superior for mouse brain mask creation (see Section 4.1.4 for more details). 
The toolbox has been cloned to this repository in *NIfTI_toolbox* for convenience.

    *Supported Matlab Operating Systems:* Matlab software is supported in Windows (10, 11, and Server 2019) as well as MacOS and Linux (i.e., Ubuntu, Debian, RedHat, SUSE). For the full Linux system requirements, please refer to the [official documentation](https://www.mathworks.com/support/requirements/matlab-linux.html). If installing WSL using a Linux distribution other than Ubuntu or Debian as described in *SoftwareInstallation_fsl_afni_ants.txt*, replace all `apt` and `apt-get` commands with the equivalent command for your OS package manager (e.g., [zypper](https://en.opensuse.org/SDB:Zypper_usage) for SUSE).

    *Running Without Matlab Support:* By default, in *preproc_script_1.sh*, if WSL isn't detected, the default Matlab directory is set to `matlab`. Override this by passing a `--matlab_dir` argument in 
the CLI. To run the first script without Matlab or PCNN3D, set the `--matlab_dir` argument to `NA`.

p2data=['./Input/' data '.mat']; % path of data file which must includes following parameters:
% D0:   a nsbj X nscn cell matrix. Each cell has a nroi X ntimepoints 
%       matrix of EPI timeseries
% MotionInf: a nsbj X nscn cell matrix. Each cell >=1 segments of
%       timepoints without significant motions.
% nY:   # of total networks
% G2Y:  a (nroi X 1) vector including the network# of each ROI
% ibY:  a (nY+1 X 1) vector including the last ROI label of each 
%       functional network;  always set ibY(1)=0.
% iG2Y: a (nY X 1) cell vector including the list of ROI labels for each network.
% YLB: a (nY X 1) cell vector including the shorthand label for each network
% %%%The last 5 variables can be generated from the original data and the
% atlas-network file (see ./resources)




<a name="section-2"></a>
## 2. Main Pipeline
<a name="section-2-1"></a>
### 2.1 (Step 1) Run 'st1_ParamsSet.m'
This is for setting up the intial parameters for the QPP analysis in step 2. The following variables will be predefined, and a parameter file Params_`data`\_`ext`.mat will be generated after running this script.
|      Purpose     |  Variable name | Description | Note   | 
|------------------|-----------------|--------|-------------|
|  Filepath  		|`data`   | the input filename |The input should has the filename `data`.mat |
|                  	|`ext`    | filename extension for the parameter file| The parameter filename will be  Params_`data`\_`ext`.mat |
|  QPP global parameters|`nP`     | total # of QPPs to detect (nP<=5)| If nP=1, only detect the primary QPP (QPP1); if nP=2, detect both QPP1 & QPP2; etc.|
|		   	|`PL`     | a (nP X 1) vector of QPP window length | ~20s for humans (e.g., PL(ip)=20/TR), |
|  QPP detection	|`cth13` & `cth34`     | a 2D vector of correlation threshold for QPP1-QPP3 (`cth13`) & for QPP4-QPP5 (`cth45`)| If you do not need to detect QPP4-QPP5, please assign `cth34` a random number (e.g., `cth34`=[0, 0]).|
|  QPP phase adjustment	|`cthph` | similarity threshold when phase-adjusting (phadj) a QPP |Default value: cthph=0.88|
|		   	|`s`     | control for strict phase adjustment (`s`=1) or relaxed phase adjustment (`s`=0)||
|		   	|`sdph`     | a (nP X 1) cell array of reference parcels| Each cell may include >=1 parcel IDs. The phase adjusted QPP waveform will start from rising positive values for the selected parcels.|
|  Functional connectivity (FC) analysis|`fz` | control for the output matrix `FCr` to be the pearson correlation (`fz`=1) or to be the Fisher Z-Transformation of the pearson correlaion (`fz`=1).|
<a name="section-2-2"></a>
### 2.2 (Step 2) Run 'st2_QPPanalysis.m'
This is for detecting and analyzing QPPs based on the parameters setup in step 1. 
<a name="section-2-1-1"></a>
#### 2.2.1 Parameters to be prespecified        
The following three parameters needs to be prespecified at the beginning of this script.
|      Purpose     |  Variable name  | Description | Note   | 
|------------------|-----------------|-------------|--------|
|  Filepath  		|`dataext`   | parameter filename |The parameter .mat file generated from step 1, which has the filename Param_`dataext`.mat |
|  Data concatenation method |`runM`     | control the way to concatenate the data| If `runM`=1, concatenate all D{i,j} as a whole group and detect group QPP; if `runM`=2, concatenate all D{i,:} and detect QPP from all scans of each subject; if `runM`=3, concatenate all D{:,:} and detect QPP from all subjects of each scan.|
| QPP detection	method|`rbstScrn`     | control for fast QPP dectection (`rbstScrn`=0) or robust QPP detection (`rbstScrn`=1)|The fast QPP detection selectes a limited number of starting points which was used in XXXX, whereas the robust detection selects all possible starting points which was used in (XXX).|

<a name="section-2-1-1"></a>
#### 2.2.2 Automated QPP analysis
The following analytical procedures will be executed.


<a name="section-3"></a>
## 3. Input File(./Inputs/)
<a name="section-3-1"></a>	
### 3.1 Input variables (`data`.mat)
The input file should has the filename `data`.mat, which includes the following variables:


### 3.2 (optional step) run 'st0_ROIreOrg.m' for variable generations
<a name="section-4-3-3"></a>
#### 4.3.3 Nuisance regressions: 26 possible regressors (Chuang et al., 2018) & a user specified file of regressors
    a. 3 for detrends: constant, linear, and quadratic trends
    b. 10 PCs from non-brain tissues
    c. 6 motion regressors (based on motion correction results) 
    d. 6 motion derivative regressors: the temporal derivative of each motion regressor
    e. csf or/and wmcsf signal(s)
    f. one *.txt file containing user specified regressors (e.g., task patterns to be regressed)
The script also generates a default outputs (0EPI_\*) which only regresses out the 3 trends (a). One can specify any combinations of above regressors in the command line. If you are preprocessing a group dataset, the same combination of regressors will be applied to all data.

<a name="section-4-3-4"></a>
#### 4.3.4 Normalization & temporal filtering
    a. Normalize the regressed signals
    b. Bandpass filter: bandwidth depends on the use of anesthesia
    	e.g., 0.01–0.1Hz for iso and 0.01–0.25Hz for dmed (Pan et al., 2013)
Output: \_mc_c_norm_fil

<a name="section-4-3-5"></a>
#### 4.3.5 EPI template registration & spatial smoothing & seed extraction
    a. EPI template registration: transform cleaned-up data to template space by the transformation matrix estimated in (2.a)
    b. Use Gaussian kernel for spatial smoothing. Set the FWHM value in mm in command line options
    c. Extract the averaged timeseries based on atlas.
Output: \_mc_c_norm_fil_reg_sm, \_mc_c_norm_fil_reg_sm_seed.txt

In the data sample folder, the functional connectivity map (FC.tif) generated by Matlab in our post analysis using the preprocessed timeseries is also provided.

<a name="section-5"></a>
## 5. References
Avants, B., Tustison, N. J., & Song, G. (2022). Advanced Normalization Tools: V1.0. The Insight Journal. https://doi.org/10.54294/UVNHIN

Barrière, D. A., Magalhães, R., Novais, A., Marques, P., Selingue, E., Geffroy, F., Marques, F., Cerqueira, J., Sousa, J. C., Boumezbeur, F., Bottlaender, M., Jay, T. M., Cachia, A., Sousa, N., & Mériaux, S. (2019). The SIGMA rat brain templates and atlases for multimodal MRI data analysis and visualization. Nature Communications, 10(1), 1–13. https://doi.org/10.1038/s41467-019-13575-7

Chou, N., Wu, J., Bai Bingren, J., Qiu, A., & Chuang, K. H. (2011). Robust automatic rodent brain extraction using 3-D pulse-coupled neural networks (PCNN). IEEE Transactions on Image Processing : A Publication of the IEEE Signal Processing Society, 20(9), 2554–2564. https://doi.org/10.1109/TIP.2011.2126587

Chuang, K.-H., Lee, H.-L., Li, Z., Chang, W.-T., Nasrallah, F. A., Yeow, L. Y., & Singh, K. K. D. /O. R. (2018). Evaluation of nuisance removal for functional MRI of rodent brain. NeuroImage. https://doi.org/10.1016/J.NEUROIMAGE.2018.12.048

Cox, R. W. (1996). AFNI: Software for analysis and visualization of functional magnetic resonance neuroimages. Computers and Biomedical Research, 29(3), 162–173. https://doi.org/10.1006/cbmr.1996.0014

Cox, R. W., & Hyde, J. S. (1997). Software tools for analysis and visualization of fMRI data. NMR in Biomedicine, 10(4–5), 171–178. https://doi.org/10.1002/(SICI)1099-1492(199706/08)10:4/5<171::AID-NBM453>3.0.CO;2-L

Jenkinson, M., Beckmann, C. F., Behrens, T. E. J., Woolrich, M. W., & Smith, S. M. (2012). FSL. NeuroImage, 62(2), 782–790. https://doi.org/10.1016/j.neuroimage.2011.09.015

Lee, H. L., Li, Z., Coulson, E. J., & Chuang, K. H. (2019). Ultrafast fMRI of the rodent brain using simultaneous multi-slice EPI. NeuroImage, 195, 48–58. https://doi.org/10.1016/j.neuroimage.2019.03.045

Lein, E. S., Hawrylycz, M. J., Ao, N., Ayres, M., Bensinger, A., Bernard, A., Boe, A. F., Boguski, M. S., Brockway, K. S., Byrnes, E. J., Chen, L., Chen, L., Chen, T.-M., Chi Chin, M., Chong, J., Crook, B. E., Czaplinska, A., Dang, C. N., Datta, S., … Jones, A. R. (2006). Genome-wide atlas of gene expression in the adult mouse brain. Nature 2006 445:7124, 445(7124), 168–176. https://doi.org/10.1038/nature05453

Pan, W.-J., Thompson, G. J., Magnuson, M. E., Jaeger, D., & Keilholz, S. (2013). Infraslow LFP correlates to resting-state fMRI BOLD signals. Neuroimage, 74(0), 288–297. https://doi.org/10.1016/j.neuroimage.2013.02.035

Thompson, G. J., Pan, W. J., Magnuson, M. E., Jaeger, D., & Keilholz, S. D. (2014). Quasi-periodic patterns (QPP): Large-scale dynamics in resting state fMRI that correlate with local infraslow electrical activity. NeuroImage, 84, 1018–1031. https://doi.org/10.1016/j.neuroimage.2013.09.029

Chen, J. (2022). Tools for NIfTI and ANALYZE image. MATLAB Central File Exchange. https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

Xu, N., LaGrow, T. J., Anumba, N., Lee, A., Zhang, X., Yousefi, B., Bassil, Y., Clavijo, G. P., Khalilzad Sharghi, V., Maltbie, E., Meyer-Baese, L., Nezafati, M., Pan, W.-J., & Keilholz, S. (2022). Functional Connectivity of the Brain Across Rodents and Humans. Frontiers in Neuroscience, 0, 272. https://doi.org/10.3389/FNINS.2022.816331

