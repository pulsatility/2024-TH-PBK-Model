# Spatially Dependent Tissue Distribution of Thyroid Hormones by Plasma Thyroid Hormone Binding Proteins
Anish D. Bagga<sup>1</sup>, Brian P. Johnson<sup>2</sup>, and Qiang Zhang<sup>3</sup>

Published in Pflügers Archiv - European Journal of Physiology (2025) https://doi.org/10.1007/s00424-024-03060-6

1. Emory College of Arts and Sciences, Emory University, Atlanta, GA 30322, USA

2. Department of Pharmacology and Toxicology, Michigan State University, East Lansing, MI 48824, USA

3. Gangarosa Department of Environmental Health, Rollins School of Public Health, Emory University, GA 30322, USA

  
**Abstract:**
Plasma thyroid hormone (TH) binding proteins (THBPs), including thyroxine-binding globulin (TBG), transthyretin (TTR), and albumin (ALB), carry THs to extrathyroidal sites, where THs are unloaded locally and then taken up via membrane transporters into the tissue proper. The respective roles of THBPs in supplying THs for tissue uptake are not completely understood. To investigate this, we developed a spatial human physiologically based kinetic (PBK) model of THs, which produces several novel findings. (1) Contrary to postulations that TTR and/or ALB are the major local T4 contributors, the three THBPs may unload comparable amounts of T4 in Liver, a rapidly perfused organ; however, their contributions in slowly perfused tissues follow the order of abundances of T4TBG, T4TTR, and T4ALB. The T3 amounts unloaded from or loaded onto THBPs in a tissue acting as a T3 sink or source respectively follow the order of abundance of T3TBG, T3ALB, and T3TTR regardless of perfusion rate. (2) Any THBP alone is sufficient to maintain spatially uniform TH tissue distributions. (3) The TH amounts unloaded by each THBP species are spatially dependent and nonlinear in a tissue, with ALB being the dominant contributor near the arterial end but conceding to TBG near the venous end. (4) Spatial gradients of TH transporters and metabolic enzymes may modulate these contributions, producing spatially invariant or heterogeneous TH tissue concentrations depending on whether the blood-tissue TH exchange operates in near-equilibrium mode. In summary, our modeling provides novel insights into the differential roles of THBPs in local TH tissue distribution.

**Keywords:** Thyroid hormone, Thyroxine-binding globulin, Transthyretin, Albumin, Tissue distribution, Gradient 



#  MATLAB Code
- TH_PBK_Spatial_CMD.m: Main MATLAB code to run the models and generate figures and tables.
- TH_PBK_Spatial_ODE.m: MATLAB ODE code to be called by TH_PBK_Spatial_CMD.m.
- TH_PBK_Spatial_ODE_Fig8.m: MATLAB ODE code to be called by TH_PBK_Spatial_CMD.m to generate Fig. 8.
- TH_PBK_Spatial_ODE_gradient.m: MATLAB ODE code to be called by TH_PBK_Spatial_CMD.m to generate Fig. 10-12, S8-S11.
- get_gradient_value.m: MATLAB code to be called by TH_PBK_Spatial_CMD.m to generate Fig. 10-12, S8-S11.
- Matlab_code_FigS1.m: MATLAB code to generate Fig. S1.
