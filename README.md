# SEAICEtools.jl

## Julia Package to work with remote sensing data of Sea Ice.

<p align="center">
<img width="466" height="480" alt="mosaic0152_20191115T03271573_20191118T05301574_LF" src="https://github.com/user-attachments/assets/c91c4345-8907-4c2c-a2a0-d5cbe49cc6be" />

 *Figure of Sea ice lead fraction (LF) during the MOSAiC expedition. The red-star represent the RV Polarstern and the black line is the dift trajectory during the previos days. Adapted from [5].*
</p>

At this point only the following data providers are supported:
- Institute of Environmental Physics (IEP) at the University of Bremen, [1,2]:
  * AMSRE -AMSR2: SIC Arctic-3.125km @ 5 km
  * Merged AMSR2-MODIS: SIC @ 1 km
- SENTINEL-1A by ESA provided by [3]:
  * lead fraction (LF) and sea ice divergence (DIV) @ 0.7 km
- OSI SAF EUMETSAT
  * SSMIS, AMSR2 SIC @ 10 km

[To Uni Bremen data repository](https://seaice.uni-bremen.de/sea-ice-concentration/amsre-amsr2/)

### References to data sources where the package can be used:
[1] Spreen, G., Kaleschke, L., and Heygster, G.: Sea ice remote sensing using AMSR-E 89-GHz channels, Journal of Geophysical Research:
Oceans, 113, https://doi.org/https://doi.org/10.1029/2005JC003384, 2008.

[2] Ludwig, V., Spreen, G., and Pedersen, L. T.: Evaluation of a New Merged Sea-Ice Concentration Dataset at 1 km Resolution from Thermal Infrared and Passive Microwave Satellite Data in the Arctic, Remote Sens.-Basel, 12, 3183, https://doi.org/10.3390/rs12193183, 2020.

[3] von Albedyll, L., Hendricks, S., Hutter, N., Murashkin, D., Kaleschke, L., Willmes, S., Thielke, L., Tian-Kunze, X., Spreen, G., and Haas, C.: Lead fractions from SAR-derived sea ice divergence during MOSAiC, The Cryosphere Discuss. [preprint], https://doi.org/10.5194/tc-2023-123, in review, 2023.

### The package has been used or for the following research:
[4] Saavedra Garfias, P., and H. Kalesse-Los: Observed modulation of wintertime Western Arctic mixed-phase cloud properties by sea ice conditions, their long-term variabilities and trends, EGUsphere, DOI: 10.5194/egusphere-2025-2327, 2025.

[5] Saavedra Garfias, P., H. Kalesse-Los, L. von Albedyll, H. Griesche, and G. Spreen: Asymmetries in cloud microphysical properties ascribed to sea ice leads via water vapour transport in the central Arctic, Atmos. Chem. Phys., DOI:10.5194/acp-23-14521-2023, 2023.

[6] Saavedra Garfias, P., H. Kalesse-Los:Wintertime Arctic cloud properties coupled to sea ice leads during MOSAiC expedition, Harvard Dataverse, V1, DOI:10.7910/DVN/DZSUV7, 2023.

[7] Saavedra Garfias, P. et al., "Long-term statistical analysis of wintertime cloud thermodynamic phase and micro-physical properties in relation to sea ice condition at NSA Utqiaǵvik site". American Geophysical Union 2023, 11-15 December 2023, San Francisco, USA. Presentation [DOI](https://essopenarchive.org/users/552406)

[8] Saavedra Garfias, P. et al., "Variations of cloud properties ascribed by sea ice states in the Central and Western Arctic". XXVIII General Assembly of the International Union of Geodesy and Geophysics (IUGG) 2023, 12-19 July 2023, Berlin, Germany. Abstract DOI:10.57757/IUGG23-3160, Presentation DOI:10.22541/essoar.169008271.12504472/v1

[9] Saavedra Garfias, P. et al., "Climatology of clouds containing supercooled liquid in the Western and Central Arctic". American Geophysical Union 2021, 13-17 December 2021, New Orleans, USA. DOI10.1002/essoar.10509918.1

---
--<br>
(c) 2021, Pablo Saavedra Garfias<br>
pablo.saavedra@uni-leipzig.de

University of Leipzig<br>
Faculty of Physics and Geosciences<br>
LIM<br>

See LICENSE

