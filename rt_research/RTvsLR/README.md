# LR vs RT vs xTB



In this project we will compare spectrum obtained from three different methods; LR-TDDFT, RT-TDDFT, and xtb-TD



Several small dye molecules with a wide range of optical response have been choosen to explore effect of different variables on the final spectrum 



Dye list: benzothiazole, indigo, DL-Tryptophan, ...

- Add some other small dyes with different optical characters (\pi-\pi, D-A, n-\sigma, ...excitations) to above list
- Use 3 different categories of xc to calculate these optical responses (PBE as pure, B3LYP as Hybrid, CAM-B3LYP as range Seperated )

### Question to Answer

- What is the difference between first peak of spectrum corresponding to different XCs ( We can expect smaller difference in the case of RT-TD because of dipole key role in their response?) 

 

## To Do

- Calculate all optical responses LR using g09 or NWchem (N dye and 3 XC )
- Perform all RT-TDDFT calculations employing NWchem  (It is important to find proper dt, t, and padding and broadening related parameter  ...It would be tricky a bit)
- Calculate all above using x-tb  
- categorize dyes into D-A organic, Amino, ....families and argue in each case and find possible correlation between dye family and above comparisons ..





