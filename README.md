# The computational package for calculating CSS (Circadian Sleep Sufficiency)
This repository contains the matlab code for calculating CSS (Circadian Sleep Sufficiency) via calculating Homeostatic sleep pressure and circadian sleep by tracked sleep patterns.
## The description of files in the package
1. CSS_package.m
> The main file of this package. It calculates CSS and circadian necessary sleep for each sleep episodes of users by using the sleep patterns and the light profiles of users provided as Input1_sleep_light.csv and Input2_WASO_Main.csv. 
2. Input1_sleep_light.csv
> This file is the input of CSS_package.m which need to consist of the sleep pattern of user as the 1st column and the light profile of user as the 2nd column. See [the example for details](Input1_sleep_light2.csv) 
3. Input2_WASO_Main.csv
> This file is the input of CSS_package.m which need to consist of the WASO of each sleep episode and whether each sleep epsiode is an main sleep as the 1st and the 2nd column, respectively. See [the example for details](Input1_WASO_main2.csv)  
4. Sleep_make.m
> This function create the Maked_sleep.csv which can be used as 1st column of Input1_sleep_light.csv. Users need to fill in date of sleep onset (year-month-day), sleep onset (h), date of sleep offset and sleep offset in the 1st, 2nd, 3rd and 4th column of ‘Input1_sleep.csv’, respectively which is the input of this function. See [the example for details](Input1_sleep.csv) 
5. WASO_make.m
>

