# cuzick
Perform the Cuzick's test on trend.<br/>
This function provides a Wilcoxon-type test for trend across a group of
three or more independent random samples.
Assumptions:
- Data must be at least ordinal
- Groups must be selected in a meaningful order i.e. ordered
If you do not choose to enter your own group scores then scores are
allocated uniformly (1 ... n) in order of selection of the n groups.
The null hypothesis of no trend across the groups T will have mean E(T),
variance var(T) and the null hypothesis is tested using the normalised
test statistic z.
A logistic distribution is assumed for errors. Please note that this test
is more powerful than the application of the Wilcoxon rank-sum /
Mann-Whitney test between more than two groups of data.
Cuzick J. A Wilcoxon-Type Test for Trend. Statistics in Medicine
1985;4:87-89.

           Created by Giuseppe Cardillo
           giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2008) Cuzick's test: A Wilcoxon-Type Test for Trend
http://www.mathworks.com/matlabcentral/fileexchange/22059
