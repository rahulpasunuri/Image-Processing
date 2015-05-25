make
#./p1 2.1 lena.jpg greyLena.jpg
#./p1 2.2a greyLena.jpg bwLena_a.jpg
#./p1 2.2b greyLena.jpg bwLena_b.jpg
#./p1 3.1 greyLena.jpg uniformNoiseLena.jpg
#./p1 3.2 greyLena.jpg gaussianNoiseLena.jpg 0.5
#./p1 3.3 greyLena.jpg salAndPepperNoiseLena.jpg
#valgrind --leak-check=full ./p1 4.2a greyLena.jpg meanFilterLena.jpg 3
#./p1 4.2a greyLena.jpg meanFilterLena2.jpg 21

#./p1 4.2b greyLena.jpg gaussianFilterLena2.jpg 4.0
#./p1 4.3 greyLena.jpg  medianFilterLena.jpg 101
#valgrind --leak-check=full ./p1 4.3 salAndPepperNoiseLena.jpg medianFilterLena.jpg 11
#./p1 5 lena.jpg greyLena.jpg


####################
#MSE 5:- 
#rm -rf MSE/*.jpg
#./p1 2.1 lena.jpg MSE/m_greyLena.jpg
#./p1 3.1 MSE/m_greyLena.jpg MSE/m_uniformNoiseLena.jpg
#./p1 3.2 MSE/m_greyLena.jpg MSE/m_gaussianNoiseLena.jpg 0.5
#./p1 3.3 MSE/m_greyLena.jpg MSE/m_salAndPepperNoiseLena.jpg

#./p1 4.2a MSE/m_uniformNoiseLena.jpg MSE/m_uniform_meanFilterLena.jpg 21
#./p1 4.2a MSE/m_gaussianNoiseLena.jpg MSE/m_gauss_meanFilterLena.jpg 21
#./p1 4.2a MSE/m_salAndPepperNoiseLena.jpg MSE/m_salt_meanFilterLena.jpg 21

#./p1 4.2b MSE/m_uniformNoiseLena.jpg MSE/m_uniform_gaussianFilterLena.jpg 1
#./p1 4.2b MSE/m_gaussianNoiseLena.jpg MSE/m_gauss_gaussianFilterLena.jpg 1
#./p1 4.2b MSE/m_salAndPepperNoiseLena.jpg MSE/m_salt_gaussianFilterLena.jpg 1

#./p1 4.3 MSE/m_uniformNoiseLena.jpg MSE/m_uniform_medianFilterLena.jpg 101
#./p1 4.3 MSE/m_gaussianNoiseLena.jpg MSE/m_gauss_medianFilterLena.jpg 101
#./p1 4.3 MSE/m_salAndPepperNoiseLena.jpg MSE/m_salt_medianFilterLena.jpg 101

#./p1 5 MSE/m_greyLena.jpg MSE/m_uniform_meanFilterLena.jpg
#./p1 5 MSE/m_greyLena.jpg MSE/m_gauss_meanFilterLena.jpg
#./p1 5 MSE/m_greyLena.jpg MSE/m_salt_meanFilterLena.jpg

#./p1 5 MSE/m_greyLena.jpg MSE/m_uniform_gaussianFilterLena.jpg
#./p1 5 MSE/m_greyLena.jpg MSE/m_gauss_gaussianFilterLena.jpg
#./p1 5 MSE/m_greyLena.jpg MSE/m_salt_gaussianFilterLena.jpg

#./p1 5 MSE/m_greyLena.jpg MSE/m_uniform_medianFilterLena.jpg
#./p1 5 MSE/m_greyLena.jpg MSE/m_gauss_medianFilterLena.jpg
#./p1 5 MSE/m_greyLena.jpg MSE/m_salt_medianFilterLena.jpg

#echo "mse_mean_unif" >> "report.txt"
################################

#6.1
#./p1 6.1 greyLena.jpg meanFilterLena3.jpg 21 0
#./p1 6.1 greyLena.jpg gaussianFilterLena4.jpg 21 1
#./p1 6.2 greyLena.jpg test.jpg 15
#./p1 6.3 greyLena.jpg test.jpg 1
#./p1 6.3 greyLena.jpg fastGaussian.jpg 3
#./p1 4.2b greyLena.jpg simpleGaussian.jpg 3
#./p1 6.3 goku.jpg gokuGaussian.jpg 7

#./p1 2.1 beck.jpeg beckg.jpeg
./p1 2.1 cosmo1.jpg cosmo1g.jpg
./p1 2.1 cosmo2.jpg cosmo2g.jpg
./p1 2.1 tattoo.jpg tattoog.jpg

