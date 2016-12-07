#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>

#include "./functions.c"

using namespace std;


void print_signal(const vector<double>& x)
{
	for (size_t i=0;i<x.size();i++)
		cout << "x[" << i << "]= " << x[i] << endl;
	cout << endl;
}

void save_signal(const vector<double>& x, const string& filename) {
	FILE *fp=fopen(filename.c_str(),"wt");
	for (size_t i=0;i<x.size();i++) {
	fprintf(fp,"%f\n",x[i]);
	}
	fclose(fp);
}

void read_signal(vector<double>& x, int size, const string& filename) {
	FILE *fp=fopen(filename.c_str(),"rt");
	x = vector<double>(size);
	for (int i=0;i<size;i++) {
		float t;
		fscanf(fp,"%f\n",&t);
		x[i] = static_cast<double>(t);
	}
	fclose(fp);
}

/* TP1 */
/* 1 - Interpolation d'un facteur 2 */
void interpolation_2(vector<double>& x)
{
	int p = x.size();
	vector<double> y (p);
	for(int i = 0; i < p / 2; i++)
	{
		y[2 * i] = x[i];
		y[2 * i +1] = 0; 
	}

	x = y;
}

/* 2 - Décimation d'un facteur 2 */
void decimation_2(vector<double>& x)
{
	int p = x.size();
	vector<double> y(p); // y.size() doit être égale à x.size() (concervation de la taille)
	fill(y.begin(), y.end(), 0);
	for(int i = 0; i < p / 2; i++)
		y[i] = x[2 * i];
	
	x = y;
}

/* 3 - Filtrage */
void filtre(vector<double>& x, const vector<double>& _h)
{
	int p = x.size();
	int q = _h.size();
	vector<double> y(p);
	for(int i = 0; i < p; i++)
	{
		double sum = 0;
		for(int j = 0; j < q; j++)
		{
			int k = j - q / 2;
			int n = i - k;
			double xValue;
			if(n < 0 || n >= p)
				xValue = x[i + k];
			else
				xValue = x[n];

			sum += _h[j] * xValue;
		}
		y[i] = sum;
	}

	x = y;
}

/* 4 - Banc de filtres d'analyse de Haar */
void half_analyse(const vector<double>& x, const vector<double>& _h, vector<double>& res)
{
	res = x;
	filtre(res, _h);
	decimation_2(res);
}

void analyse(vector<double>& x, const vector<double>& _h0, const vector<double>& _h1)
{
	int p = x.size();
	vector<double> xbd = x;
	half_analyse(x, _h0, xbd);

	vector<double> xhd = x;
	half_analyse(x, _h1, xhd);

	// Concatenation
	for(int i = 0; i <= p/2; i++)
	{
		x[i] = xbd[i];
		if(i + p/2 < p)
			x[i + p/2] = xhd[i];
	}
}

void analyse_haar(vector<double>& x)
{
	vector<double> _h0 = {1/sqrt(2), 1/sqrt(2), 0}; // passe bas
	vector<double> _h1 = {1/sqrt(2), -1/sqrt(2), 0};	// passe haut
	analyse(x, _h0, _h1);
}

/* 5 - Reconstruction */
void half_synthese(vector<double>& x, const vector<double>& _h)
{
	interpolation_2(x);
	filtre(x, _h);
}

void synthese(vector<double>& x, const vector<double>& _g0, const vector<double>& _g1)
{
	int p = x.size();
	vector<double> yb(p);
	vector<double> yh(p);

	// Init ybd et yhd
	for(int i = 0; i <= p/2; i++)
	{
		yb[i] = x[i];
 
		int halfIndex = i + p/2;
		if(halfIndex < p)
		{
			yh[i] = x[halfIndex];
			yh[halfIndex] = 0;
			yb[halfIndex] = 0;
		}
	}
	
	// Synthese
	half_synthese(yb, _g0);
	half_synthese(yh, _g1);
	
	for(int i = 0; i < p; i++)
		x[i] = yb[i] + yh[i];
}

void synthese_haar(vector<double>& x)
{
	vector<double> _g0 = {0, 1/sqrt(2), 1/sqrt(2)}; // passe bas
	vector<double> _g1 = {0, -1/sqrt(2), 1/sqrt(2)};	// passe haut
	synthese(x, _g0, _g1);
}

/* 6 - Banc de filtres biorthogonaux 9/7 */
void analyse_97(vector<double>& x)
{
	vector<double> _h0 = {0.037828455507, -0.023849465019, -0.110624404418, 0.377402855613,
	0.852698679009, 0.377402855613, -0.110624404418, -0.023849465019, 0.037828455507}; // passe bas
	vector<double> _h1 = {0.064538882629, -0.040689417610, -0.418092273222, 0.788485616406,
	-0.418092273222, -0.040689417610, 0.064538882629, 0, -0};	// passe haut
	
	analyse(x, _h0, _h1);
}

void synthese_97(vector<double>& x)
{
	vector<double> _g0 = {-0.064538882629, -0.040689417610, 0.418092273222, 0.788485616406,
	0.418092273222, -0.040689417610, -0.064538882629}; // passe bas
	vector<double> _g1 = {0, -0, 0.037828455507, 0.023849465019, -0.110624404418, -0.377402855613,
	0.852698679009, -0.377402855613, -0.110624404418, 0.023849465019, 0.037828455507};	// passe haut
	
	synthese(x, _g0, _g1);
}

/* TP2 */
int mirror_increment(int index, int increment, int size) 
{
	/*
	if((index + increment) < 0)
		return -(index + increment);
	else if((index + increment) >= size)
		return (size - index) + (size - (index + increment + 1));
	else
		return increment;
	*/
	return ((index + increment) < 0 || (index + increment) >= size) ? (-increment) : increment;
}

void analyse_97_lifting_prediction(vector<double>& x, double a)
{
	size_t p = x.size();
	size_t limit = p * 0.5;
	for(size_t i = 0; i < limit; ++i)
	{
		x[2 * i + 1] =  a * x[2 * i] + x[2 * i + mirror_increment(2 * i, 1, p)] + a * x[2 * i + mirror_increment(2 * i, 2, p)];
	}
}

void analyse_97_lifting_update(vector<double>& x, double a)
{
	size_t p = x.size();
	size_t limit = p * 0.5;
	for(size_t i = 0; i < limit; ++i)
	{
		x[2 * i] = a * x[2 * i + mirror_increment(2 * i, -1, p) ]+ x[2 * i] + a * x[2 * i + mirror_increment(2 * i, 1, p)];
	}
}

/* 1 - Décomposition 9/7 par lifting */
void analyse_97_lifting(vector<double>& x)
{
	size_t p = x.size();

	// Prediction 1
	analyse_97_lifting_prediction(x, -1.586134342);

	// Update 1
	analyse_97_lifting_update(x, -0.05298011854);

	// Prediction 2
	analyse_97_lifting_prediction(x, 0.8829110762);

	// Update 2	
	analyse_97_lifting_update(x, 0.4435068522);

	// Scaling
	float a = 1 / 1.149604398;
	for(size_t i = 0; i < p; i += 2)
	{
		x[i] /= a;
		x[i + 1] *= a;
	}

	// 
	vector<double> y(p);
	size_t halfP = p * 0.5;
	for(size_t i = 0; i < halfP; ++i)
	{
		y[i] = x[2 * i];
		y[i + halfP] = x[2 * i + 1];
	}
	x = y;
}


/* Reconstruction par lifting */
void synthese_97_lifting(vector<double>& x)
{
	size_t p = x.size();

	{
		vector<double> y(p);
		size_t halfP = p * 0.5;
		for(size_t i = 0; i < halfP; ++i)
		{
			y[2 * i] = x[i];
			y[2 * i + 1] = x[i + halfP];
		}
		x = y;
	}

	// Scaling
	float a = 1.149604398;
	for(size_t i = 0; i < p; i += 2)
	{
		x[i] /= a;
		x[i + 1] *= a;
	}

	// Update 2
	analyse_97_lifting_update(x, -0.4435068522);

	// Prediction 2
	analyse_97_lifting_prediction(x, -0.8829110762);

	// Update 1
	analyse_97_lifting_update(x, 0.05298011854);

	// Prediction 1
	analyse_97_lifting_prediction(x, 1.586134342);
}

/* Error */
double error(const vector<double>& x , const vector<double>& y)
{
	size_t p = x.size();

	double error = 0;
	for(size_t i = 0; i < p; ++i)
	{
		error += (x[i] - y[i]) * (x[i] - y[i]); 
	}

	return error;
}

/* TP 3 */
	/* AMR */
void amr(std::vector<double>& x, int level)
{
	if(level > 0)
	{
		analyse_97_lifting(x);
		
		// Split
		vector<double> xa;
		xa.insert(xa.begin(), x.begin(), x.begin() + (0.5 * x.size()));  
		vector<double> xd;
		xd.insert(xd.begin(), x.begin() + (0.5 * x.size()), x.end());

		amr(xa, level-1);
		
		// Merge		
		xa.insert( xa.end(), xd.begin(), xd.end() );
		x = xa;
	}
}

void iamr(std::vector<double>& x, int level)
{
	if(level > 0)
	{
		// Split
		vector<double> xa;
		xa.insert(xa.begin(), x.begin(), x.begin() + (0.5 * x.size()));  
		vector<double> xd;
		xd.insert(xd.begin(), x.begin() + (0.5 * x.size()), x.end());
		
		iamr(xa, level-1);

		// Merge
		xa.insert( xa.end(), xd.begin(), xd.end() );
		synthese_97_lifting(xa);
		x = xa;
	}
}

double computeAverage(const std::vector<double>& x)
{
	return std::accumulate(x.begin(), x.end(), 0) / static_cast<double>(x.size());
}

double computeVariance(const std::vector<double> x)
{
	double average = computeAverage(x);
	double sum = std::accumulate(x.begin(), x.end(), 0,
		[average](double sum, double b) {
		return sum + (b - average) * (b - average);
	});
	return sum / static_cast<double>(x.size());
}

void coutValues(const std::vector<double>& x)
{
	std::cout << "Size: " << x.size() << std::endl;
	std::cout << "Min: " << *std::min_element(x.begin(), x.end()) << std::endl;
	std::cout << "Max: " << *std::max_element(x.begin(), x.end()) << std::endl;
	double average = computeAverage(x);
	std::cout << "Average: " << average << std::endl;
	double variance = computeVariance(x);
	/*double varianceSum = 0;
	for(auto it : x)
		varianceSum += (it - average)*(it - average);
	*/
	std::cout << "Variance: " << variance << std::endl;
}

void subband(std::vector<double>& x, int level)
{
	if(level > 0)
	{
		vector<double> xa;
		xa.insert(xa.begin(), x.begin(), x.begin() + (0.5 * x.size()));  
		vector<double> xd;
		xd.insert(xd.begin(), x.begin() + (0.5 * x.size()), x.end());
		std::cout << "xd" << std::endl;
		coutValues(xd);
		subband(xa, level-1);
	}
	else
	{
		std::cout << "xa" << std::endl;
		coutValues(x);
	}
}

	/* 2D */
std::vector<double> getLine(const std::vector<double>& image, size_t width, size_t index)
{
	std::vector<double> line;
	/*for(size_t i = 0; i < width; ++i)
		line.push_back(image[index * width + i]);*/
	line.insert(line.begin(), image.begin() + index * width, image.begin() + (index + 1) * width /*- 1*/);
	return line;
}

void setLine(std::vector<double>& image, const std::vector<double>& line, size_t width, size_t index)
{
	for(size_t i = 0; i < width; ++i)
		image[index * width + i] = line[i];
}

std::vector<double> getColonne(const std::vector<double>& image, size_t width, size_t height, size_t index)
{
	std::vector<double> colonne;
	for(size_t i = 0; i < height; ++i)
		colonne.push_back(image[i*width + index]);

	return colonne;
}

void setColonne(std::vector<double>& image, const std::vector<double>& colonne, size_t width, size_t height, size_t index)
{
	for(size_t i = 0; i < height; ++i)
		image[i* width + index] = colonne[i];
}

void analyse2D_97(vector<double>& image, size_t width, size_t height)
{
	std::cout << "\tLine Processing" << std::endl;
	for(size_t i = 0; i < height; i++)
	{
		std::vector<double> line = getLine(image, width, i);
		analyse_97_lifting(line);
		setLine(image, line, width, i);
	}

	std::cout << "\tColonne Processing" << std::endl;
	for(size_t i = 0; i < width; i++)
	{
		std::vector<double> colonne = getColonne(image, width, height, i);
		analyse_97_lifting(colonne);
		setColonne(image, colonne, width, height, i);
	}
}

void synthese2D_97(vector<double>& image, size_t width, size_t height)
{
	for(size_t i = 0; i < width; i++)
	{
		std::vector<double> colonne = getColonne(image, width, height, i);
		synthese_97_lifting(colonne);
		setColonne(image, colonne, width, height, i);
	}

	for(size_t i = 0; i < height; i++)
	{
		std::vector<double> line = getLine(image, width, i);
		synthese_97_lifting(line);
		setLine(image, line, width, i);
	}
}

void amr2D_97(vector<double>& image, size_t width, size_t height, int level)
{
	if(level > 0)
	{
		analyse2D_97(image, width, height);
		
		// Split
		vector<double> imageA;
		size_t widthA = 0.5 * width;
		size_t heightA = 0.5 * height;
		for(size_t i = 0; i < heightA; ++i)
		{
			for(size_t j = 0; j < widthA; ++j)
				imageA.push_back(image[i * width + j]);
		}
		amr2D_97(imageA, widthA, heightA, level-1);
		
		// Merge		
		for(size_t i = 0; i < heightA; ++i)
		{
			for(size_t j = 0; j < widthA; ++j)
				image[i * width + j] = imageA[i * widthA + j];
		}
	}
}

void iamr2D_97(vector<double>& image, size_t width, size_t height, int level)
{
	if(level > 0)
	{
		vector<double> imageA;
		size_t widthA = 0.5 * width;
		size_t heightA = 0.5 * height;
		for(size_t i = 0; i < heightA; ++i)
		{
			for(size_t j = 0; j < widthA; ++j)
				imageA.push_back(image[i * width + j]);
		}

		iamr2D_97(imageA, widthA, heightA, level-1);

				// Merge		
		for(size_t i = 0; i < heightA; ++i)
		{
			for(size_t j = 0; j < widthA; ++j)
				image[i * width + j] = imageA[i * widthA + j];
		}

		synthese2D_97(image, width, height);
	}
}

vector<double> subPicture(const vector<double>& image, size_t width, size_t subWidth, size_t subHeight, size_t indexWidth, size_t indexHeight)
{
	//std::cout << "size: " << image.size() << " -- width: " << width << " -- subWidth: " << subWidth << " -- subHeight: " << subHeight << " -- indexWidth: " << indexWidth << " -- indexHeight: " << indexHeight << std::endl;
	vector<double> subImage;
	for(size_t i = indexHeight; i < (indexHeight + subHeight); ++i)
		for(size_t j = indexWidth; j < (indexWidth + subWidth); ++j)
			subImage.push_back(image[i * width + j]);

	//std::cout << "subImage.size(): " << subImage.size() << std::endl;

	return subImage;
}

vector<double> subband2D(const vector<double>& image, size_t width, size_t height, int level)
{
	std::vector<double> variances;
	if(level > 0)
	{
		size_t halfWidth = width / 2;
		size_t halfHeight = height / 2;

		vector<double> imageA = subPicture(image, width, halfWidth, halfHeight, 0, 0);
		variances = subband2D(imageA, halfWidth, halfHeight, level-1);

		std::cout << "\nLevel: " << level << std::endl;
		vector<double> imageDH = subPicture(image, width, halfWidth, halfHeight, 0, halfHeight);
		std::cout << "\tImage DH: " << std::endl;
		coutValues(imageDH);
		variances.push_back(computeVariance(imageDH));

		vector<double> imageDV = subPicture(image, width, halfWidth, halfHeight, halfWidth, 0);
		std::cout << "\tImage DV: " << std::endl;
		coutValues(imageDV);
		variances.push_back(computeVariance(imageDV));

		vector<double> imageDD = subPicture(image, width, halfWidth, halfHeight, halfWidth, halfHeight);
		std::cout << "\tImage DD: " << std::endl;
		coutValues(imageDD);
		variances.push_back(computeVariance(imageDD));
	}
	else
	{
		std::cout << "\nLevel: " << level << std::endl;
		std::cout << "\tImage A: " << std::endl;
		coutValues(image);

		variances.push_back(computeVariance(image));
	}
	return variances;
}

void debitBand(const std::vector<double>& variances, float debit, int level)
{
	double product = 0;
	for(size_t i = 0; i < variances.size(); ++i)
	{
		/*int index = (i == 0)? 0 : i - 1;
		int*/ 
	}
}

void image_processing()
{
	string filePath = "./lena.bmp";
	uint32_t dim = 512;

	vector<double> imageOriginal;

	double* data = charge_bmp256(filePath.c_str(), &dim, &dim);
	for(size_t i = 0; i < dim * dim; ++i)
		imageOriginal.push_back(data[i]);

	vector<double> image = imageOriginal;
	analyse2D_97(image, dim, dim);

	for(size_t i = 0; i < dim*dim ; ++i)
		data[i] = image[i];

	string exitPath = "./analyse_lifting_lena.bmp";
	ecrit_bmp256(exitPath.c_str(), dim, dim, data);

	synthese2D_97(image, dim, dim);

	for(size_t i = 0; i < dim*dim ; ++i)
		data[i] = image[i];

	exitPath = "./synthese_lifting_lena.bmp";
	ecrit_bmp256(exitPath.c_str(), dim, dim, data);	


	/* AMR */
	image = imageOriginal;
	amr2D_97(image, dim, dim, 3);
	for(size_t i = 0; i < dim*dim ; ++i)
		data[i] = image[i];

	exitPath = "./amr2D_97_lena.bmp";
	ecrit_bmp256(exitPath.c_str(), dim, dim, data);

	std::cout << "Subband 2D" << std::endl;
	vector<double> variances = subband2D(image, dim, dim, 3);

	/* DEBIT */


	iamr2D_97(image, dim, dim, 3);
	for(size_t i = 0; i < dim*dim ; ++i)
		data[i] = image[i];

	exitPath = "./iamr2D_97_lena.bmp";
	ecrit_bmp256(exitPath.c_str(), dim, dim, data);	
}

/* Main */
int main (int argc, char* argv[])
{
/* Rampe */
	cout << "Rampe" << std::endl;
	vector<double> rampeOriginal;
	for(int i = 0; i < 255; i++)
		rampeOriginal.push_back(i);
	save_signal(rampeOriginal,"rampe.txt");

	/* Haar */
		/* Analyse de Haar */
	vector<double> rampe = rampeOriginal;
	cout <<"\n Analyse de Haar" << endl;
	analyse_haar(rampe);
	save_signal(rampe,"./ouput_data/rampe_analyse_haar.txt");

		/* Reconstitution de Haar */
	cout <<"\n Reconstitution de Haar" << endl;
	synthese_haar(rampe);
	save_signal(rampe,"./ouput_data/rampe_synthese_haar.txt");
	/*
		Reconstitution parfaite.
	*/

	/* biorthogonaux 9/7 */
	rampe = rampeOriginal;
		/* Analyse biorthogonaux 9/7 */
	cout <<"\n Analyse de biorthogonaux 9/7" << endl;
	analyse_97(rampe);
	save_signal(rampe,"./ouput_data/rampe_analyse_97.txt");

		/* Reconstition biorthogonaux 9/7 */
	cout <<"\n Synthese de biorthogonaux 9/7" << endl;
	synthese_97(rampe);
	save_signal(rampe,"./ouput_data/rampe_synthese_97.txt");
	/*
		Reconstitution proche du signal. Différence sur les dernières valeurs. 
		Les valeurs reconstituées sont différentes des valeurs d'origines, probablement dû 
		à une approximation utilisée dans les valeurs des filtres passe bas et passe haut.
	*/

/* Leleccum */
	cout << "Leleccum" << std::endl;
	vector<double> leleccumOriginal;
	read_signal(leleccumOriginal, 4096, "./leleccum.txt");

	vector<double> leleccum = leleccumOriginal;

	/* Haar */
		/* Analyse de Haar */
	cout <<"\n Analyse de Haar" << endl;
	analyse_haar(leleccum);
	save_signal(leleccum,"./ouput_data/leleccum_analyse_haar.txt");

		/* Reconstitution de Haar */
	cout <<"\n Reconstitution de Haar" << endl;
	synthese_haar(leleccum);
	save_signal(leleccum,"./ouput_data/leleccum_synthese_haar.txt");
	cout << "Error: " << error(leleccum, leleccumOriginal) <<endl;

	/* biorthogonaux 9/7 */
	leleccum = leleccumOriginal;
		/* Analyse biorthogonaux 9/7 */
	cout <<"\n Analyse de biorthogonaux 9/7" << endl;
	analyse_97(leleccum);
	save_signal(leleccum,"./ouput_data/leleccum_analyse_97.txt");

		/* Reconstition biorthogonaux 9/7 */
	cout <<"\n Synthese de biorthogonaux 9/7" << endl;
	synthese_97(leleccum);
	save_signal(leleccum,"./ouput_data/leleccum_synthese_97.txt");
	cout << "Error: " << error(leleccum, leleccumOriginal) <<endl;

	/* Leleccum */
	leleccum = leleccumOriginal;

	cout <<"\n Lifting de biorthogonaux 9/7" << endl;
	analyse_97_lifting(leleccum);
	save_signal(leleccum,"./ouput_data/leleccum_analyse_97_lifting.txt");

	synthese_97_lifting(leleccum);
	save_signal(leleccum,"./ouput_data/leleccum_synthese_97_lifting.txt");
	cout << "Error: " << error(leleccum, leleccumOriginal) <<endl;

/* Test */
	cout << "Test" << std::endl;
	vector<double> testOriginal;
	read_signal(testOriginal, 512, "./test.txt");

	vector<double> test = testOriginal;

	/* Haar */
		/* Analyse de Haar */
	cout <<"\n Analyse de Haar" << endl;
	analyse_haar(test);
	save_signal(test,"./ouput_data/test_analyse_haar.txt");

		/* Reconstitution de Haar */
	cout <<"\n Reconstitution de Haar" << endl;
	synthese_haar(test);
	save_signal(test,"./ouput_data/test_synthese_haar.txt");
	cout << "Error: " << error(test, testOriginal) <<endl;

	/* biorthogonaux 9/7 */
		/* Analyse biorthogonaux 9/7 */
	test = testOriginal;
	cout <<"\n Analyse de biorthogonaux 9/7" << endl;
	analyse_97(test);
	save_signal(test,"./ouput_data/test_analyse_97.txt");

		/* Reconstition biorthogonaux 9/7 */
	cout <<"\n Synthese de biorthogonaux 9/7" << endl;
	synthese_97(test);
	save_signal(test,"./ouput_data/test_synthese_97.txt");
	cout << "Error: " << error(test, testOriginal) <<endl;

	/* Leleccum */
	test = testOriginal;

	cout <<"\n Lifting de biorthogonaux 9/7" << endl;
	analyse_97_lifting(test);
	save_signal(test,"./ouput_data/test_analyse_97_lifting.txt");

	synthese_97_lifting(test);
	save_signal(test,"./ouput_data/test_synthese_97_lifting.txt");
	cout << "Error: " << error(test, testOriginal) <<endl;


/* AMR */
	test = testOriginal;
	int level = 2;//log2(test.size());

	cout <<"\n AMR synthese de lifting 9/7" << endl;
	amr(test, level);
	save_signal(test,"./ouput_data/test_amr_lifting.txt");

	subband(test, level);

	iamr(test, level);
	save_signal(test,"./ouput_data/test_iamr_lifting.txt");
	cout << "Error: " << error(test, testOriginal) <<endl;


/* 2D */
	std::cout << "\nImage Processing" << std::endl;
	image_processing();

	return 1;
}