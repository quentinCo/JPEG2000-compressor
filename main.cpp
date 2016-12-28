#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>

#include "./functions.c"
#include "./quantlm.cpp"

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

	return error / p;
}

double psnr(const vector<double>& x , const vector<double>& y)
{
	return (10 * log10(65025 / error(x,y))); // 65 025 = 255^2 
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
	double sum = std::accumulate(x.begin(), x.end(), 0);
	return ( sum / static_cast<double>(x.size()));
}

double computeVariance(const std::vector<double> x)
{
	double average = computeAverage(x);
	double sum = std::accumulate(x.begin(), x.end(), 0,
		[average](double sum, double b) {
		return sum + (b - average) * (b - average);
	});
	return (sum / static_cast<double>(x.size()));
}

void coutValues(const std::vector<double>& x)
{
	std::cout << "Size: " << x.size() << std::endl;
	std::cout << "Min: " << *std::min_element(x.begin(), x.end()) << std::endl;
	std::cout << "Max: " << *std::max_element(x.begin(), x.end()) << std::endl;
	double average = computeAverage(x);
	std::cout << "Average: " << average << std::endl;
	double variance = computeVariance(x);
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
		analyse_97(line);
		setLine(image, line, width, i);
	}

	std::cout << "\tColonne Processing" << std::endl;
	for(size_t i = 0; i < width; i++)
	{
		std::vector<double> colonne = getColonne(image, width, height, i);
		analyse_97(colonne);
		setColonne(image, colonne, width, height, i);
	}
}

void synthese2D_97(vector<double>& image, size_t width, size_t height)
{
	for(size_t i = 0; i < width; i++)
	{
		std::vector<double> colonne = getColonne(image, width, height, i);
		synthese_97(colonne);
		setColonne(image, colonne, width, height, i);
	}

	for(size_t i = 0; i < height; i++)
	{
		std::vector<double> line = getLine(image, width, i);
		synthese_97(line);
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

/* TD 4 */
vector<double> getSubPicture(const vector<double>& image, size_t width, size_t subWidth, size_t subHeight, size_t indexWidth, size_t indexHeight)
{
	vector<double> subImage;
	for(size_t i = indexHeight; i < (indexHeight + subHeight); ++i)
		for(size_t j = indexWidth; j < (indexWidth + subWidth); ++j)
			subImage.push_back(image[i * width + j]);

	return subImage;
}

void setSubPicture(vector<double>& image, vector<double>& subPicture, size_t width, size_t subWidth, size_t subHeight, size_t indexWidth, size_t indexHeight)
{
	for(size_t i = 0; i < subHeight; ++i)
		for(size_t j = 0; j < subWidth; ++j)
			image[(indexHeight + i) * width + (indexWidth + j)] = subPicture[i * subWidth + j];
}

vector<double> subband2D(const vector<double>& image, size_t width, size_t height, int level)
{
	std::vector<double> variances;
	if(level > 0)
	{
		size_t halfWidth = width / 2;
		size_t halfHeight = height / 2;

		vector<double> imageA = getSubPicture(image, width, halfWidth, halfHeight, 0, 0);
		variances = subband2D(imageA, halfWidth, halfHeight, level-1);

		std::cout << "\nLevel: " << level << std::endl;
		vector<double> imageDH = getSubPicture(image, width, halfWidth, halfHeight, 0, halfHeight);
		std::cout << "\tImage DH: " << std::endl;
		coutValues(imageDH);
		variances.push_back(computeVariance(imageDH));

		vector<double> imageDV = getSubPicture(image, width, halfWidth, halfHeight, halfWidth, 0);
		std::cout << "\tImage DV: " << std::endl;
		coutValues(imageDV);
		variances.push_back(computeVariance(imageDV));

		vector<double> imageDD = getSubPicture(image, width, halfWidth, halfHeight, halfWidth, halfHeight);
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

vector<double> debitBand(const std::vector<double>& variances, float debit, int level, size_t imgSize)
{
	vector<double> debitsPerBands;
	double product = 1;
	for(size_t i = 0; i < variances.size(); ++i)
	{
		product *= pow(variances[i],double(1 / pow(4,level)));
		std::cout << "pow : " << pow(variances[i],double(1 / pow(2,level+1))) << std::endl;
		std::cout << "Nj/N ["<< i << "] : " << pow(2,level+1) << std::endl;
		std::cout << "variance ["<< i << "] : " << variances[i] << std::endl;
		std::cout << "product ["<< i << "] : " << product << std::endl;
		
		if((i > 0) && ((i % 3) == 0))
		{
			level -= 1;
		}
	}
	std::cout << "Prod: " << product << std::endl;
	for(auto variance : variances)
	{
    	double dpb = debit + 0.5 * log2(variance /  product);
		debitsPerBands.push_back(dpb);
		std::cout << "Debit per band : " << dpb << std::endl;
	}
	return debitsPerBands;
}

/* TD 5 */
void quantifier(vector<double>& image, const vector<double>& debits, size_t width, size_t height, int level)
{
	std::vector<double> variances;
	if(level > 0)
	{
		int indexDebit = 3 * (level - 1);

		size_t halfWidth = width / 2;
		size_t halfHeight = height / 2;

		vector<double> imageA = getSubPicture(image, width, halfWidth, halfHeight, 0, 0);
		quantifier(imageA, debits, halfWidth, halfHeight, level-1);
		setSubPicture(image, imageA, width, halfWidth, halfHeight, 0, 0);		

		vector<double> imageDH = getSubPicture(image, width, halfWidth, halfHeight, 0, halfHeight); //TODO : DV
		quantlm(imageDH.data(), imageDH.size(), ceil(pow(2,debits[indexDebit + 1])));
		setSubPicture(image, imageDH, width, halfWidth, halfHeight, 0, halfHeight);	


		vector<double> imageDV = getSubPicture(image, width, halfWidth, halfHeight, halfWidth, 0);
		quantlm(imageDV.data(), imageDV.size(), ceil(pow(2,debits[indexDebit + 2])));
		setSubPicture(image, imageDV, width, halfWidth, halfHeight, halfWidth, 0);	

		vector<double> imageDD = getSubPicture(image, width, halfWidth, halfHeight, halfWidth, halfHeight);
		quantlm(imageDD.data(), imageDD.size(), ceil(pow(2,debits[indexDebit + 3])));
		setSubPicture(image, imageDD, width, halfWidth, halfHeight, halfWidth, halfHeight);
	}
	else
	{
		quantlm(image.data(), image.size(), ceil(pow(2,debits[0])));
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
	std::cout << "" << std::endl;
	vector<double> debitsPerBand = debitBand(variances, 1, 3, image.size());

	std::cout << "Quantification" << std::endl;
	quantifier(image, debitsPerBand, dim, dim, 3);

	for(size_t i = 0; i < dim*dim ; ++i)
		data[i] = image[i];

	exitPath = "./amr2D_97_quantification_lena.bmp";
	ecrit_bmp256(exitPath.c_str(), dim, dim, data);

	iamr2D_97(image, dim, dim, 3);
	for(size_t i = 0; i < dim*dim ; ++i)
		data[i] = image[i];

	exitPath = "./iamr2D_97_lena.bmp";
	ecrit_bmp256(exitPath.c_str(), dim, dim, data);	

	std::cout << "PSNR: " << psnr(imageOriginal, image) << std::endl;
}

void executeCompression(const string& methodName, const string& signalName, vector<double>& signal, void (*analyse)(vector<double>&), void (*synthese)(vector<double>&))
{
	string fileName = "";
	/* Analyse */
	vector<double> processedSignal = signal;
	cout <<"\tAnalyse de " << methodName << endl;
	analyse(processedSignal);
	fileName = "./ouput_data/" + signalName +"_analyse_" + methodName + ".txt";
	save_signal(processedSignal, fileName.c_str());

	/* Reconstitution */
	cout <<"\tReconstitution de " << methodName << endl;
	synthese(processedSignal);
	fileName = "./ouput_data/" + signalName +"_synthese_" + methodName + ".txt";
	save_signal(processedSignal, fileName.c_str());

	cout << "\t\tReconstitution Error : " << error(processedSignal, signal) <<endl;
}

void testFiltringProcess(vector<double>& signal, const string& signalName)
{
	/* Haar */
	executeCompression("haar", signalName, signal, &analyse_haar, &synthese_haar);

	/* biorthogonaux 9/7 */
	executeCompression("biorthogonaux_9_7", signalName, signal, &analyse_97, &synthese_97);

	/* Lifting de biorthogonaux 9/7 */
	if(signalName != "ramp")
		executeCompression("biorthogonaux_9_7", signalName, signal, &analyse_97_lifting, &synthese_97_lifting);
}


/* Main */
int main (int argc, char* argv[])
{
/* Rampe */
	cout << "Rampe Signal" << std::endl;
	vector<double> rampe;
	for(int i = 0; i < 255; i++)
		rampe.push_back(i);
	save_signal(rampe,"rampe.txt");

	testFiltringProcess(rampe, "rampe");

/* Leleccum */
	cout << "Leleccum" << std::endl;
	vector<double> leleccum;
	read_signal(leleccum, 4096, "./leleccum.txt");

	testFiltringProcess(leleccum, "leleccum");

/* Test */
	cout << "Signal Test" << std::endl;
	vector<double> test;
	read_signal(test, 512, "./test.txt");

	testFiltringProcess(test, "signal_test_file");


/* AMR */
	vector<double> processedSignal = test;
	int level = 2;	//log2(test.size());

	cout <<"\n AMR synthese de lifting 9/7" << endl;
	amr(processedSignal, level);
	save_signal(processedSignal,"./ouput_data/test_amr_lifting.txt");

	subband(processedSignal, level);

	iamr(processedSignal, level);
	save_signal(processedSignal,"./ouput_data/test_iamr_lifting.txt");
	cout << "Error: " << error(processedSignal, test) <<endl;


/* 2D */
	std::cout << "\nImage Processing" << std::endl;
	image_processing();

	return 1;
}