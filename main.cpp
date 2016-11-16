#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;


void print_signal(const vector<double>& x)
{
	for (int i=0;i<x.size();i++)
		cout << "x[" << i << "]= " << x[i] << endl;
	cout << endl;
}

void save_signal(const vector<double>& x,char* filename) {
  int i;
  FILE *fp=fopen(filename,"wt");
  for (i=0;i<x.size();i++) {
   fprintf(fp,"%f\n",x[i]);
  }
  fclose(fp);
}

void read_signal(vector<double>& x, int size, char* filename) {
	FILE *fp=fopen(filename,"rt");
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
	return (index < 0 || index >= size) ? (-increment) : increment;
}

void analyse_97_lifting_prediction(vector<double>& x, float a)
{
	size_t p = x.size();
	size_t limit = p * 0.5;
	for(size_t i = 0; i < limit; ++i)
	{
		x[2 * i + 1] =  a * x[2 * i] + x[2 * i + mirror_increment(i, 1, p)] + a * x[2 * i + mirror_increment(i, 2, p)];
	}
}

void analyse_97_lifting_update(vector<double>& x, float a)
{
	size_t p = x.size();
	size_t limit = p * 0.5;
	for(size_t i = 0; i < limit; ++i)
	{
		x[2 * i] = a * x[2 * i + mirror_increment(i, -1, p)]+x[2 * i] + a * x[2 * i + mirror_increment(i, 1, p)];
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
float error(const vector<double>& x , const vector<double>& y)
{
	size_t p = x.size();

	float error = 0;
	for(size_t i = 0; i < p; ++i)
	{
		error += (x[i] - y[i]) * (x[i] - y[i]); 
	}

	return error;
}
/* Main */
int main (int argc, char* argv[])
{
/* Rampe */
	vector<double> rampe;
	for(int i = 0; i < 255; i++)
		rampe.push_back(i);
	save_signal(rampe,"rampe.txt");

	/* Haar */
		/* Analyse de Haar */
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

	/* biorthogonaux 9/7 */
		/* Analyse biorthogonaux 9/7 */
	cout <<"\n Analyse de biorthogonaux 9/7" << endl;
	analyse_97(leleccum);
	save_signal(leleccum,"./ouput_data/leleccum_analyse_97.txt");

		/* Reconstition biorthogonaux 9/7 */
	cout <<"\n Synthese de biorthogonaux 9/7" << endl;
	synthese_97(leleccum);
	save_signal(leleccum,"./ouput_data/leleccum_synthese_97.txt");

/* Leleccum */
	leleccum = leleccumOriginal;

	cout <<"\n Lifting de biorthogonaux 9/7" << endl;
	analyse_97_lifting(leleccum);
	save_signal(leleccum,"./ouput_data/leleccum_analyse_97_lifting.txt");

	synthese_97_lifting(leleccum);
	save_signal(leleccum,"./ouput_data/leleccum_synthese_97_lifting.txt");

	return 1;
}