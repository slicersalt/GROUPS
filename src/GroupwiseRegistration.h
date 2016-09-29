/*************************************************
*	GroupwiseRegistration.h
*
*	Release: Sep 2016
*	Update: Sep 2016
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#pragma once
#include <algorithm>
#include <iostream>
#include <vector>
#include "Mesh.h"
#include "AABB.h"

using namespace std;

class GroupwiseRegistration
{
public:
	GroupwiseRegistration(void);
	GroupwiseRegistration(const char **sphere, int nSubj, const char **property, int nProperties, const char **output, const float *weight, int deg = 5, const char **landmark = NULL, float weightLoc = 0, const char **coeff = NULL, const char **surf = NULL, int maxIter = 50000);
	~GroupwiseRegistration(void);
	void run(void);
	void saveCoeff(const char *filename, int id);
	float cost(float *coeff, int statusStep = 10);

private:
	// class members for initilaization
	void init(const char **sphere, const char **property, const float *weight, const char **landmark, float weightLoc, const char **coeff, const char **surf, int samplingDegree = 3);
	void initSphericalHarmonics(int subj, const char **coeff);
	void initTriangleFlipping(int subj);
	void initProperties(int subj, const char **property, int nHeaderLines);
	void initLandmarks(int subj, const char **landmark);
	int icosahedron(int degree);

	// entropy computation
	void optimization(void);
	void updateLandmark(void);
	void updateLandmarkMedian(void);
	void updateProperties(void);
	void eigenvalues(float *M, int dim, float *eig);
	float entropy(void);
	float propertyInterpolation(float *refMap, int index, float *coeff, Mesh *mesh);
	int testTriangleFlip(Mesh *mesh, const bool *flip);

	// deformation field reconstruction
	void updateDeformation(int subject);
	bool updateCoordinate(const float *v0, float *v1, const float *Y, const float **coeff, float degree, const float *pole);
	
private:
	struct point
	{
		float p[3];
		float id;
		float *Y;
		int subject;
	};
	struct spharm
	{
		int degree;
		float **coeff;
		float **coeff_prev_step;
		float pole[3];
		vector<point *> vertex;
		AABB *tree;
		Mesh *sphere;
		Mesh *surf;
		float *property;
		int *tree_cache;
		float *meanProperty;
		float *maxProperty;
		float *minProperty;
		float *sdevProperty;
		vector<point *> landmark;
		bool *flip;
	};

	int m_nSubj;
	int m_csize;
	int m_nProperties;
	int m_nSurfaceProperties;
	int m_maxIter;
	int m_degree;
	int m_degree_inc;	// incremental degree
	
	float *m_coeff;
	float *m_coeff_prev_step;	// previous coefficients
	bool *m_updated;
	spharm *m_spharm;
	vector<float *> m_propertySamples;
	
	float m_mincost;
	
	// work space for the entire procedure
	float *m_cov;
	float *m_feature;
	float *m_feature_weight;
	float *m_eig;
	float *m_work;	// for lapack eigenvalue computation
	
	// tic
	int nIter;

	// output list
	const char **m_output;
};

class cost_function
{
public:
    cost_function (GroupwiseRegistration *instance)
    {
        m_instance = instance;
    }

    double operator () (float *arg)
    {
		float cost = m_instance->cost(arg);
        return (double)cost;
    }

private:
	GroupwiseRegistration *m_instance;
};
