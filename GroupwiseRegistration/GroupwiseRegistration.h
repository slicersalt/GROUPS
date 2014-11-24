#pragma once
#include <algorithm>
#include <iostream>
#include <time.h>
#include <assert.h>
#include <vector>
#include "Mesh.h"
#include "newuoa.h"
#include "AABB.h"

using namespace std;

class GroupwiseRegistration
{
public:
	GroupwiseRegistration(void);
	GroupwiseRegistration(char *sphere, char **tmpDepth, char **subjDepth, char **coeff, char **correspondence, int nSubj, int deg = 5, char *coeffLog = NULL, int nProperties = 1, int maxIter = 30000);
	GroupwiseRegistration(char *sphere, char **tmpDepth, char **subjDepth, int nSubj, int deg = 5, int nProperties = 1, bool propLoc = false, char *coeffLog = NULL, char **coeff = NULL, int maxIter = 30000);
	~GroupwiseRegistration(void);
	void saveLDeformation(char *filename);
	void saveLCoeff(char *filename, int id);
	float cost(float *coeff, int statusStep = 10);

private:
	void init_multi(char *sphere, char **tmpDepth, char **subjDepth, char **coeff, int nSubj, int deg = 5, int nProperties = 1, bool propLoc = false);
	void init(char *sphere, char **tmpDepth, char **subjDepth, char **coeff, char **correspondence, int nSubj, int deg = 5, int nProperties = 1);
	void reconsCoord(const float *v0, float *v1, float *Y, float **coeff, float degree, float *pole);
	void updateDeformation(int subject);
	void optimization(void);
	void eigenvalues(float *M, int dim, float *eig);
	int icosahedron(int degree);
	float landmarkEntropyMulti(void);
	float landmarkEntropy(void);
	float landmarkEntropyMedian(void);
	float depthVariance(void);
	float depthInterpolation(float *refMap, int index, float *coeff, Mesh *mesh);

private:
	struct point
	{
		float p[3];
		float id;
		float *Y;
		int subject;
	};
	struct sphereVertex
	{
		float v[3];
		float *Y;
	};
	struct spharm
	{
		int degree;
		float **coeff;
		float pole[3];
		vector<sphereVertex *> vertex;
		AABB *tree;
		Mesh *sphere;
		float *depth;
		int *depth_cache;
		float *meanDepth;
		float *maxDepth;
		float *minDepth;
	};
	struct entropy
	{
		bool occupied;
		vector<point *> samples;
		float p[3];
	};

	int m_nSubj;
	int m_csize;
	int m_nProperties;
	int m_maxIter;
	float *m_depth;
	float *m_meanDepth;
	float *m_maxDepth;
	float *m_minDepth;
	float *m_cov_weight;
	bool *m_updated;
	spharm *m_spharm;
	Mesh *m_sphere;
	float *m_coeff;
	AABB *m_tree;
	vector<point *> m_point;
	vector<float *> m_depthvar;
	vector<int> m_entropyList;
	entropy *m_entropy;
	char **m_coeff_fn;

	// work space
	float *m_cov_depth;
	float *m_pointList;
	float *m_cov_eig;

	// tic
	int nIter;

	// log
	FILE *m_clfp;
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
