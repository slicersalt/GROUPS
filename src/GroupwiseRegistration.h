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
	GroupwiseRegistration(char *sphere, char **tmpDepth, char **subjDepth, char **coeff, char **correspondence, char **asphere, int nSubj, const float *weight, int deg = 5, char *coeffLog = NULL, int nProperties = 1, int maxIter = 30000);
	GroupwiseRegistration(char *sphere, char **tmpDepth, char **subjDepth, int nSubj, char **landmark, char **asphere, const float *weight, int deg = 5, int nProperties = 1, float propLoc = 0, char *tmpSurf = NULL, char **surf = NULL, char *coeffLog = NULL, char **coeff = NULL, int maxIter = 30000, char **output = NULL);
	~GroupwiseRegistration(void);
	void saveLDeformation(char *filename);
	void saveLCoeff(char *filename, int id);
	float cost(float *coeff, int statusStep = 10);

private:
	void init_multi(char *sphere, char **tmpDepth, char **subjDepth, char **coeff, char **asphere, int nSubj, const float *weight, int deg = 5, int nProperties = 1, float propLoc = 0, char *tmpSurf = NULL, char **surf = NULL);
	void init(char *sphere, char **tmpDepth, char **subjDepth, char **coeff, char **correspondence, char **asphere, int nSubj, const float *weight, int deg = 5, int nProperties = 1, float propLoc = 0, char *tmpSurf = NULL, char **surf = NULL);
	void reconsCoord(const float *v0, float *v1, float *Y, float **coeff, float degree, float *pole);
	void updateDeformation(int subject);
	void optimization(void);
	void eigenvalues(float *M, int dim, float *eig);
	int testFlip(Mesh *mesh, const bool *flip);
	int icosahedron(int degree);
	float landmarkEntropyMulti(void);
	float landmarkEntropy(void);
	float landmarkEntropyMedian(void);
	float depthVariance(void);
	float depthInterpolation(float *refMap, int index, float *coeff, Mesh *mesh);
	float flipCost(void);

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
		Mesh *surf;
		float *depth;
		int *depth_cache;
		float *meanDepth;
		float *maxDepth;
		float *minDepth;
		float *sdevDepth;
		float *edge_var;
		bool *flip;
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
	float *m_sdevDepth;
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
	int m_nEdges;
	float m_minscore;
	bool m_prop_only;
	
	// work space
	float *m_cov_depth;
	float *m_pointList;
	float *m_cov_eig;
	float *m_edgeVar;

	// tic
	int nIter;

	// log
	FILE *m_clfp;
	
	// output list
	char **m_output;
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
