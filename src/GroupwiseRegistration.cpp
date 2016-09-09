/*************************************************
*	GroupwiseRegistration.cpp
*
*	Release: Sep 2016
*	Update: Sep 2016
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include <cstring>
#include <float.h>
#include "GroupwiseRegistration.h"
#include "SphericalHarmonics.h"
#include <lapacke.h>
#include "newuoa.h"

GroupwiseRegistration::GroupwiseRegistration(void)
{
	m_maxIter = 0;
	m_nSubj = 0;
	m_mincost = FLT_MAX;
	m_nProperties = 0;
	m_nSurfaceProperties = 0;
	m_output = NULL;
	m_degree = 0;
	m_degree_inc = 1;	// starting degree for the incremental optimization
}

GroupwiseRegistration::GroupwiseRegistration(const char **sphere, int nSubj, const char **property, int nProperties, const char **output, const float *weight, int deg, const char **landmark, float weightLoc, const char **coeff, const char **surf, int maxIter)
{
	m_maxIter = maxIter;
	m_nSubj = nSubj;
	m_mincost = FLT_MAX;
	m_nProperties = nProperties;
	m_nSurfaceProperties = (weightLoc > 0)? 3: 0;
	m_output = output;
	m_degree = deg;
	m_degree_inc = 3;	// starting degree for the incremental optimization
	init(sphere, property, weight, landmark, weightLoc, coeff, surf, 4);
}

GroupwiseRegistration::~GroupwiseRegistration(void)
{
	delete [] m_cov;
	delete [] m_feature_weight;
	delete [] m_eig;
	delete [] m_feature;
	delete [] m_updated;
	delete [] m_work;
	delete [] m_coeff;
	delete [] m_coeff_prev_step;
	for (int subj = 0; subj < m_nSubj; subj++)
	{
		delete m_spharm[subj].tree;
		delete m_spharm[subj].surf;
		delete m_spharm[subj].sphere;
		delete [] m_spharm[subj].coeff;
		delete [] m_spharm[subj].coeff_prev_step;
		delete [] m_spharm[subj].tree_cache;
		delete [] m_spharm[subj].meanProperty;
		delete [] m_spharm[subj].maxProperty;
		delete [] m_spharm[subj].minProperty;
		delete [] m_spharm[subj].sdevProperty;
		delete [] m_spharm[subj].property;
		delete [] m_spharm[subj].flip;
	}
	delete [] m_spharm;
	for (int i = 0; i < m_propertySamples.size(); i++)
		delete m_propertySamples[i];
}

void GroupwiseRegistration::run(void)
{
	cout << "Optimization\n";
	optimization();

	// write the solutions
	for (int subj = 0; subj < m_nSubj; subj++)
	{
		saveCoeff(m_output[subj], subj);
	}
	cout << "All done!\n";
}

void GroupwiseRegistration::init(const char **sphere, const char **property, const float *weight, const char **landmark, float weightLoc, const char **coeff, const char **surf, int samplingDegree)
{
	m_spharm = new spharm[m_nSubj];	// spharm info
	m_updated = new bool[m_nSubj];	// AABB tree cache
	m_eig = new float[m_nSubj];		// eigenvalues
	m_work = new float[m_nSubj * 3 - 1];	// workspace for eigenvalue computation
	m_csize = (m_degree + 1) * (m_degree + 1) * m_nSubj; // total # of coefficients
	m_coeff = new float[m_csize * 2];	// how many coefficients are required: the sum of all possible coefficients
	m_coeff_prev_step = new float[m_csize * 2];	// the previous coefficients

	// set all the coefficient to zeros
	memset(m_coeff, 0, sizeof(float) * m_csize * 2);
	memset(m_coeff_prev_step, 0, sizeof(float) * m_csize * 2);
	memset(m_updated, 0, sizeof(bool) * m_nSubj);
	
	cout << "Initialzation of subject information\n";

	for (int subj = 0; subj < m_nSubj; subj++)
	{
		cout << "Subject " << subj << " - " << sphere[subj] << endl;
		
		// spehre and surface information
		cout << "-Sphere information\n";
		m_spharm[subj].sphere = new Mesh();
		if (sphere != NULL)
		{
			m_spharm[subj].sphere->openFile(sphere[subj]);
			// make sure a unit sphere
			m_spharm[subj].sphere->centering();
			for (int i = 0; i < m_spharm[subj].sphere->nVertex(); i++)
			{
				Vertex *v = (Vertex *)m_spharm[subj].sphere->vertex(i);	// vertex information on the sphere
				const float *v0 = v->fv();
				Vector V(v0); V.unit();
				v->setVertex(V.fv());
			}
		}
		else cout << " Fatal error: No sphere mapping is provided!\n";
		
		// previous spherical harmonic deformation fields
		cout << "-Spherical harmonics information\n";
		initSphericalHarmonics(subj, coeff);
		updateDeformation(subj);	// deform the sphere for efficient AABB tree creation
		
		if (m_nSurfaceProperties > 0)
		{
			cout << "-Location information\n";
			m_spharm[subj].surf = new Mesh();
			m_spharm[subj].surf->openFile(surf[subj]);
		}
		else m_spharm[subj].surf = NULL;
		
		if (m_nProperties + m_nSurfaceProperties > 0)
		{
			// AABB tree construction for speedup computation
			cout << "-AABB tree construction\n";
			m_spharm[subj].tree = new AABB(m_spharm[subj].sphere);
		}
		else m_spharm[subj].tree = NULL;
		
		// triangle flipping
		cout << "-Triangle flipping\n";
		initTriangleFlipping(subj);
		
		// property information
		cout << "-Property information\n";
		initProperties(subj, property, 3);
		
		// landmarks
		if (landmark != NULL)
		{
			cout << "-Landmark information\n";
			initLandmarks(subj, landmark);
		}
		cout << "----------" << endl;
	}

	// icosahedron subdivision for evaluation on properties: this generates uniform sampling points over the sphere - m_propertySamples
	if (m_nProperties + m_nSurfaceProperties > 0) icosahedron(samplingDegree);
	// landmark information - the number of landamrks should be the same across subjects
	if (landmark != NULL)
	{
		int nLandmark = m_spharm[0].landmark.size();
		for (int subj = 0; subj < m_nSubj; subj++)
			if (nLandmark != m_spharm[subj].landmark.size())
				cout << " Fatal error: # of landamrks should be agreed!\n";
	}

	cout << "Computing weight terms\n";
	int nLandmark = m_spharm[0].landmark.size() * 3;	// # of landmarks: we assume all the subject has the same number, which already is checked above.
	int nSamples = m_propertySamples.size();	// # of sampling points for property map agreement
	
	// weights for covariance matrix computation
	int nTotalProperties = m_nProperties + m_nSurfaceProperties;	// if location information is provided, total number = # of property + 3 -> (x, y, z location)
	m_feature_weight = new float[nLandmark + nSamples * nTotalProperties];
	float landmarkWeight = (nLandmark > 0) ? (float)nSamples / (float)nLandmark: 0;	// based on the number ratio (balance between landmark and property)
	float totalWeight = weightLoc;
	for (int n = 0; n < m_nProperties; n++) totalWeight += weight[n];
	landmarkWeight *= totalWeight;
	if (landmarkWeight == 0) landmarkWeight = 1;
	cout << "Total properties: " << nTotalProperties << endl;
	cout << "Sampling points: " << nSamples << endl;

	// assign the weighting factors
	for (int i = 0; i < nLandmark; i++) m_feature_weight[i] = landmarkWeight;
	for (int n = 0; n < m_nProperties; n++)
		for (int i = 0; i < nSamples; i++)
			m_feature_weight[nLandmark + nSamples * n + i] = weight[n];
	// weight for location information
	for (int n = 0; n < m_nSurfaceProperties; n++)
		for (int i = 0; i < nSamples; i++)
			m_feature_weight[nLandmark + nSamples * (m_nProperties + n) + i] = weightLoc;
	
	if (nLandmark > 0) cout << "Landmark weight: " << landmarkWeight << endl;
	if (m_nProperties > 0)
	{
		cout << "Property weight: ";
		for (int i = 0; i < m_nProperties; i++) cout << weight[i] << " ";
		cout << endl;
	}
	if (weightLoc > 0) cout << "Location weight: " << weightLoc << endl;
	
	cout << "Initialization of work space\n";
	m_cov = new float[m_nSubj * m_nSubj];	// convariance matrix defined in the duel space with dimensions: nSubj x nSubj
	m_feature = new float[m_nSubj * (nLandmark + nSamples * nTotalProperties)];	// the entire feature vector map for optimization

	// AABB tree cache for each subject: this stores the closest face of the sampling point to the corresponding face on the input sphere model
	for (int subj = 0; subj < m_nSubj; subj++)
	{
		if (nTotalProperties > 0)
		{
			m_spharm[subj].tree_cache = new int[nSamples];
			for (int i = 0; i < nSamples; i++)
				m_spharm[subj].tree_cache[i] = -1;	// initially, set to -1 (invalid index)
		}
		else m_spharm[subj].tree_cache = NULL;
	}

	// inital coefficients for the previous step
	memcpy(m_coeff_prev_step, m_coeff, sizeof(float) * m_csize * 2);

	cout << "Feature vector creation\n";
	// set all the tree needs to be updated
	memset(m_updated, 0, sizeof(bool) * m_nSubj);

	// feature update
	if (nLandmark > 0) updateLandmark(); // update landmark
	if (nSamples > 0) updateProperties(); // update properties

	cout << "Initialization done!" << endl;
}

void GroupwiseRegistration::initSphericalHarmonics(int subj, const char **coeff)
{
	// spherical harmonics information
	int n = (m_degree + 1) * (m_degree + 1);	// total number of coefficients (this must be the same across all the subjects at the end of this program)
	// new memory allocation for coefficients
	m_spharm[subj].coeff = new float*[n * 2];
	m_spharm[subj].coeff_prev_step = new float*[n * 2];
	
	for (int i = 0; i < n; i++)
	{
		// store coefficients by asc order (low to high frequencies)
		m_spharm[subj].coeff[i] = &m_coeff[m_nSubj * 2 * i + subj * 2];	// latitudes
		m_spharm[subj].coeff[n + i] = &m_coeff[m_nSubj * 2 * i + subj * 2 + 1];	// longitudes
		m_spharm[subj].coeff_prev_step[i] = &m_coeff_prev_step[m_nSubj * 2 * i + subj * 2];	// latitudes
		m_spharm[subj].coeff_prev_step[n + i] = &m_coeff_prev_step[m_nSubj * 2 * i + subj * 2 + 1];	// longitudes
	}

	if (coeff != NULL)	// previous spherical harmonics information
	{
		FILE *fp = fopen(coeff[subj],"r");
		fscanf(fp, "%f %f %f", &m_spharm[subj].pole[0], &m_spharm[subj].pole[1], &m_spharm[subj].pole[2]);	// optimal pole information
		fscanf(fp, "%d", &m_spharm[subj].degree);	// previous deformation field degree

		if (m_spharm[subj].degree > m_degree) m_spharm[subj].degree = m_degree;	// if the previous degree is larger than desired one, just crop it.

		// load previous coefficient information
		for (int i = 0; i < (m_spharm[subj].degree + 1) * (m_spharm[subj].degree + 1); i++)
			fscanf(fp, "%f %f", m_spharm[subj].coeff[i], m_spharm[subj].coeff[n + i]);
		fclose(fp);
	}
	else	// no spherical harmonic information is provided
	{
		// just set the pole to [0, 0, 1]
		m_spharm[subj].pole[0] = 0;
		m_spharm[subj].pole[1] = 0;
		m_spharm[subj].pole[2] = 1;
	}
	m_spharm[subj].degree = m_degree;
	
	// build spherical harmonic basis functions for each vertex
	int nVertex = m_spharm[subj].sphere->nVertex();
	for (int i = 0; i < nVertex; i++)
	{
		// vertex information
		point *p = new point();	// new spherical information allocation
		Vertex *v = (Vertex *)m_spharm[subj].sphere->vertex(i);	// vertex information on the sphere
		const float *v0 = v->fv();
		p->Y = new float[(m_degree + 1) * (m_degree + 1)];
		p->p[0] = v0[0]; p->p[1] = v0[1]; p->p[2] = v0[2];
		p->id = i;
		p->subject = subj;
		SphericalHarmonics::basis(m_degree, p->p, p->Y);
		m_spharm[subj].vertex.push_back(p);
	}
}

void GroupwiseRegistration::initProperties(int subj, const char **property, int nHeaderLines)
{
	int nVertex = m_spharm[subj].sphere->nVertex();	// this is the same as the number of properties
	int nFace = m_spharm[subj].sphere->nFace();
	if (m_nProperties + m_nSurfaceProperties > 0)
	{
		m_spharm[subj].meanProperty = new float[m_nProperties + m_nSurfaceProperties];
		m_spharm[subj].maxProperty = new float[m_nProperties + m_nSurfaceProperties];
		m_spharm[subj].minProperty = new float[m_nProperties + m_nSurfaceProperties];
		m_spharm[subj].sdevProperty = new float[m_nProperties + m_nSurfaceProperties];
		m_spharm[subj].property = new float[(m_nProperties + m_nSurfaceProperties) * nVertex];
	}
	else
	{
		m_spharm[subj].meanProperty = NULL;
		m_spharm[subj].maxProperty = NULL;
		m_spharm[subj].minProperty = NULL;
		m_spharm[subj].sdevProperty = NULL;
		m_spharm[subj].property = NULL;
	}
	for (int i = 0; i < m_nProperties; i++)	// property information
	{
		int index = subj * m_nProperties + i;
		cout << "\t" << property[index] << endl;
		FILE *fp = fopen(property[index], "r");
		
		// remove header lines
		char line[1024];
		for (int j = 0; j < nHeaderLines; j++) fgets(line, sizeof(line), fp);
		
		// load property information
		for (int j = 0; j < nVertex; j++) fscanf(fp, "%f", &m_spharm[subj].property[nVertex * i + j]);
		fclose(fp);
	}
	for (int i = 0; i < m_nSurfaceProperties; i++)	// x, y, z dimensions
	{
		for (int j = 0; j < nVertex; j++)
		{
			Vertex *v = (Vertex *)m_spharm[subj].surf->vertex(j);
			const float *v0 = v->fv();
			m_spharm[subj].property[nVertex * (m_nProperties + i) + j] = v0[i];
		}
	}
	
	// find the best statistics across subjects
	for (int i = 0; i < m_nProperties + m_nSurfaceProperties; i++)
	{
		cout << "--Property " << i << endl;
		m_spharm[subj].meanProperty[i] = Statistics::mean(&m_spharm[subj].property[nVertex * i], nVertex);
		m_spharm[subj].maxProperty[i] = Statistics::max(&m_spharm[subj].property[nVertex * i], nVertex);
		m_spharm[subj].minProperty[i] = Statistics::min(&m_spharm[subj].property[nVertex * i], nVertex);
		m_spharm[subj].sdevProperty[i] = sqrt(Statistics::var(&m_spharm[subj].property[nVertex * i], nVertex));
		cout << "---Min/Max: " << m_spharm[subj].minProperty[i] << ", " << m_spharm[subj].maxProperty[i] << endl;
		cout << "---Mean/Stdev: " << m_spharm[subj].meanProperty[i] << ", " << m_spharm[subj].sdevProperty[i] << endl;
	}
}

void GroupwiseRegistration::initTriangleFlipping(int subj)
{
	int nFace = m_spharm[subj].sphere->nFace();
	m_spharm[subj].flip = new bool[nFace];
	
	// check triangle flips
	for (int i = 0; i < nFace; i++)
	{
		const float *v1 = m_spharm[subj].sphere->vertex(m_spharm[subj].sphere->face(i)->list(0))->fv();
		const float *v2 = m_spharm[subj].sphere->vertex(m_spharm[subj].sphere->face(i)->list(1))->fv();
		const float *v3 = m_spharm[subj].sphere->vertex(m_spharm[subj].sphere->face(i)->list(2))->fv();

		Vector V1(v1), V2(v2), V3(v3);
		Vector V = (V1 + v2 + V3) / 3;

		Vector U = (V2 - V1).cross(V3 - V1);

		if (V * U < 0) m_spharm[subj].flip[i] = true;
		else m_spharm[subj].flip[i] = false;
	}
}

void GroupwiseRegistration::initLandmarks(int subj, const char **landmark)
{
	FILE *fp = fopen(landmark[subj], "r");
	int i = 0;
	while (!feof(fp))
	{
		// indices for corresponding points
		int srcid;

		// landmark information
		int id;
		if (fscanf(fp, "%d", &id) == -1) break;
		const float *v = m_spharm[subj].sphere->vertex(id)->fv();
		
		float *Y = new float[(m_degree + 1) * (m_degree + 1)];
		point *p = new point();
		p->p[0] = v[0]; p->p[1] = v[1]; p->p[2] = v[2];
		SphericalHarmonics::basis(m_degree, p->p, Y);
		p->subject = subj;
		p->Y = Y;
		p->id = id;

		m_spharm[subj].landmark.push_back(p);
		
		i++;
	}
	fclose(fp);
}

bool GroupwiseRegistration::updateCoordinate(const float *v0, float *v1, const float *Y, const float **coeff, float degree, const float *pole)
{
	// spharm basis
	int n = (degree + 1) * (degree + 1);

	Vector p0(pole), axis;
	float dot;

	// fit to the equator
	float rot[9];
	Vector v(v0);
	axis = p0.cross(v);
	if (axis.norm() == 0)	// point == pole
	{
		memcpy(v1, v0, sizeof(float) * 3);
		return false;
	}

	float delta[2] = {0, 0};
	for (int i = 0; i < n; i++)
	{
		delta[0] += Y[i] * *coeff[i];
		delta[1] += Y[i] * *coeff[(m_degree + 1) * (m_degree + 1) + i];
	}

	if (delta[0] == 0 && delta[1] == 0)
	{
		memcpy(v1, v0, sizeof(float) * 3);
		return false;
	}

	dot = p0 * v;
	dot = (dot > 1) ? 1: dot;
	dot = (dot < -1) ? -1: dot;
	float deg = PI / 2 - acos(dot);
	Coordinate::rotation(axis.fv(), deg, rot);

	// rotation to the eqautor
	float rv[3];
	Coordinate::rotPoint(v0, rot, rv);

	// polar coodinate
	float phi, theta;
	Coordinate::cart2sph(rv, &phi, &theta);
	
	// displacement
	phi += delta[0];	// longitude (azimuth) change
	Coordinate::sph2cart(phi, theta, rv);
	
	// rotation to the new longitude change
	// a simple inverse rotation is not enough due to not exact location back
	Vector u(rv);
	axis = p0.cross(u);
	if (axis.norm() == 0) axis = p0;
	Coordinate::rotation(axis.fv(), -deg, rot);
	
	theta += delta[1];
	Coordinate::sph2cart(phi, theta, rv);	// locally normalized polar system

	// inverse rotation
	Coordinate::rotPoint(rv, rot, v1);
	
	return true;
}

void GroupwiseRegistration::updateDeformation(int subject)
{
	// note: the deformation happens only if the coefficients change; otherwise, nothing to do
	bool updated = m_updated[subject];
	
	// check if the coefficients change
	int n = (m_degree_inc + 1) * (m_degree_inc + 1);
	for (int i = 0; i < n && updated; i++)
		if (*m_spharm[subject].coeff[i] != *m_spharm[subject].coeff_prev_step[i] ||
			*m_spharm[subject].coeff[(m_degree + 1) * (m_degree + 1) + i] != *m_spharm[subject].coeff_prev_step[(m_degree + 1) * (m_degree + 1) + i])
			updated = false;
	
	// deform a sphere based on the current coefficients if necessary
	for (int i = 0; i < m_spharm[subject].vertex.size() && !updated; i++)
	{
		Vertex *v = (Vertex *)m_spharm[subject].sphere->vertex(i);
		float v1[3];
		const float *v0 = v->fv();
		updateCoordinate(m_spharm[subject].vertex[i]->p, v1, m_spharm[subject].vertex[i]->Y, (const float **)m_spharm[subject].coeff, m_degree_inc, m_spharm[subject].pole); // update using the current incremental degree
		{
			Vector V(v1); V.unit();
			v->setVertex(V.fv());
		}
	}
	m_updated[subject] = updated;
}

void GroupwiseRegistration::updateLandmark(void)
{
	int nLandmark = m_spharm[0].landmark.size();
	int nSamples = m_propertySamples.size();	// # of sampling points for property map agreemen

	for (int i = 0; i < nLandmark; i++)
	{
		float m[3] = {0, 0, 0};	// mean
		for (int subj = 0; subj < m_nSubj; subj++)
		{
			int id = m_spharm[subj].landmark[i]->id;
			updateCoordinate(m_spharm[subj].landmark[i]->p, &m_feature[subj * (nLandmark * 3 + nSamples * (m_nProperties + m_nSurfaceProperties)) + i * 3], m_spharm[subj].landmark[i]->Y, (const float **)m_spharm[subj].coeff, m_degree_inc, m_spharm[subj].pole);

			// mean locations
			for (int k = 0; k < 3; k++) m[k] += m_feature[subj * (nLandmark * 3 + nSamples * (m_nProperties + m_nSurfaceProperties)) + i * 3 + k];
		}
		// forcing the mean to be on the sphere
		float norm = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
		for (int k = 0; k < 3; k++) m[k] /= norm;
		
		// projection
		for (int subj = 0; subj < m_nSubj; subj++)
		{
			float newp[3];
			Coordinate::proj2plane(m[0], m[1], m[2], -1, &m_feature[subj * (nLandmark * 3 + nSamples * (m_nProperties + m_nSurfaceProperties)) + i * 3], newp);
			memcpy(&m_feature[subj * (nLandmark * 3 + nSamples * (m_nProperties + m_nSurfaceProperties)) + i * 3], newp, sizeof(float) * 3);
		}
	}
}

void GroupwiseRegistration::updateLandmarkMedian(void)
{
	int nLandmark = m_spharm[0].landmark.size();
	int nSamples = m_propertySamples.size();	// # of sampling points for property map agreemen
	
	// m-estimator
	float m1 = 0.03 * 0.01;
	float m2 = 16.8 * 0.01;
	float u = (m2 + m1) / 2;
	float sigma = ((m2 + m1) / 2 - m1) / 3;

	// deformed points
	float *x = new float[m_nSubj];
	float *y = new float[m_nSubj];
	float *z = new float[m_nSubj];
	
	for (int i = 0; i < nLandmark; i++)
	{
		float m[3] = {0, 0, 0};	// median
		
		for (int subj = 0; subj < m_nSubj; subj++)
		{
			int id = m_spharm[subj].landmark[i]->id;
			updateCoordinate(m_spharm[subj].landmark[i]->p, &m_feature[subj * (nLandmark * 3 + nSamples * (m_nProperties + m_nSurfaceProperties)) + i * 3], m_spharm[subj].landmark[i]->Y, (const float **)m_spharm[subj].coeff, m_degree_inc, m_spharm[subj].pole);

			// median locations
			x[subj] = m_feature[subj * (nLandmark * 3 + nSamples * (m_nProperties + m_nSurfaceProperties)) + i * 3 + 0];
			y[subj] = m_feature[subj * (nLandmark * 3 + nSamples * (m_nProperties + m_nSurfaceProperties)) + i * 3 + 1];
			z[subj] = m_feature[subj * (nLandmark * 3 + nSamples * (m_nProperties + m_nSurfaceProperties)) + i * 3 + 2];
		}
		m[0] = Statistics::median(x, m_nSubj);
		m[1] = Statistics::median(y, m_nSubj);
		m[2] = Statistics::median(z, m_nSubj);

		// forcing the mean to be on the sphere
		float norm = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
		for (int k = 0; k < 3; k++) m[k] /= norm;

		// projection
		for (int subj = 0; subj < m_nSubj; subj++)
		{
			float newp[3];
			Coordinate::proj2plane(m[0], m[1], m[2], -1, &m_feature[subj * (nLandmark * 3 + nSamples * (m_nProperties + m_nSurfaceProperties)) + i * 3], newp);

			// m-estimator
			float len = sqrt((newp[0] - m[0]) * (newp[0] - m[0]) + (newp[1] - m[1]) * (newp[1] - m[1]) + (newp[2] - m[2]) * (newp[2] - m[2]));
			float slen = Statistics::normal_cdf(len, u, sigma);
			slen = (len < m2) ? len: m2;
			float ratio = 1.0f;
			if (len > 0) ratio = slen / len;

			newp[0] = m[0] + (newp[0] - m[0]) * ratio;
			newp[1] = m[1] + (newp[1] - m[1]) * ratio;
			newp[2] = m[2] + (newp[2] - m[2]) * ratio;

			memcpy(&m_feature[subj * (nLandmark * 3 + nSamples * (m_nProperties + m_nSurfaceProperties)) + i * 3], newp, sizeof(float) * 3);
		}
	}
	
	// free
	delete [] x;
	delete [] y;
	delete [] z;
}

void GroupwiseRegistration::updateProperties(void)
{
	int nLandmark = m_spharm[0].landmark.size();
	int nSamples = m_propertySamples.size();	// # of sampling points for property map agreement
	
	float err = 0;
	for (int subj = 0; subj < m_nSubj; subj++)
	{
		if (!m_updated[subj])
		{
			m_updated[subj] = true;
			m_spharm[subj].tree->update();
		}
		else continue;	// don't compute again since tree is the same as the previous. The feature vector won't be changed
		for (int i = 0; i < nSamples; i++)
		{
			int fid = -1;
			float coeff[3];
			if (m_spharm[subj].tree_cache[i] != -1)	// if previous cache is available
			{
				Face *f = (Face *)m_spharm[subj].sphere->face(m_spharm[subj].tree_cache[i]);
				Vertex *a = (Vertex *)f->vertex(0);
				Vertex *b = (Vertex *)f->vertex(1);
				Vertex *c = (Vertex *)f->vertex(2);

				// bary centric
				Coordinate::cart2bary((float *)a->fv(), (float *)b->fv(), (float *)c->fv(), m_propertySamples[i], coeff);

				fid = (coeff[0] >= err && coeff[1] >= err && coeff[2] >= err) ? m_spharm[subj].tree_cache[i]: -1;
			}
			if (fid == -1)	// if no closest face is found
			{
				fid = m_spharm[subj].tree->closestFace(m_propertySamples[i], coeff);
				if (fid == -1)	// something goes wrong
					cout << "Fatal error: no closest point found!\n";
			}

			int nVertex = m_spharm[subj].sphere->nVertex();
			for (int k = 0; k < m_nProperties + m_nSurfaceProperties; k++)
			{
				m_feature[subj * (nLandmark * 3 + nSamples * (m_nProperties + m_nSurfaceProperties)) + nLandmark * 3 + nSamples * k + i] = propertyInterpolation(&m_spharm[subj].property[nVertex * k], fid, coeff, m_spharm[subj].sphere) / m_spharm[subj].sdevProperty[k];
			}
			m_spharm[subj].tree_cache[i] = fid;
		}
	}
}

float GroupwiseRegistration::entropy(void)
{
	int nLandmark = m_spharm[0].landmark.size() * 3;	// # of landmarks: we assume all the subject has the same number
	int nSamples = m_propertySamples.size();	// # of sampling points for property map agreement
	
	float E = 0;	// entropy
	
	memset(m_cov, 0, sizeof(float) * m_nSubj * m_nSubj);	// reset covariance matrix

	// update landmark
	if (nLandmark > 0) updateLandmark();
	
	// update properties
	if (nSamples > 0) updateProperties();
	
	// dual covariance matrix (m_nSubj x m_nSubj) of feature vector (nLandmark + nSamples * (m_nProperties + m_nSurfaceProperties) x m_nSubj)
	Statistics::wcov_trans(m_feature, m_nSubj, nLandmark + nSamples * (m_nProperties + m_nSurfaceProperties), m_cov, m_feature_weight);
	
	// entropy
	eigenvalues(m_cov, m_nSubj, m_eig);

	float alpha = 1e-5;	// avoid a degenerative case
	for (int i = 1; i < m_nSubj; i++)	// just ignore the first eigenvalue (trivial = 0)
		E += log(m_eig[i] + alpha);

	return E;
}

void GroupwiseRegistration::eigenvalues(float *M, int dim, float *eig)
{
	int n = dim;
	int lwork = dim * 3 - 1;	// dimension of the work array
	int lda = n;			// lda: leading dimension
	int info;				// information (0 for successful exit)
	
	char jobz[] = "N";	// eigenvalue only
	char uplo[] = "L"; // Lower triangle
	ssyev_(jobz, uplo, &n, M, &lda, eig, m_work, &lwork, &info);
}

float GroupwiseRegistration::propertyInterpolation(float *refMap, int index, float *coeff, Mesh *mesh)
{
	float property = 0;

	if (index != -1)
	{
		Face *f = (Face *)mesh->face(index);
		Vertex *a = (Vertex *)f->vertex(0);
		Vertex *b = (Vertex *)f->vertex(1);
		Vertex *c = (Vertex *)f->vertex(2);
		property = refMap[a->id()] * coeff[0] + refMap[b->id()] * coeff[1] + refMap[c->id()] * coeff[2];
	}

	return property;
}

float GroupwiseRegistration::cost(float *coeff, int statusStep)
{
	// update defomation fields
	for (int i = 0; i < m_nSubj; i++) updateDeformation(i);
	
	// how many flips are detected
	int nFolds = 0;
	for (int i = 0; i < m_nSubj; i++)
		nFolds += testTriangleFlip(m_spharm[i].sphere, m_spharm[i].flip);

	float fcost = (nFolds == 0) ? 0: (nFolds + 1) * fabs(m_mincost);
	float ecost = (nFolds == 0) ? entropy(): m_mincost;

	float cost = ecost + fcost;
	if (m_mincost > cost)
	{
		m_mincost = cost;
		// write the current optimal solutions
		for (int subj = 0; subj < m_nSubj; subj++)
		{
			saveCoeff(m_output[subj], subj);
		}
	}
	
	if (nIter % statusStep == 0)
	{
		cout << "[" << nIter << "] " << cost << " (" << ecost << " + " << fcost << ")" << " " << m_mincost << endl;
	}
	nIter++;

	// copy previous coefficients
	memcpy(m_coeff_prev_step, m_coeff, sizeof(float) * m_csize * 2);
	
	return cost;
}

int GroupwiseRegistration::testTriangleFlip(Mesh *mesh, const bool *flip)
{
	int nFolds = 0;
	for (int i = 0; i < mesh->nFace(); i++)
	{
		const float *v1 = mesh->vertex(mesh->face(i)->list(0))->fv();
		const float *v2 = mesh->vertex(mesh->face(i)->list(1))->fv();
		const float *v3 = mesh->vertex(mesh->face(i)->list(2))->fv();

		Vector V1(v1), V2(v2), V3(v3);
		Vector V = (V1 + v2 + V3) / 3;

		Vector U = (V2 - V1).cross(V3 - V1);

		if ((V * U < 0 && !flip[i]) || (V * U > 0 && flip[i])) nFolds++;
		
		if (nFolds > 0) break;	// do not allow any flips!
	}
	return nFolds;
}

void GroupwiseRegistration::optimization(void)
{
	cost_function costFunc(this);
	int prev = 0;
	int step = 1;
	
	int n1 = (m_degree_inc + 1) * (m_degree_inc + 1) * m_nSubj * 2;
	int n2 = m_csize * 2 - n1;

	while (m_degree_inc < m_degree)
	{
		nIter = 0;
		int n = (m_degree_inc + 1) * (m_degree_inc + 1) * m_nSubj * 2 - prev;
		min_newuoa(n, &m_coeff[prev], costFunc, 1.0f, 1e-5f, m_maxIter);
		prev = (m_degree_inc + 1) * (m_degree_inc + 1) * m_nSubj * 2;
		m_degree_inc = min(m_degree_inc + step, m_degree);
	}
	
	// the entire optimization together
	nIter = 0;
	min_newuoa(m_csize * 2, m_coeff, costFunc, 1.0f, 1e-6f, m_maxIter);
}

int GroupwiseRegistration::icosahedron(int degree)
{
	// http://www.1activeserverpagesstreet.com/vb/scripts/ShowCode.asp?txtCodeId=9814&lngWId=3
	vector<Vector *> triangles;
	vector<Vector> vertices;

	float t = (1 + sqrt(5.0)) / 2.0;
	float s = sqrt(1 + t * t);

	// create the 12 vertices
	Vector v0 = Vector(t, 1.0, 0.0) / s;
	Vector v1 = Vector(-t, 1.0, 0.0) / s;
	Vector v2 = Vector(t, -1.0, 0.0) / s;
	Vector v3 = Vector(-t, -1.0, 0.0) / s;
	Vector v4 = Vector(1.0, 0.0, t) / s;
	Vector v5 = Vector(1.0, 0.0, -t) / s;
	Vector v6 = Vector(-1.0, 0.0, t) / s;
	Vector v7 = Vector(-1.0, 0.0, -t) / s;
	Vector v8 = Vector(0.0, t, 1.0) / s;
	Vector v9 = Vector(0.0, -t, 1.0) / s;
	Vector v10 = Vector(0.0, t, -1.0) / s;
	Vector v11 = Vector(0.0, -t, -1.0) / s;
    
	// create the 20 triangles
	Vector *f; 
	f = new Vector[3]; f[0] = v0; f[1] = v8; f[2] = v4; triangles.push_back(f);
	f = new Vector[3]; f[0] = v1; f[1] = v10; f[2] = v7; triangles.push_back(f);
	f = new Vector[3]; f[0] = v2; f[1] = v9; f[2] = v11; triangles.push_back(f);
	f = new Vector[3]; f[0] = v7; f[1] = v3; f[2] = v1; triangles.push_back(f);
	f = new Vector[3]; f[0] = v0; f[1] = v5; f[2] = v10; triangles.push_back(f);
	f = new Vector[3]; f[0] = v3; f[1] = v9; f[2] = v6; triangles.push_back(f);
	f = new Vector[3]; f[0] = v3; f[1] = v11; f[2] = v9; triangles.push_back(f);
	f = new Vector[3]; f[0] = v8; f[1] = v6; f[2] = v4; triangles.push_back(f);
	f = new Vector[3]; f[0] = v2; f[1] = v4; f[2] = v9; triangles.push_back(f);
	f = new Vector[3]; f[0] = v3; f[1] = v7; f[2] = v11; triangles.push_back(f);
	f = new Vector[3]; f[0] = v4; f[1] = v2; f[2] = v0; triangles.push_back(f);
	f = new Vector[3]; f[0] = v9; f[1] = v4; f[2] = v6; triangles.push_back(f);
	f = new Vector[3]; f[0] = v2; f[1] = v11; f[2] = v5; triangles.push_back(f);
	f = new Vector[3]; f[0] = v0; f[1] = v10; f[2] = v8; triangles.push_back(f);
	f = new Vector[3]; f[0] = v5; f[1] = v0; f[2] = v2; triangles.push_back(f);
	f = new Vector[3]; f[0] = v10; f[1] = v5; f[2] = v7; triangles.push_back(f);
	f = new Vector[3]; f[0] = v1; f[1] = v6; f[2] = v8; triangles.push_back(f);
	f = new Vector[3]; f[0] = v1; f[1] = v8; f[2] = v10; triangles.push_back(f);
	f = new Vector[3]; f[0] = v6; f[1] = v1; f[2] = v3; triangles.push_back(f);
	f = new Vector[3]; f[0] = v11; f[1] = v7; f[2] = v5; triangles.push_back(f);

	// subdivision
	for (int d = 0; d < degree; d++)
	{
		int nFaces = triangles.size();
		for (int i = 0 ; i < nFaces; i++)
		{
			Vector *f = triangles[i];
			Vector a = f[0], b = f[1], c = f[2];
			Vector v1 = a + b;
			Vector v2 = c + a;
			Vector v3 = b + c;
			// normalization
			v1.unit(); v2.unit(); v3.unit();
			f[0] = v1; f[1] = v3; f[2] = v2; // overwrite the original
			/*Vector f1[3] = {a, v1, v2}; triangles.push_back(f1);
			Vector f2[3] = {c, v2, v3}; triangles.push_back(f2);
			Vector f3[3] = {b, v3, v1}; triangles.push_back(f3);*/
			Vector *f1 = new Vector[3]; f1[0] = a; f1[1] = v1; f1[2] = v2; triangles.push_back(f1);
			Vector *f2 = new Vector[3]; f2[0] = c; f2[1] = v2; f2[2] = v3; triangles.push_back(f2);
			Vector *f3 = new Vector[3]; f3[0] = b; f3[1] = v3; f3[2] = v1; triangles.push_back(f3);
		}
	}
	for (int i = 0; i < triangles.size(); i++)
	{
		Vector *f = triangles[i];
		for (int j = 0; j < 3; j++)
			vertices.push_back(f[j]);
	}
	sort(vertices.begin(), vertices.end());
	vertices.erase(unique(vertices.begin(), vertices.end()), vertices.end());

	for (int i = 0; i < vertices.size(); i++)
	{
		float *p = new float[3];
		p[0] = vertices[i][0];
		p[1] = vertices[i][1];
		p[2] = vertices[i][2];
		m_propertySamples.push_back(p);
	}

	return m_propertySamples.size();
}

void GroupwiseRegistration::saveCoeff(const char *filename, int id)
{
	FILE *fp = fopen(filename, "w");
	fprintf(fp, "%f %f %f\n", m_spharm[id].pole[0], m_spharm[id].pole[1], m_spharm[id].pole[2]);
	fprintf(fp, "%d\n", m_spharm[id].degree);
	int n = (m_spharm[id].degree + 1) * (m_spharm[id].degree + 1);
	for (int i = 0; i < n; i++)
	{
		fprintf(fp, "%f %f\n", *m_spharm[id].coeff[i], *m_spharm[id].coeff[n + i]);
	}
	fclose(fp);
}
