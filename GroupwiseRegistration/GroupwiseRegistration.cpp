#include "GroupwiseRegistration.h"
#include "SphericalHarmonics.h"
#include "lapack/lapacke.h"
#include "newuoa.h"
#include "string.h"

GroupwiseRegistration::GroupwiseRegistration(void)
{
}

GroupwiseRegistration::GroupwiseRegistration(char *sphere, char **tmpDepth, char **subjDepth, int nSubj, int deg, int nProperties, bool propLoc, char *tmpSurf, char **surf, char *coeffLog, char **coeff, int maxIter)
{
	m_maxIter = maxIter;
	m_nSubj = nSubj;
	m_nProperties = nProperties;
	init_multi(sphere, tmpDepth, subjDepth, coeff, nSubj, deg, nProperties, propLoc, tmpSurf, surf);
	if (coeffLog != NULL) m_clfp = fopen(coeffLog, "w");
	else m_clfp = NULL;
	optimization();
	if (coeffLog != NULL) fclose(m_clfp);
}

GroupwiseRegistration::GroupwiseRegistration(char *sphere, char **tmpDepth, char **subjDepth, char **coeff, char **correspondence, int nSubj, int deg, char *coeffLog, int nProperties, int maxIter)
{
	m_maxIter = maxIter;
	m_nSubj = nSubj;
	m_nProperties = nProperties;
	init(sphere, tmpDepth, subjDepth, coeff, correspondence, nSubj, deg, nProperties);
	if (coeffLog != NULL) m_clfp = fopen(coeffLog, "w");
	else m_clfp = NULL;
	optimization();
	if (coeffLog != NULL) fclose(m_clfp);
}

GroupwiseRegistration::~GroupwiseRegistration(void)
{
	delete [] m_coeff;
	delete [] m_cov_depth;
	delete [] m_cov_eig;
	delete [] m_pointList;
	delete [] m_cov_weight;
	delete [] m_updated;
}

void GroupwiseRegistration::init_multi(char *sphere, char **tmpDepth, char **subjDepth, char **coeff, int nSubj, int deg, int nProperties, bool propLoc, char *tmpSurf, char **surf)
{
	// unit sphere information
	cout << "Loading unit sphere information..\n";
	m_sphere = new Mesh();
	m_sphere->openFile(sphere);

	m_coeff_fn = coeff;
	m_spharm = new spharm[nSubj];
	m_entropy = new entropy[m_sphere->nFace()];

	int nDepth = m_sphere->nVertex();
	Mesh *tSurf;
	if (propLoc)
	{
		nProperties += 3;
		tSurf = new Mesh();
		tSurf->openFile(tmpSurf);
	}
	m_depth = new float[nDepth * nProperties];
	m_meanDepth = new float[nProperties];
	m_maxDepth = new float[nProperties];
	m_minDepth = new float[nProperties];
	
	int count = 0;
	for (int n = 0; n < nProperties; n++)
	{
		if (n >= nProperties - 3 && propLoc)
		{
			for (int i = 0; i < nDepth; i++)
			{
				Vertex *v = (Vertex *)tSurf->vertex(i);
				const float *v0 = v->fv();
				m_depth[nDepth * n + i] = v0[count];
			}
			m_meanDepth[n] = Statistics::mean(m_depth, nDepth);
			m_maxDepth[n] = Statistics::max(m_depth, nDepth);
			m_minDepth[n] = Statistics::min(m_depth, nDepth);
			count++;
		}
		else
		{
			FILE *fp = fopen(tmpDepth[n], "r");
			char line[1024];
			fgets(line, sizeof(line), fp);
			fgets(line, sizeof(line), fp);
			fgets(line, sizeof(line), fp);
			for (int i = 0; i < nDepth && !feof(fp); i++)
			{
				fscanf(fp, "%f", &m_depth[nDepth * n + i]);
				//if (m_depth[i] > 0) m_depth[i] = 0;
			}
			fclose(fp);
			m_meanDepth[n] = Statistics::mean(m_depth, nDepth);
			m_maxDepth[n] = Statistics::max(m_depth, nDepth) - m_meanDepth[n];
			m_minDepth[n] = Statistics::min(m_depth, nDepth) - m_meanDepth[n];
		}
	}

	//cout << "Depth Range: [" << m_minDepth << ", " << m_maxDepth << "]\n";

	for (int i = 0; i < m_sphere->nFace(); i++)
		m_entropy[i].occupied = false;

	m_tree = new AABB(m_sphere);

	cout << "Loading initial deformation..\n";
	m_csize = 0;

	if (coeff != NULL)
	{
		for (int subj = 0; subj < nSubj; subj++)
		{
			int degree;
			int pole[3];
			FILE *fp = fopen(coeff[subj],"r");
			fscanf(fp, "%f %f %f", &pole[0], &pole[1], &pole[2]);
			fscanf(fp, "%d", &degree);
			//if (degree < deg) degree = deg;
			degree = deg;
			m_csize += (degree + 1) * (degree + 1);
			fclose(fp);
		}
	}
	else
	{
		for (int subj = 0; subj < nSubj; subj++)
		{
			int degree = deg;
			int pole[3] = {0, 0, 1};
			m_csize += (degree + 1) * (degree + 1);
		}
	}
	m_coeff = new float[m_csize * 2];
	memset(m_coeff, 0, sizeof(float) * m_csize * 2);

	for (int subj = 0; subj < nSubj; subj++)
	{
		cout << "suject " << subj << endl;
		
		// spherical harmonics information
		m_spharm[subj] = spharm();
		int degree = deg;
		if (coeff != NULL)
		{
			FILE *fp = fopen(coeff[subj],"r");
			fscanf(fp, "%f %f %f", &m_spharm[subj].pole[0], &m_spharm[subj].pole[1], &m_spharm[subj].pole[2]);
			fscanf(fp, "%d", &m_spharm[subj].degree);
			/*int degree = m_spharm[subj].degree;
			if (degree < deg) degree = deg;*/
			if (m_spharm[subj].degree > deg) m_spharm[subj].degree = deg;

			int n = (degree + 1) * (degree + 1);
			m_spharm[subj].coeff = new float*[n * 2];
			for (int i = 0; i < n; i++)
			{
				m_spharm[subj].coeff[i] = &m_coeff[nSubj * 2 * i + subj * 2];
				m_spharm[subj].coeff[n + i] = &m_coeff[nSubj * 2 * i + subj * 2 + 1];
			}

			for (int i = 0; i < (m_spharm[subj].degree + 1) * (m_spharm[subj].degree + 1); i++)
				fscanf(fp, "%f %f", m_spharm[subj].coeff[i], m_spharm[subj].coeff[n + i]);
			fclose(fp);
		}
		else
		{
			m_spharm[subj].pole[0] = 0;
			m_spharm[subj].pole[1] = 0;
			m_spharm[subj].pole[2] = 1;
			m_spharm[subj].degree = deg;
			if (m_spharm[subj].degree > deg) m_spharm[subj].degree = deg;

			int n = (degree + 1) * (degree + 1);
			m_spharm[subj].coeff = new float*[n * 2];
			for (int i = 0; i < n; i++)
			{
				m_spharm[subj].coeff[i] = &m_coeff[nSubj * 2 * i + subj * 2];
				m_spharm[subj].coeff[n + i] = &m_coeff[nSubj * 2 * i + subj * 2 + 1];
			}

			for (int i = 0; i < (m_spharm[subj].degree + 1) * (m_spharm[subj].degree + 1); i++)
			{
				*m_spharm[subj].coeff[i] = 0;
				*m_spharm[subj].coeff[n + i] = 0;
			}
		}

		m_spharm[subj].degree = degree;

		// sphere information
		m_spharm[subj].sphere = new Mesh();
		m_spharm[subj].sphere->openFile(sphere);
		
		if (propLoc)
		{
			m_spharm[subj].surf = new Mesh();
			m_spharm[subj].surf->openFile(surf[subj]);
		}

		// vertices on the sphere
		for (int i = 0; i < m_spharm[subj].sphere->nVertex(); i++)
		{
			sphereVertex *p = new sphereVertex();
			Vertex *v = (Vertex *)m_spharm[subj].sphere->vertex(i);
			const float *v0 = v->fv();
			p->Y = new float[(degree + 1) * (degree + 1)];
			p->v[0] = v0[0]; p->v[1] = v0[1]; p->v[2] = v0[2];
			SphericalHarmonics::basis(degree, p->v, p->Y);
			m_spharm[subj].vertex.push_back(p);
		}

		// deformed irregular sphere
		m_updated = new bool[subj];
		cout << "Updating sphere deformation.. ";
		updateDeformation(subj);
		cout << "done\n";

		// AABB tree construction
		cout << "Constructing a search tree.. ";
		m_spharm[subj].tree = new AABB(m_spharm[subj].sphere);
		cout << "done\n";

		// depth information
		cout << "Loading depth information.. \n";
		int nDepth = m_spharm[subj].sphere->nVertex();
		m_spharm[subj].meanDepth = new float[nProperties];
		m_spharm[subj].maxDepth = new float[nProperties];
		m_spharm[subj].minDepth = new float[nProperties];
		m_spharm[subj].depth = new float[nDepth * nProperties];
		
		count = 0;
		for (int n = 0; n < nProperties; n++)
		{
			if (n >= nProperties - 3 && propLoc)
			{
				for (int i = 0; i < nDepth; i++)
				{
					Vertex *v = (Vertex *)m_spharm[subj].surf->vertex(i);
					const float *v0 = v->fv();
					m_spharm[subj].depth[nDepth * n + i] = v0[count];
				}
				count++;
				m_spharm[subj].meanDepth[n] = Statistics::mean(m_spharm[subj].depth, nDepth);
				m_spharm[subj].maxDepth[n] = Statistics::max(m_spharm[subj].depth, nDepth);
				m_spharm[subj].minDepth[n] = Statistics::min(m_spharm[subj].depth, nDepth);
				if (m_maxDepth[n] < m_spharm[subj].maxDepth[n]) m_maxDepth[n] = m_spharm[subj].maxDepth[n];
				if (m_minDepth[n] > m_spharm[subj].minDepth[n]) m_minDepth[n] = m_spharm[subj].minDepth[n];
			}
			else
			{
				cout << "\t" << subjDepth[nSubj * n + subj] << endl;
				FILE *fp = fopen(subjDepth[nSubj * n + subj], "r");
				char line[1024];
				fgets(line, sizeof(line), fp);
				fgets(line, sizeof(line), fp);
				fgets(line, sizeof(line), fp);
				//for (int i = 0; i < nDepth && !feof(fp); i++)
				for (int i = 0; i < nDepth; i++)
				{
					fscanf(fp, "%f", &m_spharm[subj].depth[nDepth * n + i]);
					//if (m_spharm[subj].depth[n * nProperties + i] > 0) m_spharm[subj].depth[n * nProperties + i] = 0;
				}
				fclose(fp);
				m_spharm[subj].meanDepth[n] = Statistics::mean(m_spharm[subj].depth, nDepth);
				m_spharm[subj].maxDepth[n] = Statistics::max(m_spharm[subj].depth, nDepth);
				m_spharm[subj].minDepth[n] = Statistics::min(m_spharm[subj].depth, nDepth);
				if (m_maxDepth[n] < m_spharm[subj].maxDepth[n]) m_maxDepth[n] = m_spharm[subj].maxDepth[n];
				if (m_minDepth[n] > m_spharm[subj].minDepth[n]) m_minDepth[n] = m_spharm[subj].minDepth[n];
			}
		}
		cout << "done\n";
	}

	for (int i = 0; i < m_sphere->nFace(); i++)
		if (m_entropy[i].occupied && m_entropy[i].samples.size() == nSubj) m_entropyList.push_back(i);

	// depth variance
	cout << "Generating sample points for depth variance.. ";
	icosahedron(3);
	cout << "done\n";
	
	// work space
	int depth = m_depthvar.size();
	m_cov_weight = new float[depth * nProperties];
	float w = 1.0f;
	for (int n = 0; n < nProperties; n++)
		for (int i = 0; i < depth; i++)
			m_cov_weight[n * depth + i] = 1.0f;

	m_cov_depth = new float[nSubj + 1];
	m_cov_eig = new float[(nSubj + 1) * (nSubj + 1)];
	m_pointList = new float[(nSubj + 1) * (depth * nProperties)];
	for (int n = 0; n < nProperties; n++)
	{
		for (int i = 0; i < depth; i++)
		{
			float coeffs[3];
			int id = m_tree->closestFace(m_depthvar[i], coeffs);
			m_pointList[nSubj * depth * nProperties + n * depth + i] = (depthInterpolation(&m_depth[nDepth * n], id, coeffs, m_sphere) - m_minDepth[n]) / (m_maxDepth[n] - m_minDepth[n]);
		}
	}
	
	// dpeth cache
	for (int subj = 0; subj < nSubj; subj++)
	{
		m_spharm[subj].depth_cache = new int[depth * nProperties];
		for (int n = 0; n < nProperties; n++)
			for (int i = 0; i < depth; i++)
				m_spharm[subj].depth_cache[depth * n + i] = -1;
	}
	
	cout << "Cache size: " << depth * nProperties << endl;
	cout << "Initialization done\n";
}

void GroupwiseRegistration::init(char *sphere, char **tmpDepth, char **subjDepth, char **coeff, char **correspondence, int nSubj, int deg, int nProperties)
{
	// unit sphere information
	cout << "Loading unit sphere information..\n";
	m_sphere = new Mesh();
	m_sphere->openFile(sphere);

	m_coeff_fn = coeff;
	m_spharm = new spharm[nSubj];
	m_entropy = new entropy[m_sphere->nFace()];

	int nDepth = m_sphere->nVertex();
	m_depth = new float[nDepth * nProperties];
	m_meanDepth = new float[nProperties];
	m_maxDepth = new float[nProperties];
	m_minDepth = new float[nProperties];

	for (int n = 0; n < nProperties; n++)
	{
		FILE *fp = fopen(tmpDepth[n], "r");
		for (int i = 0; i < nDepth && !feof(fp); i++)
		{
			fscanf(fp, "%f", &m_depth[nDepth * n + i]);
			//if (m_depth[i] > 0) m_depth[i] = 0;
		}
		fclose(fp);
		m_meanDepth[n] = Statistics::mean(m_depth, nDepth);
		m_maxDepth[n] = Statistics::max(m_depth, nDepth);
		m_minDepth[n] = Statistics::min(m_depth, nDepth);
	}

	//cout << "Depth Range: [" << m_minDepth << ", " << m_maxDepth << "]\n";

	for (int i = 0; i < m_sphere->nFace(); i++)
		m_entropy[i].occupied = false;

	m_tree = new AABB(m_sphere);

	cout << "Loading initial deformation..\n";
	m_csize = 0;
	for (int subj = 0; subj < nSubj; subj++)
	{
		int degree;
		int pole[3];
		FILE *fp = fopen(coeff[subj],"r");
		fscanf(fp, "%f %f %f", &pole[0], &pole[1], &pole[2]);
		fscanf(fp, "%d", &degree);
		//if (degree < deg) degree = deg;
		degree = deg;
		m_csize += (degree + 1) * (degree + 1);
		fclose(fp);
	}
	m_coeff = new float[m_csize * 2];
	memset(m_coeff, 0, sizeof(float) * m_csize * 2);

	for (int subj = 0; subj < nSubj; subj++)
	{
		cout << "suject " << subj << endl;
		// spherical harmonics information
		m_spharm[subj] = spharm();
		FILE *fp = fopen(coeff[subj],"r");
		fscanf(fp, "%f %f %f", &m_spharm[subj].pole[0], &m_spharm[subj].pole[1], &m_spharm[subj].pole[2]);
		fscanf(fp, "%d", &m_spharm[subj].degree);
		/*int degree = m_spharm[subj].degree;
		if (degree < deg) degree = deg;*/
		int degree = deg;
		if (m_spharm[subj].degree > deg) m_spharm[subj].degree = deg;

		int n = (degree + 1) * (degree + 1);
		m_spharm[subj].coeff = new float*[n * 2];
		for (int i = 0; i < n; i++)
		{
			m_spharm[subj].coeff[i] = &m_coeff[nSubj * 2 * i + subj * 2];
			m_spharm[subj].coeff[n + i] = &m_coeff[nSubj * 2 * i + subj * 2 + 1];
		}

		for (int i = 0; i < (m_spharm[subj].degree + 1) * (m_spharm[subj].degree + 1); i++)
			fscanf(fp, "%f %f", m_spharm[subj].coeff[i], m_spharm[subj].coeff[n + i]);
		fclose(fp);

		m_spharm[subj].degree = degree;

		// sphere information
		m_spharm[subj].sphere = new Mesh();
		m_spharm[subj].sphere->openFile(sphere);

		// vertices on the sphere
		for (int i = 0; i < m_spharm[subj].sphere->nVertex(); i++)
		{
			sphereVertex *p = new sphereVertex();
			Vertex *v = (Vertex *)m_spharm[subj].sphere->vertex(i);
			const float *v0 = v->fv();
			p->Y = new float[(degree + 1) * (degree + 1)];
			p->v[0] = v0[0]; p->v[1] = v0[1]; p->v[2] = v0[2];
			SphericalHarmonics::basis(degree, p->v, p->Y);
			m_spharm[subj].vertex.push_back(p);
		}

		// deformed irregular sphere
		m_updated = new bool[subj];
		cout << "Updating sphere deformation.. ";
		updateDeformation(subj);
		cout << "done\n";

		// AABB tree construction
		cout << "Constructing a search tree.. ";
		m_spharm[subj].tree = new AABB(m_spharm[subj].sphere);
		cout << "done\n";

		// depth information
		cout << "Loading depth information.. ";
		int nDepth = m_spharm[subj].sphere->nVertex();
		m_spharm[subj].meanDepth = new float[nProperties];
		m_spharm[subj].maxDepth = new float[nProperties];
		m_spharm[subj].minDepth = new float[nProperties];
		
		for (int n = 0; n < nProperties; n++)
		{
			m_spharm[subj].depth = new float[nDepth * nProperties];
			/*m_spharm[subj].depth_cache = new int[nDepth * nProperties];
			for (int i = 0; i < nDepth; i++) m_spharm[subj].depth_cache[n * nDepth + i] = -1;*/
			fp = fopen(subjDepth[nProperties * subj + n], "r");
			for (int i = 0; i < nDepth && !feof(fp); i++)
			{
				fscanf(fp, "%f", &m_spharm[subj].depth[n * nDepth + i]);
				//if (m_spharm[subj].depth[n * nProperties + i] > 0) m_spharm[subj].depth[n * nProperties + i] = 0;
			}
			fclose(fp);
			m_spharm[subj].meanDepth[n] = Statistics::mean(m_spharm[subj].depth, nDepth);
			m_spharm[subj].maxDepth[n] = Statistics::max(m_spharm[subj].depth, nDepth);
			m_spharm[subj].minDepth[n] = Statistics::min(m_spharm[subj].depth, nDepth);
			if (m_maxDepth[n] < m_spharm[subj].maxDepth[n]) m_maxDepth[n] = m_spharm[subj].maxDepth[n];
			if (m_minDepth[n] > m_spharm[subj].minDepth[n]) m_minDepth[n] = m_spharm[subj].minDepth[n];
		}
		cout << "done\n";

		// entropy computation
		cout << "Loading entropy samples.. ";
		FILE *fp1 = fopen(correspondence[subj * 2], "r");
		FILE *fp2 = fopen(correspondence[subj * 2 + 1], "r");
		while (!feof(fp1))
		{
			int n;
			if (fscanf(fp1, "%d", &n) == -1) break;
			for (int i = 0; i < n; i++)
			{
				// indices for corresponding points
				int srcid, tmpid;
				fscanf(fp1, "%d %d", &srcid, &tmpid);

				// correspondence information
				float src[3], tmp[3], id;
				fscanf(fp2, "%f %f %f %f %f %f %f", &src[0],&src[1],&src[2],&tmp[0],&tmp[1],&tmp[2],&id);
				float *Y = new float[(degree + 1) * (degree + 1)];
				point *p = new point();
				SphericalHarmonics::basis(degree, src, Y);
				p->p[0] = src[0]; p->p[1] = src[1]; p->p[2] = src[2];
				p->subject = subj;
				p->Y = Y;
				m_point.push_back(p);

				if (!m_entropy[tmpid].occupied)
				{
					m_entropy[tmpid].p[0] = tmp[0]; m_entropy[tmpid].p[1] = tmp[1]; m_entropy[tmpid].p[2] = tmp[2];
					m_entropy[tmpid].occupied = true;
				}
				// debug
				else assert(m_entropy[tmpid].p[0] == tmp[0] && m_entropy[tmpid].p[1] == tmp[1] && m_entropy[tmpid].p[2] == tmp[2]);

				m_entropy[tmpid].samples.push_back(p);
			}
		}
		fclose(fp1);
		fclose(fp2);
		cout << "done\n";
	}

	for (int i = 0; i < m_sphere->nFace(); i++)
		if (m_entropy[i].occupied && m_entropy[i].samples.size() == nSubj) m_entropyList.push_back(i);

	// depth variance
	cout << "Generating sample points for depth variance.. ";
	icosahedron(5);
	cout << "done\n";

	// work space
	int landmark = m_entropyList.size() * 3;
	//int depth = m_entropyList.size() * 1;
	int depth = m_depthvar.size();
	m_cov_weight = new float[landmark + depth * nProperties];
	float w = (float)(depth * nProperties) / (float)landmark;
	//float w = 1.0f;
	for (int i = 0; i < landmark; i++) m_cov_weight[i] = w;
	for (int i = 0; i < depth * nProperties; i++) m_cov_weight[landmark + i] = 1.0f;

	/*for (int i = 0; i < landmark; i++) m_cov_weight[i] = 1.0f;
	for (int i = 0; i < depth; i++) m_cov_weight[landmark + i] = 0;*/

	//cout << "Depth Range: [" << m_minDepth << ", " << m_maxDepth << "]\n";

	m_cov_depth = new float[nSubj + 1];
	//m_cov_eig = new float[(landmark + depth) * (landmark + depth)];
	m_cov_eig = new float[(nSubj + 1) * (nSubj + 1)];
	m_pointList = new float[(nSubj + 1) * (landmark + depth * nProperties)];
	for (int n = 0; n < nProperties; n++)
	{
		for (int i = 0; i < depth; i++)
		{
			float coeffs[3];
			int id = m_tree->closestFace(m_depthvar[i], coeffs);
			//m_pointList[nSubj * (landmark + depth) + landmm_nPropertiesark + i] = depthInterpolation(m_depth, id, coeffs, m_sphere);
			m_pointList[nSubj * (landmark + depth) + landmark + i] = (depthInterpolation(&m_depth[depth * n], id, coeffs, m_sphere) - m_minDepth[n]) / (m_maxDepth[n] - m_minDepth[n]);
		}
	}
	/*for (int i = 0; i < depth; i++)
	{
		float coeff[3];
		int id = m_tree->closestFace(m_entropy[m_entropyList[i]].p, coeff);
		m_pointList[nSubj * (landmark + depth) + landmark + i] = depthInterpolation(m_depth, id, coeff, m_sphere);
	}*/
	/*for (int i = 0; i < m_entropyList.size(); i++)
	{
		int id = m_entropyList[i];
		memcpy(&m_pointList[nSubj * m_entropyList.size() * 3 + i * 3], m_entropy[id].p, sizeof(float) * 3);
	}*/
	
	// dpeth cache
	for (int subj = 0; subj < nSubj; subj++)
	{
		m_spharm[subj].depth_cache = new int[depth * nProperties];
		for (int n = 0; n < nProperties; n++)
			for (int i = 0; i < depth; i++)
				m_spharm[subj].depth_cache[depth * n + i] = -1;
	}
}

void GroupwiseRegistration::reconsCoord(const float *v0, float *v1, float *Y, float **coeff, float degree, float *pole)
{
	// spharm basis
	int n = (degree + 1) * (degree + 1);

	float mat[9];
	memset(mat, 0, sizeof(float) * 9);
	Coordinate::rotation2equator((float *)v0, mat, pole);

	// rotation to the eqautor
	float rv[3];
	Coordinate::rotPoint(v0, mat, rv);
	
	// polar coodinate
	float phi, theta;
	Coordinate::cart2sph(rv, &phi, &theta);

	// displacement
	float delta[2] = {0, 0};
	for (int i = 0; i < n; i++)
	{
		delta[0] += Y[i] * *coeff[i];
		delta[1] += Y[i] * *coeff[n + i];
	}
	phi += delta[0]; theta += delta[1];
	
	Coordinate::sph2cart(phi, theta, rv);

	// inverse rotation
	Coordinate::rotPointInv(rv, mat, v1);
}

void GroupwiseRegistration::updateDeformation(int subject)
{
	for (int i = 0; i < m_spharm[subject].vertex.size(); i++)
	{
		float v1[3];
		reconsCoord(m_spharm[subject].vertex[i]->v, v1, m_spharm[subject].vertex[i]->Y, m_spharm[subject].coeff, m_spharm[subject].degree, m_spharm[subject].pole);
		Vertex *v = (Vertex *)m_spharm[subject].sphere->vertex(i);
		MathVector V(v1); V.unit();
		v->setVertex(V.fv());

		// debug
		/*float *v0 = m_spharm[subject].vertex[i]->v;
		for (int j = 0; j < 3; j++)
			cout << v0[j] << " ";
		cout << " -> ";
		for (int j = 0; j < 3; j++)
			cout << v1[j] << " ";
		cout << endl;*/
	}
	m_updated[subject] = false;
}

float GroupwiseRegistration::landmarkEntropyMulti(void)
{
	int d = m_depthvar.size();
	int nSubj = m_nSubj + 1;
	float E = 0;
	float *cov = m_cov_eig;	memset(cov, 0, sizeof(float) * nSubj * nSubj);
	float *p = m_pointList;
	float m[3];
	float *newp = new float[3];
	
	// depth information
	int fid;
	float err = 0;
	for (int n = 0; n < m_nProperties; n++)
	{
		for (int i = 0; i < d; i++)
		{
			float coeff[3];
			for (int j = 0; j < m_nSubj; j++)
			{
				if (m_spharm[j].depth_cache[d * n + i] != -1)
				{
					Face f = *m_spharm[j].sphere->face(m_spharm[j].depth_cache[d * n + i]);
					Vertex a = *f.vertex(0);
					Vertex b = *f.vertex(1);
					Vertex c = *f.vertex(2);

					// bary centric
					Coordinate::cart2bary((float *)a.fv(), (float *)b.fv(), (float *)c.fv(), m_depthvar[i], coeff);

					if (coeff[0] >= err && coeff[1] >= err && coeff[2] >= err)
					{
						fid = m_spharm[j].depth_cache[d * n + i];
					}
					else
					{
						if (!m_updated[j])
						{
							m_updated[j] = true;
							m_spharm[j].tree->update();
						}
						fid = m_spharm[j].tree->closestFace(m_depthvar[i], coeff);
					}
				}
				else
				{
					if (!m_updated[j])
					{
						m_updated[j] = true;
						m_spharm[j].tree->update();
					}
					fid = m_spharm[j].tree->closestFace(m_depthvar[i], coeff);
				}
				//p[j * (n * 3 + d * 1) + n * 3 + i] = depthInterpolation(m_spharm[j].depth, fid, coeff, m_spharm[j].sphere);
				p[j * (d * m_nProperties) + d * n + i] = (depthInterpolation(&m_spharm[j].depth[m_spharm[j].vertex.size() * n], fid, coeff, m_spharm[j].sphere) - m_minDepth[n]) / (m_maxDepth[n] - m_minDepth[n]);
				//p[j * (d * m_nProperties) + d * n + i] = depthInterpolation(&m_spharm[j].depth[m_spharm[j].vertex.size() * n], fid, coeff, m_spharm[j].sphere);

				m_spharm[j].depth_cache[d * n + i] = fid;
			}
		}
	}

	Statistics::wcov_trans(p, nSubj, d * m_nProperties, cov, m_cov_weight);
	//Statistics::cov_trans(p, nSubj, d * m_nProperties, cov);

	// entropy
	float *eig = new float[nSubj];
	eigenvalues(cov, nSubj, eig);
	for (int i = 1; i < nSubj; i++)
		E += log(eig[i]);

	return E;
}

float GroupwiseRegistration::landmarkEntropy(void)
{
	int n = m_entropyList.size();
	int d = m_depthvar.size();
	int nSubj = m_nSubj + 1;
	float E = 0;
	float *cov = m_cov_eig;	memset(cov, 0, sizeof(float) * nSubj * nSubj);
	float *p = m_pointList;
	float m[3];
	float *newp = new float[3];

	//FILE *fp=fopen("a.txt","w");
	for (int i = 0; i < n; i++)
	{
		int id = m_entropyList[i];
		int nSamples = m_entropy[id].samples.size();
		// point on the template
		memcpy(&m_pointList[nSamples * (n * 3 + d * m_nProperties) + i * 3], m_entropy[id].p, sizeof(float) * 3);
		memcpy(m, m_entropy[id].p, sizeof(float) * 3);
		// deformed points
		for (int j = 0; j < nSamples; j++)
		{
			int subj = m_entropy[id].samples[j]->subject;
			reconsCoord(m_entropy[id].samples[j]->p, &p[subj * (n * 3 + d * m_nProperties) + i * 3], m_entropy[id].samples[j]->Y, m_spharm[subj].coeff, m_spharm[subj].degree, m_spharm[subj].pole);

			// mean
			for (int k = 0; k < 3; k++) m[k] += p[subj * (n * 3 + d * m_nProperties) + i * 3 + k];
		}

		// forcing the mean to be on the sphere
		float norm = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
		for (int k = 0; k < 3; k++) m[k] /= norm;

		// projection
		for (int j = 0; j <= nSamples; j++)
		{
			Coordinate::proj2plane(m[0], m[1], m[2], -1, &p[j * (n * 3 + d * 1) + i * 3], newp);
			memcpy(&p[j * (n * 3 + d * m_nProperties) + i * 3], newp, sizeof(float) * 3);
			//fprintf(fp, "%f %f %f %f %f %f %f\n", m[0], m[1], m[2], newp[0],newp[1],newp[2],(float)j);
		}
	}

	// depth information
	int fid;
	float err = 0;
	for (int k = 0; k < m_nProperties; k++)
	{
		for (int i = 0; i < d; i++)
		{
			float coeff[3];
			for (int j = 0; j < m_nSubj; j++)
			{
				if (m_spharm[j].depth_cache[d * n + i] != -1)
				{
					Face f = *m_spharm[j].sphere->face(m_spharm[j].depth_cache[d * k + i]);
					Vertex a = *f.vertex(0);
					Vertex b = *f.vertex(1);
					Vertex c = *f.vertex(2);

					// bary centric
					Coordinate::cart2bary((float *)a.fv(), (float *)b.fv(), (float *)c.fv(), m_depthvar[i], coeff);

					if (coeff[0] >= err && coeff[1] >= err && coeff[2] >= err)
					{
						fid = m_spharm[j].depth_cache[d * n + i];
					}
					else
					{
						if (!m_updated[j])
						{
							m_updated[j] = true;
							m_spharm[j].tree->update();
						}
						fid = m_spharm[j].tree->closestFace(m_depthvar[i], coeff);
					}
				}
				else
				{
					if (!m_updated[j])
					{
						m_updated[j] = true;
						m_spharm[j].tree->update();
					}
					fid = m_spharm[j].tree->closestFace(m_depthvar[i], coeff);
				}
				//p[j * (n * 3 + d * 1) + n * 3 + i] = depthInterpolation(m_spharm[j].depth, fid, coeff, m_spharm[j].sphere);
				p[j * (n * 3 + d * m_nProperties) + n * 3 + d * k + i] = (depthInterpolation(&m_spharm[j].depth[d * k], fid, coeff, m_spharm[j].sphere) - m_minDepth[k]) / (m_maxDepth[k] - m_minDepth[k]);

				m_spharm[j].depth_cache[d * k + i] = fid;
			}
		}
	}

	Statistics::wcov_trans(p, nSubj, n * 3 + d * m_nProperties, cov, m_cov_weight);
	//Statistics::cov_trans(p, nSubj, n * 3 + d * m_nProperties, cov);

	// entropy
	float *eig = new float[nSubj];
	eigenvalues(cov, nSubj, eig);
	for (int i = 1; i < nSubj; i++)
		E += log(eig[i]);

	return E;
}

float GroupwiseRegistration::landmarkEntropyMedian(void)
{
	int n = m_entropyList.size();
	int d = m_depthvar.size();
	int nSubj = m_nSubj + 1;
	float E = 0;
	float *cov = m_cov_eig;	memset(cov, 0, sizeof(float) * nSubj * nSubj);
	float *p = m_pointList;
	float m[3];
	float *newp = new float[3];

	// scaling
	float d1 = 0.03 * 0.01;
    float d2 = 16.8 * 0.01;
	float u = (d2 + d1) / 2;
	float sigma = ((d2 + d1) / 2 - d1) / 3;

	//FILE *fp=fopen("lm.txt","w");
	for (int i = 0; i < n; i++)
	{
		int id = m_entropyList[i];
		int nSamples = m_entropy[id].samples.size();
		// point on the template
		memcpy(&m_pointList[nSamples * (n * 3 + d * m_nProperties) + i * 3], m_entropy[id].p, sizeof(float) * 3);
		//memcpy(m, m_entropy[id].p, sizeof(float) * 3);

		// deformed points
		float *x = new float[nSamples];
		float *y = new float[nSamples];
		float *z = new float[nSamples];

		for (int j = 0; j < nSamples; j++)
		{
			int subj = m_entropy[id].samples[j]->subject;
			reconsCoord(m_entropy[id].samples[j]->p, &p[subj * (n * 3 + d * m_nProperties) + i * 3], m_entropy[id].samples[j]->Y, m_spharm[subj].coeff, m_spharm[subj].degree, m_spharm[subj].pole);

			// median
			x[j] = p[subj * (n * 3 + d * m_nProperties) + i * 3 + 0];
			y[j] = p[subj * (n * 3 + d * m_nProperties) + i * 3 + 1];
			z[j] = p[subj * (n * 3 + d * m_nProperties) + i * 3 + 2];
		}

		m[0] = Statistics::median(x, nSamples);
		m[1] = Statistics::median(y, nSamples);
		m[2] = Statistics::median(z, nSamples);

		delete [] x;
		delete [] y;
		delete [] z;

		// forcing the mean to be on the sphere
		float norm = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
		for (int k = 0; k < 3; k++) m[k] /= norm;

		// projection
		for (int j = 0; j <= nSamples; j++)
		{
			Coordinate::proj2plane(m[0], m[1], m[2], -1, &p[j * (n * 3 + d * m_nProperties) + i * 3], newp);

			// scaling
			float len = sqrt((newp[0] - m[0]) * (newp[0] - m[0]) + (newp[1] - m[1]) * (newp[1] - m[1]) + (newp[2] - m[2]) * (newp[2] - m[2]));
			//float slen = Statistics::normal_cdf(len, u, sigma);
			float slen = (len < d2) ? len: d2;
			float ratio = 1.0f;
			if (len > 0) ratio = slen / len;

			newp[0] = m[0] + (newp[0] - m[0]) * ratio;
			newp[1] = m[1] + (newp[1] - m[1]) * ratio;
			newp[2] = m[2] + (newp[2] - m[2]) * ratio;

			memcpy(&p[j * (n * 3 + d * m_nProperties) + i * 3], newp, sizeof(float) * 3);
			//fprintf(fp, "%f %f %f %f %f %f %f\n", m[0], m[1], m[2], newp[0],newp[1],newp[2],(float)j);
			/*fprintf(fp, "%f ", p[j * (n * 3 + d * 1) + i * 3]);
			fprintf(fp, "%f ", p[j * (n * 3 + d * 1) + i * 3+1]);
			fprintf(fp, "%f ", p[j * (n * 3 + d * 1) + i * 3+2]);*/
		}
		//fprintf(fp, "\n");
	}
	//fclose(fp);

	// depth information
	//fp=fopen("depth.txt","w");
	int fid;
	float err = 0;
	for (int k = 0; k < m_nProperties; k++)
	{
		for (int i = 0; i < d; i++)
		{
			float coeff[3];
			for (int j = 0; j < m_nSubj; j++)
			{
				if (m_spharm[j].depth_cache[d * n + i] != -1)
				{
					Face f = *m_spharm[j].sphere->face(m_spharm[j].depth_cache[d * k + i]);
					Vertex a = *f.vertex(0);
					Vertex b = *f.vertex(1);
					Vertex c = *f.vertex(2);

					// bary centric
					Coordinate::cart2bary((float *)a.fv(), (float *)b.fv(), (float *)c.fv(), m_depthvar[i], coeff);

					if (coeff[0] >= err && coeff[1] >= err && coeff[2] >= err)
					{
						fid = m_spharm[j].depth_cache[d * n + i];
					}
					else
					{
						if (!m_updated[j])
						{
							m_updated[j] = true;
							m_spharm[j].tree->update();
						}
						fid = m_spharm[j].tree->closestFace(m_depthvar[i], coeff);
					}
				}
				else
				{
					if (!m_updated[j])
					{
						m_updated[j] = true;
						m_spharm[j].tree->update();
					}
					fid = m_spharm[j].tree->closestFace(m_depthvar[i], coeff);
				}
				//p[j * (n * 3 + d * 1) + n * 3 + i] = depthInterpolation(m_spharm[j].depth, fid, coeff, m_spharm[j].sphere);
				p[j * (n * 3 + d * m_nProperties) + n * 3 + d * k + i] = (depthInterpolation(&m_spharm[j].depth[d * k], fid, coeff, m_spharm[j].sphere) - m_minDepth[k]) / (m_maxDepth[k] - m_minDepth[k]);
				//fprintf(fp, "%f(%f %f %f) ", p[j * (n * 3 + d * 1) + n * 3 + i], coeff[0], coeff[1], coeff[2]);

				m_spharm[j].depth_cache[d * k + i] = fid;
			}
			//fprintf(fp, "\n");
		}
		//fclose(fp);
	}

	Statistics::wcov_trans(p, nSubj, n * 3 + d * m_nProperties, cov, m_cov_weight);

	//fp = fopen("col.txt", "w");
	//for (int i = 0; i < n * 3 + d * 1; i++)
	//{
	//	//fprintf(fp, "%d: ", i);
	//	for (int j = 0; j < nSubj; j++)
	//		fprintf(fp, "%f ", p[i + (n * 3 + d * 1) * j]);
	//	fprintf(fp, "\n");
	//}
	//fclose(fp);

	/*fp = fopen("w.txt", "w");
	for (int i = 0; i < n+d; i++)
	{
		fprintf(fp, "%f\n", m_cov_weight[i]);
	}
	fclose(fp);*/

	/*FILE *fp = fopen("wcov.txt", "w");
	for (int i = 0; i < nSubj; i++)
	{
		for (int j = 0; j < nSubj; j++)
			fprintf(fp, "%f ", cov[i * nSubj + j]);
		fprintf(fp, "\n");
	}
	fclose(fp);*/

	/*float *mm = new float[n*3+d*1];
	
	memset(mm, 0, sizeof(float) * (n*3+d*1));
	Statistics::sum(p, nSubj, n*3+d*1, mm);

	fp = fopen("sum.txt", "w");
	for (int i=0;i<n*3+d*1;i++)
		fprintf(fp, "%f\n", mm[i]);
	fclose(fp);

	memset(mm, 0, sizeof(float) * (n*3+d*1));
	Statistics::mean(p, nSubj, n*3+d*1, mm);

	fp = fopen("mean.txt", "w");
	for (int i=0;i<n*3+d*1;i++)
		fprintf(fp, "%f\n", mm[i]);
	fclose(fp);*/

	// entropy
	float *eig = new float[nSubj];
	eigenvalues(cov, nSubj, eig);
	for (int i = 1; i < nSubj; i++)
		if (eig[i] > 0) E += log(eig[i]);

	return E;
}

void GroupwiseRegistration::eigenvalues(float *M, int dim, float *eig)
{
	int n = dim;
	int lwork = dim * 3 - 1;	// dimension of the work array
	int lda = n;			// lda: leading dimension
	int info;				// information (0 for successful exit)
	float *work = new float[lwork];

	ssyev_("N", "L", &n, M, &lda, eig, work, &lwork, &info);

	/*int neg = 0;
	int nnz = 0;
	for (int i = 0; i < dim; i++)
	{
		if (fabs(eig[i]) > 0)
		{
			nnz++;
			if (eig[i] < 0)
				neg++;
		}
	}
	if (nnz == neg)
		for (int i = 0; i < dim; i++)
			eig[i] *= -1;*/

	delete [] work;
}

float GroupwiseRegistration::depthVariance(void)
{
	float coeff[3];
	float var = 0;

	for (int i = 0; i < m_depthvar.size(); i++)
	{
		//m_cov_depth[0] = m_depth[i] - m_meanDepth;
		m_cov_depth[0] = (m_depth[i] - m_minDepth[0]) / (m_maxDepth[0] - m_minDepth[0]);
		for (int j = 0; j < m_nSubj; j++)
		{
			int id = m_spharm[j].tree->closestFace(m_depthvar[i], coeff);
			float depth = depthInterpolation(&m_spharm[j].depth[0], id, coeff, m_spharm[j].sphere);
			//m_cov_depth[j + 1] = depth - m_spharm[j].meanDepth;
			m_cov_depth[j + 1] = (depth - m_minDepth[0]) / (m_maxDepth[0] - m_minDepth[0]);
			//cout << "Qry: "<< m_depthvar[i][0] << ", " << m_depthvar[i][1] << ", " << m_depthvar[i][2] << " => ";
			//cout << id << " (" << coeff[0] << ", " << coeff[1] << ", " << coeff[2] << ")" << endl;
		}
		var += Statistics::var(m_cov_depth, m_nSubj + 1);
	}

	return var;
}

float GroupwiseRegistration::depthInterpolation(float *refMap, int index, float *coeff, Mesh *mesh)
{
	float depth = 0;
	float err = 0; // numerical error
	bool success = false;

	Face f = *mesh->face(index);
	Vertex a = *f.vertex(0);
	Vertex b = *f.vertex(1);
	Vertex c = *f.vertex(2);
	depth = refMap[a.id()] * coeff[0] + refMap[b.id()] * coeff[1] + refMap[c.id()] * coeff[2];

	return depth;
}

float GroupwiseRegistration::cost(float *coeff, int statusStep)
{
	for (int i = 0; i < m_nSubj; i++)
	{
		updateDeformation(i);
		//m_spharm[i].tree->update();
	}
	//float cost = landmarkEntropy();
	//float cost = landmarkEntropyMedian();
	float cost = landmarkEntropyMulti();
	if (nIter % statusStep == 0)
	{
		cout << "Cost: " << cost << endl;
		//printf("%0#8d\t%f\n", nIter, cost);
		if (m_clfp != NULL)
		{
			fprintf(m_clfp, "%0#8d\t%f\t", nIter, cost);
			//for (int id = 0; id < m_nSubj; id++)
			//{
			//	fprintf(m_clfp, "%f\t%f\t", *m_spharm[id].coeff[0], *m_spharm[id].coeff[1]);
			//}
			//for (int id = 0; id < m_nSubj; id++)
			//{
			//	//int n = (m_spharm[id].degree + 1) * (m_spharm[id].degree + 1);
			//	int n = 36;
			//	fprintf(m_clfp, "%f\t%f\t", *m_spharm[id].coeff[n], *m_spharm[id].coeff[n + 1]);
			//}
			/*fprintf(m_clfp, "\n");
			for (int i = 0; i < 36; i++)
			{
				fprintf(m_clfp, "%f\t", *m_spharm[0].coeff[i]);
			}
			fprintf(m_clfp, "\n");
			for (int i = 0; i < 36; i++)
			{
				fprintf(m_clfp, "%f\t", *m_spharm[0].coeff[36 + i]);
			}
			fprintf(m_clfp, "\n");*/
			fprintf(m_clfp, "\n");
			fflush(m_clfp);

			saveLDeformation("group.txt");
			for (int i = 0; i < 42; i++)
			{
				char fn[1024];
				sprintf(fn, "%s", &m_coeff_fn[i][29]);
				saveLCoeff(fn, i);
			}
		}
		//fflush(stdout);
	}
	nIter++;

	return cost;
}

void GroupwiseRegistration::optimization(void)
{
	cost_function costFunc(this);
	//cout << "var: " << depthVariance() << endl;
	int prev = 0;
	int deg = 2;
	int step = 1;
	
	int n1 = (deg + 1) * (deg + 1) * m_nSubj * 2;
	int n2 = m_csize * 2 - n1;
	//memset(&m_coeff[n1], 0, n2 * sizeof(float));

	while (deg < m_spharm[0].degree)
	{
		nIter = 0;
		int n = (deg + 1) * (deg + 1) * m_nSubj * 2 - prev;
		min_newuoa(n, &m_coeff[prev], costFunc, 0.01f, 1e-6f, 30000);
		prev = (deg + 1) * (deg + 1) * m_nSubj * 2;
		deg = min(deg + step, m_spharm[0].degree);
	}
	
	nIter = 0;
	min_newuoa(m_csize * 2, m_coeff, costFunc, 0.001f, 1e-6f, 30000);
}

int GroupwiseRegistration::icosahedron(int degree)
{
	// http://www.1activeserverpagesstreet.com/vb/scripts/ShowCode.asp?txtCodeId=9814&lngWId=3
	vector<MathVector *> triangles;
	vector<MathVector> vertices;

	float t = (1 + sqrt(5.0)) / 2.0;
    float s = sqrt(1 + t * t);

    // create the 12 vertices
    MathVector v0 = MathVector(t, 1.0, 0.0) / s;
    MathVector v1 = MathVector(-t, 1.0, 0.0) / s;
    MathVector v2 = MathVector(t, -1.0, 0.0) / s;
    MathVector v3 = MathVector(-t, -1.0, 0.0) / s;
    MathVector v4 = MathVector(1.0, 0.0, t) / s;
    MathVector v5 = MathVector(1.0, 0.0, -t) / s;
    MathVector v6 = MathVector(-1.0, 0.0, t) / s;
    MathVector v7 = MathVector(-1.0, 0.0, -t) / s;
    MathVector v8 = MathVector(0.0, t, 1.0) / s;
    MathVector v9 = MathVector(0.0, -t, 1.0) / s;
    MathVector v10 = MathVector(0.0, t, -1.0) / s;
    MathVector v11 = MathVector(0.0, -t, -1.0) / s;
    
    // create the 20 triangles
	MathVector *f; 
	f = new MathVector[3]; f[0] = v0; f[1] = v8; f[2] = v4; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v1; f[1] = v10; f[2] = v7; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v2; f[1] = v9; f[2] = v11; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v7; f[1] = v3; f[2] = v1; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v0; f[1] = v5; f[2] = v10; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v3; f[1] = v9; f[2] = v6; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v3; f[1] = v11; f[2] = v9; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v8; f[1] = v6; f[2] = v4; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v2; f[1] = v4; f[2] = v9; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v3; f[1] = v7; f[2] = v11; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v4; f[1] = v2; f[2] = v0; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v9; f[1] = v4; f[2] = v6; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v2; f[1] = v11; f[2] = v5; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v0; f[1] = v10; f[2] = v8; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v5; f[1] = v0; f[2] = v2; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v10; f[1] = v5; f[2] = v7; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v1; f[1] = v6; f[2] = v8; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v1; f[1] = v8; f[2] = v10; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v6; f[1] = v1; f[2] = v3; triangles.push_back(f);
	f = new MathVector[3]; f[0] = v11; f[1] = v7; f[2] = v5; triangles.push_back(f);

	// subdivision
	for (int d = 0; d < degree; d++)
	{
		int nFaces = triangles.size();
		for (int i = 0 ; i < nFaces; i++)
		{
			MathVector *f = triangles[i];
			MathVector a = f[0], b = f[1], c = f[2];
			MathVector v1 = a + b;
			MathVector v2 = c + a;
			MathVector v3 = b + c;
			// normalization
			v1.unit(); v2.unit(); v3.unit();
			f[0] = v1; f[1] = v3; f[2] = v2; // overwrite the original
			/*MathVector f1[3] = {a, v1, v2}; triangles.push_back(f1);
			MathVector f2[3] = {c, v2, v3}; triangles.push_back(f2);
			MathVector f3[3] = {b, v3, v1}; triangles.push_back(f3);*/
			MathVector *f1 = new MathVector[3]; f1[0] = a; f1[1] = v1; f1[2] = v2; triangles.push_back(f1);
			MathVector *f2 = new MathVector[3]; f2[0] = c; f2[1] = v2; f2[2] = v3; triangles.push_back(f2);
			MathVector *f3 = new MathVector[3]; f3[0] = b; f3[1] = v3; f3[2] = v1; triangles.push_back(f3);
		}
		// debug
		//printf("phase #%d - %d\n", d, triangles.size());
		////for (int i = triangles.size() - 10; i < triangles.size(); i++)
		////for (int i = 0 ; i < triangles.size(); i++)
		//for (int i = 0 ; i < 10; i++)
		//{
		//	printf("[%0#2d] ", i);
		//	MathVector *f = triangles[i];
		//	for (int j = 0; j < 3; j++)
		//		printf("[%.3f %.3f %.3f] ", f[j][0], f[j][1], f[j][2]);
		//	printf("\n");
		//}
	}
	for (int i = 0; i < triangles.size(); i++)
	{
		MathVector *f = triangles[i];
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
		m_depthvar.push_back(p);
	}

	return m_depthvar.size();
}

void GroupwiseRegistration::saveLDeformation(char *filename)
{
	//float *q = m_entropy[m_entropyList[0]].p;
	FILE *fp = fopen(filename, "w");
	for (int i = 0; i < m_entropyList.size(); i++)
	{
		int id = m_entropyList[i];
		int nSamples = m_entropy[id].samples.size();
		for (int j = 0; j < nSamples; j++)
		{
			float v[3], *p = m_entropy[id].p;
			int subj = m_entropy[id].samples[j]->subject;
			reconsCoord(m_entropy[id].samples[j]->p, v, m_entropy[id].samples[j]->Y, m_spharm[subj].coeff, m_spharm[subj].degree, m_spharm[subj].pole);
			fprintf(fp, "%f %f %f %f %f %f %f\n", p[0], p[1], p[2], v[0], v[1], v[2], (float)(m_entropy[id].samples[j]->subject + 1));
			//float v2[3];Coordinate::proj2plane(q[0],q[1],q[2],-1,v,v2);
			//fprintf(fp, "%f %f %f %f %f %f %f\n", p[0], p[1], p[2], v2[0], v2[1], v2[2], (float)(m_entropy[id].samples[j]->subject + 1));
		}
	}
	fclose(fp);
}

void GroupwiseRegistration::saveLCoeff(char *filename, int id)
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
