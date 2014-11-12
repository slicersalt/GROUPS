#include "VtkIO.h"
#include "fstream"
#include "string.h"
#include "stdlib.h"

using namespace std;

VtkIO::VtkIO(void)
{
}

VtkIO::VtkIO(const char *filename): MeshIO()
{
	read(filename);
}

VtkIO::~VtkIO(void)
{
}

void VtkIO::read(const char *filename)
{
	int nVertex = 0;
	int nFace = 0;
	int nNormal = 0;
	char buf[255];
	
	ifstream fin(filename);
	while (!fin.eof())
	{
		fin.getline(buf, sizeof(buf));
		if (strlen(buf) > 0)
		{
			char *tokens;
			char *ptr = strtok(buf, " ");
			if (!strcasecmp(ptr, "points"))
			{
				ptr = strtok(tokens, " ");
				nVertex = atoi(ptr);
			}
			else if (!strcasecmp(ptr, "polygons"))
			{
				ptr = strtok(tokens, " ");
				nFace = atoi(ptr);
			}
			else if (!strcasecmp(ptr, "point_data"))
			{
				ptr = strtok(tokens, " ");
				nNormal = atoi(ptr);
			}
		}
	}
	if (nNormal != nVertex) nNormal = nVertex;
	else m_hasNormal = true;

	m_nVertex = nVertex;
	m_nFace = nFace;
	m_nNormal = nNormal;

	initArray();

	nVertex = 0;
	nFace = 0;
	nNormal = 0;

	fin.clear();
	fin.seekg(0, ios::beg);
	while (!fin.eof())
	{
		fin.getline(buf, sizeof(buf));
		if (strlen(buf) > 0)
		{
			char *tokens;
			char *ptr = strtok(buf, " ");
			if (!strcasecmp(ptr, "points"))
			{
				while (nVertex < m_nVertex * 3)
				{
					fin.getline(buf, sizeof(buf));
					ptr = strtok(buf, " \r\n");
					while (ptr != NULL)
					{
						m_vertex[nVertex++] = atof(ptr);
						ptr = strtok(tokens, " \r\n");
					}
				}
			}
			else if (!strcasecmp(ptr, "polygons"))
			{
				for (nFace = 0; nFace < m_nFace; nFace++)
				{
					fin.getline(buf, sizeof(buf));
					if (m_hasNormal)
					{
						sscanf(buf, "%*[0-9] %d %d %d", &m_face[nFace * 6], &m_face[nFace * 6 + 1], &m_face[nFace * 6 + 2]);
						m_face[nFace * 6 + 3] = m_face[nFace * 6];
						m_face[nFace * 6 + 4] = m_face[nFace * 6 + 1];
						m_face[nFace * 6 + 5] = m_face[nFace * 6 + 2];
					}
					else
					{
						sscanf(buf, "%*[0-9] %d %d %d", &m_face[nFace * 3], &m_face[nFace * 3 + 1], &m_face[nFace * 3 + 2]);
					}
				}
			}
			else if (!strcasecmp(ptr, "NORMALS"))
			{
				while (nNormal < m_nNormal * 3)
				{
					fin.getline(buf, sizeof(buf));
					ptr = strtok(buf, " ");
					while (ptr != NULL)
					{
						m_normal[nNormal++] = atof(ptr);
						ptr = strtok(tokens, " ");
					}
				}
			}
		}
	}
	
	fin.close();
}

void VtkIO::save(const char *filename, Mesh *mesh, bool normal, bool binary)
{
	ofstream fout;
	if (!binary) fout.open(filename, ios::out);
	else fout.open(filename, ios::out | ios::binary);
	fout << "# vtk DataFile Version 3.0" << endl;
	fout << "vtk_output" << endl;
	if (!binary) fout << "ASCII" << endl;
	else fout << "BINARY" << endl;
	fout << "DATASET POLYDATA" << endl;

	fout << "POINTS " << mesh->nVertex() << " float" << endl;
	for (int i = 0; i < mesh->nVertex(); i++)
	{
		Vertex v = *mesh->vertex(i);
		if (binary)
		{
			fout.write((char *)v.fv(), sizeof(float) * 3);
		}
		else
		{
			fout << v.fv()[0] << " " << v.fv()[1] << " " << v.fv()[2];
			fout << endl;
		}
	}
	fout << endl;

	fout << "POLYGONS " << mesh->nFace() << " " << mesh->nFace() * 4 << endl;
	for (int i = 0; i < mesh->nFace(); i++)
	{
		Face f = *mesh->face(i);
		if (binary)
		{
			unsigned char list = 3;
			fout.write((char *)&list, sizeof(unsigned char));
			Face f = *mesh->face(i);
			fout.write((char *)f.list(), sizeof(int) * 3);
		}
		else
		{
			fout << "3 " << f.list()[0] << " " << f.list()[1] << " " << f.list()[2] << endl;
		}
	}
	fout << endl;

	fout << "POINT_DATA " << mesh->nVertex() << endl;
	fout << "NORMALS normals float" << endl;
	for (int i = 0; i < mesh->nVertex(); i++)
	{
		Vertex v = *mesh->vertex(i);
		if (binary)
		{
			Normal vn = *mesh->normal(i);
			fout.write((char *)vn.fv(), sizeof(float) * 3);
		}
		else
		{
			Normal vn = *mesh->normal(i);
			fout << vn.fv()[0] << " " << vn.fv()[1] << " " << vn.fv()[2];
			fout << endl;
		}
	}
	fout << endl;

	fout.close();
}
