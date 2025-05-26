#pragma once

#include <iostream>
#include "Polyhedron.hpp"

using namespace std;

namespace PolyhedronLibrary
{
/// Import the polygonal mesh and test if the mesh is correct
/// mesh: a PolygonalMesh struct
/// return the result of the reading true if is success, false otherwise
bool ImportPolyhedronMesh(PolyhedronMesh& polyhedron, const string& InputFileDirectory);

bool ImportCell0Ds(PolyhedronMesh& polyhedron, const string& InputFileDirectory);
bool ImportCell1Ds(PolyhedronMesh& polyhedron, const string& InputFileDirectory);
bool ImportCell2Ds(PolyhedronMesh& polyhedron, const string& InputFileDirectory);


// Importazione da file
bool ImportPolyhedronMesh(PolyhedronMesh& polyhedron, const std::string& InputFileDirectory);

// Triangolazione geodetica classe I
bool TriangolaClasseI(int p, int q, int b,
                      const PolyhedronMesh& mesh_input,
                      PolyhedronMesh& mesh_output,int c);
					  // Triangolazione geodetica classe I
bool TriangolaClasse2(int p, int q, int b_input,
                      const PolyhedronMesh& mesh_input,
                      PolyhedronMesh& mesh_output,int c_input);

// Esportazione dei file .txt
bool EsportaMeshSuFile(const PolyhedronMesh& mesh, const std::string& outputDirectory);

// Funzione per trovare un vertice giÃ  presente (con tolleranza)
int TrovaVertice(const Eigen::Vector3d& v, const std::vector<Eigen::Vector3d>& lista, double toll = 1e-8);
//funzione che mi genera esagoni da cui prendo facce per triangolazione b=c_input
void GeneraEsagono(
    const std::vector<Vector3d>& nuovi_vertici,
    std::vector<std::vector<int>>& nuove_facce,
    int id_1, int id_2, int id_3, int id_4, int id_5, int id_6,
    int& id_baricentro,
    std::vector<Vector3d>& vertici_out);
}