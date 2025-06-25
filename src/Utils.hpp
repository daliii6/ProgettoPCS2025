#pragma once

#include <iostream>
#include "Polyhedron.hpp"

using namespace std;

namespace PolyhedronLibrary
{
// Legge le celle 0D (nodi) da file e le inserisce nella mesh.
bool ImportCell0Ds(PolyhedronMesh& polyhedron, const string& InputFileDirectory);

// Legge le celle 1D (spigoli) da file e le inserisce nella mesh.
bool ImportCell1Ds(PolyhedronMesh& polyhedron, const string& InputFileDirectory);

// Legge le celle 2D (facce) da file e le inserisce nella mesh.
bool ImportCell2Ds(PolyhedronMesh& polyhedron, const string& InputFileDirectory);

// Importa una mesh poliedrale da file e verifica che sia corretta.
// Restituisce true se la lettura ha successo, false altrimenti.
bool ImportPolyhedronMesh(PolyhedronMesh& polyhedron, const string& InputFileDirectory);

// Costruisce la lista di adiacenza della mesh.
// Ogni nodo ha una lista di nodi adiacenti con il peso associato al collegamento.
vector<list<pair<unsigned int, double>>> ListaAdiacenza(PolyhedronMesh& mesh);

// Applica l'algoritmo di Dijkstra per trovare il cammino minimo tra due nodi.
// Restituisce la sequenza di nodi che formano il cammino minimo e la lunghezza totale.
vector<unsigned int> DijkstraCamminoMinimo(const vector<list<pair<unsigned int, double>>>& LA,
                                           unsigned int id_start,
                                           unsigned int id_end,
                                           double& lunghezzaTotale,
                                           vector<int>& predecessori);

// A partire dal cammino minimo (sequenza di nodi), estrae i lati della mesh attraversati.
vector<unsigned int> latiCamminoMinimo(PolyhedronMesh& polyhedron,
                                       const vector<unsigned int>& cammino,
                                       const string& outputDirectory);

// Esegue la triangolazione geodetica di Classe I sulla mesh di input.
// Il parametro b determina la suddivisione, e la mesh triangolata è scritta in output.
bool TriangolaClasseI(int b,
                      const PolyhedronMesh& mesh_input,
                      PolyhedronMesh& mesh_output);

// Esegue la triangolazione geodetica di Classe II sulla mesh di input.
bool TriangolaClasseII(int b_input,
                      const PolyhedronMesh& mesh_input,
                      PolyhedronMesh& mesh_output);

// Esporta la mesh su file in formato testuale, salvando punti, lati e facce nella directory specificata.
bool EsportaMeshSuFile(const PolyhedronMesh& mesh, const std::string& outputDirectory);

// Cerca un vertice all'interno di una lista di vertici, tenendo conto di una tolleranza su coordinate.
// Se il vertice è già presente nella lista, ne restituisce l'indice. Altrimenti restituisce -1.
int TrovaVertice(const Eigen::Vector3d& v, const std::vector<Eigen::Vector3d>& lista, double toll = 1e-8);

// Genera un esagono definito da sei vertici dati.
// Restituisce nuove facce e aggiornamenti sugli ID delle celle coinvolte.
void GeneraEsagono(
    const std::vector<Vector3d>& /*nuovi_vertici*/,
    std::vector<std::vector<int>>& nuove_facce,
    int id_1, int id_2, int id_3, int id_4, int id_5, int id_6,
    int& id_b_low, int& it_es_tri);

// Ordina le facce attorno a un dato vertice.
// Utile ad esempio per costruire mesh duali.
vector<int> OrdinaFacceAttornoAlVertice(int vertex_id, const PolyhedronMesh& mesh);

// Costruisce la mesh duale a partire da una mesh di input.
bool CostruisciDualMesh(const PolyhedronMesh& StartPolyhedron, PolyhedronMesh& DualPolyhedron);

}
