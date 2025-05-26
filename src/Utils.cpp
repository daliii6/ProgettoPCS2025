#include "Utils.hpp"
#include "Eigen/Eigen"
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <set>
#include <vector>
#include <map>

using namespace std;
using namespace Eigen;

namespace PolyhedronLibrary {

// Importa Cell0Ds da CSV
bool ImportCell0Ds(PolyhedronMesh& polyhedron, const string& InputFileDirectory)
{
    ifstream file(InputFileDirectory + "/Cell0Ds.csv");
    if(file.fail()) return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);
    file.close();
    listLines.pop_front();  // Rimuove l'header

    polyhedron.NumCell0Ds = listLines.size();
    if (polyhedron.NumCell0Ds == 0) return false;

    polyhedron.Cell0DsId.reserve(polyhedron.NumCell0Ds);
    polyhedron.Cell0DsCoordinates = MatrixXd::Zero(3, polyhedron.NumCell0Ds);

    for (string& line : listLines)
    {
        replace(line.begin(), line.end(), ';', ' ');
        istringstream converter(line);

        unsigned int id, marker;
        converter >> id >> marker
                  >> polyhedron.Cell0DsCoordinates(0, id)
                  >> polyhedron.Cell0DsCoordinates(1, id)
                  >> polyhedron.Cell0DsCoordinates(2, id);
        polyhedron.Cell0DsId.push_back(id);
    }

    return true;
}

// Importa Cell1Ds da CSV
bool ImportCell1Ds(PolyhedronMesh& polyhedron, const string& InputFileDirectory)
{
    ifstream file(InputFileDirectory + "/Cell1Ds.csv");
    if(file.fail()) return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);
    file.close();
    listLines.pop_front();

    polyhedron.NumCell1Ds = listLines.size();
    if (polyhedron.NumCell1Ds == 0) return false;

    polyhedron.Cell1DsId.reserve(polyhedron.NumCell1Ds);
    polyhedron.Cell1DsExtrema = MatrixXi(2, polyhedron.NumCell1Ds);

    for (string& line : listLines)
    {
        replace(line.begin(), line.end(), ';', ' ');
        istringstream converter(line);

        unsigned int id, marker;
        converter >> id >> marker >> polyhedron.Cell1DsExtrema(0, id) >> polyhedron.Cell1DsExtrema(1, id);
        polyhedron.Cell1DsId.push_back(id);
    }

    return true;
}

// Importa Cell2Ds da CSV
bool ImportCell2Ds(PolyhedronMesh& polyhedron, const string& InputFileDirectory)
{
    ifstream file(InputFileDirectory + "/Cell2Ds.csv");
    if (file.fail()) return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);
    file.close();
    listLines.pop_front();

    polyhedron.NumCell2Ds = listLines.size();
    if (polyhedron.NumCell2Ds == 0) return false;

    polyhedron.Cell2DsId.reserve(polyhedron.NumCell2Ds);
    polyhedron.Cell2DsVertices.reserve(polyhedron.NumCell2Ds);
    polyhedron.Cell2DsEdges.reserve(polyhedron.NumCell2Ds);

    for (string& line : listLines)
    {
        replace(line.begin(), line.end(), ';', ' ');
        istringstream converter(line);

        unsigned int id, marker, numVertices, numEdges;
        converter >> id >> marker >> numVertices;

        vector<unsigned int> vertices(numVertices);
        for (unsigned int i = 0; i < numVertices; i++) converter >> vertices[i];

        converter >> numEdges;
        vector<unsigned int> edges(numEdges);
        for (unsigned int i = 0; i < numEdges; i++) converter >> edges[i];

        polyhedron.Cell2DsId.push_back(id);
        polyhedron.Cell2DsVertices.push_back(vertices);
        polyhedron.Cell2DsEdges.push_back(edges);
    }

    return true;
}

/* Cerca un vettore v all’interno di lista
Restituisce l’indice del primo vettore "vicino" (entro toll)
Se non trova nulla, restituisce -1
Usata per evitare duplicati nella triangolazione */
int TrovaVertice(const Vector3d& v, const std::vector<Vector3d>& lista, double toll) {
    for (int i = 0; i < (int)lista.size(); i++)  // int lista.size si puo fare???
        if ((lista[i] - v).norm() < toll)
            return i; //se è stato trovato un vertice uguale
    return -1;  // se NON è stato trovato un vertice uguale
}


/* Importa l'intero poliedro:
Chiama le tre funzioni di importazione sopra
Ritorna true solo se tutte le letture hanno successo */
bool ImportPolyhedronMesh(PolyhedronMesh& polyhedron, const string& InputFileDirectory)
{
    return ImportCell0Ds(polyhedron, InputFileDirectory)
        && ImportCell1Ds(polyhedron, InputFileDirectory)
        && ImportCell2Ds(polyhedron, InputFileDirectory);
}

/* Triangolazione geodetica di classe I
Itera sulle facce triangolari del solido platonico
Costruisce una griglia baricentrica usando parametri (b,0)
Normalizza ogni punto sulla sfera (perché è un solido geodetico)
Evita duplicati usando TrovaVertice
Costruisce facce triangolari locali (nuove_facce) b^2(CREDO)
Costruisce tutti gli spigoli univoci
Assegna ID a vertici, spigoli e facce
Riempie le strutture di mesh_output*/
bool TriangolaClasseI(int p, int q, int b,
                      const PolyhedronMesh& mesh_input,
                      PolyhedronMesh& mesh_output, int c)
{
	/* //Enrico: ho fatto i due casi dove uno tra b,c in input è zero e ho definito b,c compatibilmente
	if (b_input==0 && c_input>=1)
	{
		int b=c_input;
		int c=b_input;
	}
	if(c_input==0 && b_input>=1)
	{
		int b=b_input;
		int c=c_input;
	}
	 */
    //inizializzo i nuovi contenitori usando le stesse strutture di polyhedron.hpp
    vector<Vector3d> nuovi_vertici;
    vector<vector<int>> nuove_facce;

    //ciclo sulle facce del solido platonico
    for (size_t f = 0; f < mesh_input.Cell2DsVertices.size(); f++) {  //perche usiamo size_t? perchè fare .size() restituisce un size_t
        const auto& face = mesh_input.Cell2DsVertices[f]; //legge i suoi 3 vertici
        if (face.size() != 3) {
            cerr << "Errore: faccia non triangolare.\n";
            return false;
        }

        // salvo i 3 vertici della faccia
        Vector3d v0 = mesh_input.Cell0DsCoordinates.col(face[0]); //prende le 3 coordinate del vertice face[0]
        Vector3d v1 = mesh_input.Cell0DsCoordinates.col(face[1]);
        Vector3d v2 = mesh_input.Cell0DsCoordinates.col(face[2]);

        /* griglia baricentrica interna:
        index_grid[i][j] conterrà l’ID del vertice 
        generato nella posizione {i,j} nella griglia.
        Nella triangolazione baricentrica, un punto 
        interno a un triangolo è ottenuto come combinazione 
        convessa dei vertici p = l0*v0+l1*v1+l2*v2, con l0+l1+l2 = 1
        con l0=k/b l1=i/b l2 = j/b*/

        vector<vector<int>> index_grid(b + 1); 
		for (int i = 0; i <= b; i++) {  //i+j+k=b, i, j, k >=0 (genera tutte le triple possibili)
            index_grid[i].resize(b - i + 1);
            for (int j = 0; j <= b - i; j++) {
                int k = b - i - j;
                Vector3d punto = (double(i)/b)*v1 + (double(j)/b)*v2 + (double(k)/b)*v0; //punto in coordinate cartesiane
                punto.normalize();  //proiezione sulla sfera unitaria centrata nell'origine

                //cerca se il punto esiste già, altrimenti aggiunge un id
                int id = TrovaVertice(punto, nuovi_vertici);
                /*Se trova un vertice già “abbastanza uguale”, 
                restituisce il suo indice
                Se non trova nulla, restituisce -1*/
                if (id == -1) {
                    id = nuovi_vertici.size();
                    nuovi_vertici.push_back(punto); 
                }
                index_grid[i][j] = id;
            }
        }

        //cicla sulle celle della griglia
        for (int i = 0; i < b; i++) {
            for (int j = 0; j < b - i; j++) {
                int a = index_grid[i][j]; //vertice in basso a sx
                int b1 = index_grid[i + 1][j]; //vertice in basso a dx
                int c = index_grid[i][j + 1]; //vertice in alto a sx
                nuove_facce.push_back({a, b1, c}); //triangolo a sx della cella
                // se sei dentro la griglia aggiungi il secondo triangolo che completa la cella
                if (i + j < b - 1) {
                    int d = index_grid[i + 1][j + 1]; //vertice in alto a dx
                    nuove_facce.push_back({b1, d, c});
                }
            }
        }
    }

    // popolo il Cell0Ds
    mesh_output.NumCell0Ds = nuovi_vertici.size();
    mesh_output.Cell0DsCoordinates = MatrixXd(3, nuovi_vertici.size());
    mesh_output.Cell0DsId.resize(nuovi_vertici.size());
    for (int i = 0; i < (int)nuovi_vertici.size(); i++) {
        mesh_output.Cell0DsId[i] = i;
        mesh_output.Cell0DsCoordinates.col(i) = nuovi_vertici[i];
    }

    /*set che contine tutte le coppie di vertici che 
    formano gli spigoli.
    Set garantisce che ogni spigolo sia presente una sola volta.
    Si aggiungeranno solo spigoli nuovi,
    ignorando i duplicati (anche se invertiti)*/
   set<pair<int, int>> spigoli_set;  

    // costruiamo gli spigoli dai triangoli
    for (const auto& faccia : nuove_facce) { //scorro tutte le facce
        for (int i = 0; i < 3; ++i) { //scorro tutti i vertici
            int u = faccia[i]; //vertice di partenza del lato
            int v = faccia[(i + 1) % 3]; //vertice di arrivo del lato successivo

            // se esiste già l'inverso (v → u), non inserire (u → v)
            if (spigoli_set.count({v, u}) == 0) {
                spigoli_set.insert({u, v});  // salva lo spigolo come appare nella faccia
            }
            /*questo permette di rispettare la condizione
            faces.vertices[e] == faces.edges[e].origin*/
        }
    }


    mesh_output.NumCell1Ds = spigoli_set.size();
    mesh_output.Cell1DsId.resize(spigoli_set.size());
    mesh_output.Cell1DsExtrema = MatrixXi(2, spigoli_set.size());

    int idx = 0; // assegna id allo spigolo
    map<pair<int, int>, int> spigolo_id_map; //spigolo_id_map mappa ogni coppia {u, v} all’ID assegnato
    for (const auto& e : spigoli_set) {
        mesh_output.Cell1DsId[idx] = idx; //assegna idx allo spigolo corrente
        mesh_output.Cell1DsExtrema(0, idx) = e.first;
        mesh_output.Cell1DsExtrema(1, idx) = e.second;
        spigolo_id_map[e] = idx; //Mappa {u, v} → idx per poter poi associare rapidamente le facce ai loro lati
        idx++;
    }

    mesh_output.NumCell2Ds = nuove_facce.size();
    mesh_output.Cell2DsId.resize(nuove_facce.size());
    mesh_output.Cell2DsVertices.resize(nuove_facce.size());
    mesh_output.Cell2DsEdges.resize(nuove_facce.size());

    for (int i = 0; i < (int)nuove_facce.size(); i++) {//Scorre su tutte le facce i generate dalla triangolazione baricentrica
        mesh_output.Cell2DsId[i] = i; //assegna id = i alla faccia
        //Copia la tripla {a, b, c} della faccia corrente in Cell2DsVertices[i]
        mesh_output.Cell2DsVertices[i] = std::vector<unsigned int>(nuove_facce[i].begin(), nuove_facce[i].end()); 

        vector<unsigned int> edges;
        //estraggo i 3 lati della faccia
        for (int j = 0; j < 3; j++) {
            int u = nuove_facce[i][j];
            int v = nuove_facce[i][(j + 1) % 3];

            // cerca lo spigolo nel verso giusto, altrimenti nel verso opposto
            /* Prima prova {u, v} → verso coerente con la faccia
            Se non lo trova, prova {v, u} → verso opposto
            In entrambi i casi, recupera l’ID numerico e lo aggiunge alla 
            lista degli spigoli della faccia*/
            if (spigolo_id_map.count({u, v})) {
                edges.push_back(spigolo_id_map[{u, v}]);
            } else if (spigolo_id_map.count({v, u})) {
                edges.push_back(spigolo_id_map[{v, u}]);
            } else {
                cerr << "Errore: spigolo non trovato tra " << u << " e " << v << endl;
                return false;
            }
        }
        mesh_output.Cell2DsEdges[i] = edges; //salva la lista degli spigoli per la faccia
    }
    return true;
}
//funzione che da esagono mi prede facce triangolari
/*void GeneraEsagono(
    const std::vector<Vector3d>& nuovi_vertici,
    std::vector<std::vector<int>>& nuove_facce,
    int id_1, int id_2, int id_3, int id_4, int id_5, int id_6,
    int& id_baricentro,
    std::vector<Vector3d>& vertici_out)
 {
    // Calcolo del baricentro
    Vector3d bar = (
        nuovi_vertici[id_1] + nuovi_vertici[id_2] +
        nuovi_vertici[id_3] + nuovi_vertici[id_4] +
        nuovi_vertici[id_5] + nuovi_vertici[id_6]
    ) / 6.0;
    bar.normalize();

    // Verifica se esiste già, altrimenti aggiungi
    id_baricentro = TrovaVertice(bar, vertici_out);
    if (id_baricentro == -1) {
        id_baricentro = vertici_out.size();
        vertici_out.push_back(bar);
    }

    // Crea le 6 facce triangolari dell’esagono
    nuove_facce.push_back({id_1, id_2, id_baricentro});
    nuove_facce.push_back({id_2, id_3, id_baricentro});
    nuove_facce.push_back({id_3, id_4, id_baricentro});
    nuove_facce.push_back({id_4, id_5, id_baricentro});
    nuove_facce.push_back({id_5, id_6, id_baricentro});
    nuove_facce.push_back({id_6, id_1, id_baricentro});
}*/





// Esporta i file .txt
bool EsportaMeshSuFile(const PolyhedronMesh& mesh, const string& outputDirectory)
{
    namespace fs = std::filesystem;
    fs::create_directories(outputDirectory);

    ofstream file0(outputDirectory + "/Cell0Ds.txt");
    file0 << "Id;Marker;X;Y;Z\n";
    for (size_t i = 0; i < mesh.NumCell0Ds; ++i)
        file0 << i << ";0;" << mesh.Cell0DsCoordinates(0, i) << ";"
             << mesh.Cell0DsCoordinates(1, i) << ";" << mesh.Cell0DsCoordinates(2, i) << "\n";
    file0.close();

    ofstream file1(outputDirectory + "/Cell1Ds.txt");
    file1 << "Id;Marker;Origine;Fine\n";
    for (size_t i = 0; i < mesh.NumCell1Ds; ++i)
        file1 << i << ";0;" << mesh.Cell1DsExtrema(0, i) << ";" << mesh.Cell1DsExtrema(1, i) << "\n";
    file1.close();

    ofstream file2(outputDirectory + "/Cell2Ds.txt");
    file2 << "Id;Marker;NumVertici;Vertici...;NumSpigoli;Spigoli...\n";
    for (size_t i = 0; i < mesh.NumCell2Ds; ++i) {
        file2 << i << ";0;" << mesh.Cell2DsVertices[i].size();
        for (int v : mesh.Cell2DsVertices[i]) file2 << ";" << v;
        file2 << ";" << mesh.Cell2DsEdges[i].size();
        for (int e : mesh.Cell2DsEdges[i]) file2 << ";" << e;
        file2 << "\n";
    }
    file2.close();

    ofstream file3(outputDirectory + "/Cell3Ds.txt");
    file3 << "Id;NumVertici;NumSpigoli;NumFacce;Vertici...;Spigoli...;Facce...\n";
    file3 << "0;" << mesh.NumCell0Ds << ";" << mesh.NumCell1Ds << ";" << mesh.NumCell2Ds;
    for (size_t i = 0; i < mesh.NumCell0Ds; ++i) file3 << ";" << i;
    for (size_t i = 0; i < mesh.NumCell1Ds; ++i) file3 << ";" << i;
    for (size_t i = 0; i < mesh.NumCell2Ds; ++i) file3 << ";" << i;
    file3 << "\n";
    file3.close();

    return true;
}

} // namespace PolyhedronLibrary
