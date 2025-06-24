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
#include <queue>
#include <utility>
#include <algorithm>

#include <limits>
#include "UCDUtilities.hpp"
using namespace Gedim;


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

        unsigned int id;
        converter >> id
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

        unsigned int id;
        converter >> id >> polyhedron.Cell1DsExtrema(0, id) >> polyhedron.Cell1DsExtrema(1, id);
        				
        // test per verificare che nessun edge abbia lunghezza zero
        if(polyhedron.Cell1DsExtrema(0, id)  == polyhedron.Cell1DsExtrema(1, id)){
            cerr << "at least one edge has zero length";
            return false;
        }
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

        unsigned int id, numVertices, numEdges;
        converter >> id >> numVertices;

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

/* Cerca un vettore v all’interno di una lista
Restituisce l’indice del primo vettore "vicino" (entro toll)
Se non trova nulla, restituisce -1
Usata per evitare duplicati nella triangolazione */
int TrovaVertice(const Vector3d& v, const std::vector<Vector3d>& lista, double toll) {
    for (int i = 0; i < (int)lista.size(); i++)  
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

//funzione che mi costruisce lista di adiacenza prende input mesh
vector<list<pair<unsigned int, double>>> ListaAdiacenza(PolyhedronMesh& mesh)
{
	unsigned int numVertices = mesh.NumCell0Ds;
    vector<list<pair<unsigned int, double>>> LA(numVertices); 
//itero spigoli	
    for (int i = 0; i < mesh.Cell1DsExtrema.cols(); i++) 
	{
	// estraggo id estremi arco
        unsigned int id_v1 = mesh.Cell1DsExtrema(0, i);
        unsigned int id_v2 = mesh.Cell1DsExtrema(1, i);

        Vector3d coord1 = mesh.Cell0DsCoordinates.col(id_v1);
        Vector3d coord2 = mesh.Cell0DsCoordinates.col(id_v2);
		//calcolo lunghezza arco
		Vector3d coord_diff=coord1-coord2;
		double dist = sqrt(coord_diff[0]*coord_diff[0]+coord_diff[1]*coord_diff[1]+coord_diff[2]*coord_diff[2]);
		pair<unsigned int, double> arco1(id_v2, dist);
		pair<unsigned int, double> arco2(id_v1, dist);
	//push della coppia con id vertice confinante e lunghezza arco
        LA[id_v1].push_back(arco1);
        LA[id_v2].push_back(arco2);  
    }

    return LA;
	
}
//Algoritmo djkstra per il calcolo cammino minimo tra id_start e id_end	
vector<unsigned int> DijkstraCamminoMinimo(const vector<list<pair<unsigned int, double>>>& LA,
                                           unsigned int id_start,
                                           unsigned int id_end,
                                           double& lunghezzaTotale,
                                           vector<int>& predecessori)
{
	// LA è lista di adiacenza: vettore di liste, ogni lista rappresenta un id ed è fatta di pair contenenti id vertice adiacente e lunghezza arco
    unsigned int numVertici = LA.size();
	//  distanza:=vettore distanze inizializzato a +inf
    vector<double> distanza(numVertici, numeric_limits<double>::infinity());
	//predecessori:= vettore per vedere i predecessori di ogni nodo del grafo che avrà predecessori di ogni nodo in base a cammino minimo
    predecessori.assign(numVertici, -1);
	//visitato:= vettore booleano
    vector<bool> visitato(numVertici, false);
	//coda con priorità che estrae elemento minima distanza, primo parametro elemento secondo contenitore terzo funzione comparazione
    priority_queue<pair<double, unsigned int>, vector<pair<double, unsigned int>>, greater<>> coda;

    distanza[id_start] = 0.0;
    coda.push({0.0, id_start});

    while (!coda.empty()) 
	{
        auto [d, u] = coda.top();
        coda.pop();
		//se u visitato skippo iterazione,altrimenti lo marko
        if (visitato[u]) continue;
        visitato[u] = true;
		//se raggiungo id_end ho finito
        if (u == id_end) break;
		//scorro adiacenza e aggiorno distanze per i vertici adiacenti a u
        for (const auto& [v, peso] : LA[u]) 
		{
		//guardo se passando da u ho distanza minore di v da id start
            if (distanza[v] > distanza[u] + peso) 
			{
                distanza[v] = distanza[u] + peso;
                predecessori[v] = u;
                coda.push({distanza[v], v});
            }
        }
    }
	//distanza minima da id end e id start
    lunghezzaTotale = distanza[id_end];

    vector<unsigned int> cammino;
	//non andato a buon fine
    if (predecessori[id_end] == -1) return cammino;
	//ricostruisco all'indietro il cammino minimo, dove v assumera i valori degli id di interesse dei predecessori
    for (int v = id_end; v != -1; v = predecessori[v])
        cammino.push_back(v);
	//inverto perchè voglio da id start a id end
    reverse(cammino.begin(), cammino.end());	
    return cammino;
}
// funzione che restituisce vettore dei lati  che costituiscono cammino minimo trovato

vector<unsigned int> latiCamminoMinimo(PolyhedronMesh& polyhedron,
                                       const vector<unsigned int>& cammino,
                                       const string& outputDirectory)
{
    Eigen::MatrixXi estremi = polyhedron.Cell1DsExtrema;
    unsigned int n = polyhedron.NumCell1Ds;
    unsigned int m = cammino.size();
    vector<unsigned int> lati;
	//ciclo coppie di vertici id cammino minimo
    for (unsigned int i = 0; i < m - 1; i++) {
        unsigned int v1 = cammino[i];
        unsigned int v2 = cammino[i + 1];
	    // cerco lato corrispondente alla coppia di indici , controllo bidirezionale perchè grafo non orientato, mi va bene ogni verso 
        for (unsigned int j = 0; j < n; j++) {
            if ((static_cast<unsigned int>(estremi(0, j)) == v1 && static_cast<unsigned int>(estremi(1, j)) == v2) ||
                (static_cast<unsigned int>(estremi(0, j)) == v2 && static_cast<unsigned int>(estremi(1, j)) == v1)) {
                lati.push_back(j);
                break;
            }
        }
    }

    // === ESPORTAZIONE UNIFICATA ===
    vector<double> PointData(polyhedron.NumCell0Ds, 0.0);
    for (auto id : cammino) PointData[id] = 1.0;

    vector<double> EdgeData(polyhedron.NumCell1Ds, 0.0);
    for (auto id : lati) EdgeData[id] = 1.0;

    Gedim::UCDProperty<double> PropPoints;
    PropPoints.Label = "ShortPath";
    PropPoints.UnitLabel = "";
    PropPoints.Size = PointData.size();
    PropPoints.NumComponents = 1;
    PropPoints.Data = PointData.data();

    Gedim::UCDProperty<double> PropEdges;
    PropEdges.Label = "ShortPath";
    PropEdges.UnitLabel = "";
    PropEdges.Size = EdgeData.size();
    PropEdges.NumComponents = 1;
    PropEdges.Data = EdgeData.data();

    Gedim::UCDUtilities utils;
    utils.ExportPoints(outputDirectory + "/Cell0Ds.inp", polyhedron.Cell0DsCoordinates, { PropPoints });
    utils.ExportSegments(outputDirectory + "/Cell1Ds.inp", polyhedron.Cell0DsCoordinates, polyhedron.Cell1DsExtrema, {}, { PropEdges });

    return lati;
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
bool TriangolaClasseI(int b,
                      const PolyhedronMesh& mesh_input,
                      PolyhedronMesh& mesh_output)
{
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

    // costruiamo gli spigoli dei triangoli
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


//funzione che da esagono mi prende facce triangolari
void GeneraEsagono(
    const std::vector<Vector3d>& /*nuovi_vertici*/,
    std::vector<std::vector<int>>& nuove_facce,
    int id_1, int id_2, int id_3, int id_4, int id_5, int id_6,
    int& id_b_low, int& it_es_tri)
{
    nuove_facce.push_back({id_1, id_2, id_b_low});
    nuove_facce.push_back({id_2, id_3, id_b_low});
    nuove_facce.push_back({id_3, id_4, id_b_low});
    nuove_facce.push_back({id_4, id_5, id_b_low});
    nuove_facce.push_back({id_5, id_6, id_b_low});
    nuove_facce.push_back({id_6, id_1, id_b_low});
    it_es_tri = it_es_tri + 6;
}




// triangolo classe 2 che mi restituisce un oggetto PolyhedronMesh 
bool TriangolaClasseII(int b_input,
                      const PolyhedronMesh& mesh_input,
                      PolyhedronMesh& mesh_output)
{
    int b = b_input;

    std::vector<Vector3d> nuovi_vertici;
    std::vector<std::vector<int>> nuove_facce;
	//itero su ogni faccia
    for (size_t f = 0; f < mesh_input.Cell2DsVertices.size(); f++) 
    {
		
        const auto& face = mesh_input.Cell2DsVertices[f];
        if (face.size() != 3) {
            std::cerr << "Errore: faccia non triangolare.\n";
            return false;
        }

        Vector3d v0 = mesh_input.Cell0DsCoordinates.col(face[0]);
        Vector3d v1 = mesh_input.Cell0DsCoordinates.col(face[1]);
        Vector3d v2 = mesh_input.Cell0DsCoordinates.col(face[2]);
	//costruzione griglia baricentrica uguale a cl 1
        std::vector<std::vector<int>> index_grid(b + 1);

        for (int i = 0; i <= b; i++) {
            index_grid[i].resize(b - i + 1);
            for (int j = 0; j <= b - i; j++) {
                int k = b - i - j;
                Vector3d punto = (double(i)/b)*v1 + (double(j)/b)*v2 + (double(k)/b)*v0;
                punto.normalize();
                int id = TrovaVertice(punto, nuovi_vertici);
                if (id == -1) {
                    id = nuovi_vertici.size();
                    nuovi_vertici.push_back(punto);
                }
                index_grid[i][j] = id;
            }
        }
	//caso b=c=1 semplice, prendo solo baricentro triangolo e tutti i punti medi dei lati
        if (b == 1) {
            int id_p1 = index_grid[0][0];
            int id_p2 = index_grid[0][1];
            int id_p3 = index_grid[1][0];

            Vector3d p1 = nuovi_vertici[id_p1];
            Vector3d p2 = nuovi_vertici[id_p2];
            Vector3d p3 = nuovi_vertici[id_p3];

            Vector3d b_up = ((p1 + p2 + p3) / 3.0).normalized();
			
            int id_bar_up = TrovaVertice(b_up, nuovi_vertici);
            if (id_bar_up == -1) {
                id_bar_up = nuovi_vertici.size();
                nuovi_vertici.push_back(b_up);
            }

            Vector3d m12_up = ((p1 + p2) / 2.0).normalized();
            int id_m12 = TrovaVertice(m12_up, nuovi_vertici);
            if (id_m12 == -1) {
                id_m12 = nuovi_vertici.size();
                nuovi_vertici.push_back(m12_up);
            }

            Vector3d m23_up = ((p2 + p3) / 2.0).normalized();
            int id_m23 = TrovaVertice(m23_up, nuovi_vertici);
            if (id_m23 == -1) {
                id_m23 = nuovi_vertici.size();
                nuovi_vertici.push_back(m23_up);
            }

            Vector3d m31_up = ((p3 + p1) / 2.0).normalized();
            int id_m31 = TrovaVertice(m31_up, nuovi_vertici);
            if (id_m31 == -1) {
                id_m31 = nuovi_vertici.size();
                nuovi_vertici.push_back(m31_up);
            }

            nuove_facce.push_back({id_p1, id_m12, id_bar_up});
            nuove_facce.push_back({id_m12, id_p2, id_bar_up});
            nuove_facce.push_back({id_p2, id_m23, id_bar_up});
            nuove_facce.push_back({id_m23, id_p3, id_bar_up});
            nuove_facce.push_back({id_p3, id_m31, id_bar_up});
            nuove_facce.push_back({id_m31, id_p1, id_bar_up});
        } 
        else 
        {
			int it=0;
			int it_facce_lat=0;
			int it_esa=0;
			int it_es_tri=0;
		// doppio for in cui seleziono due triangoli consecutivi orizzontalmente e quello sopra e scorro lungo righe e dopo salto in alto verso cima
            for (int a = 0; a < b; a++) {
                int k = index_grid[a].size();
				
                for (int m = 0; m < k - 1; m++) {
                    /* if (a + 1 >= index_grid.size() || 
                        m >= index_grid[a+1].size() ||
                        m+1 >= index_grid[a].size() ||
                        m+1 >= index_grid[a+1].size()) {
                        continue; // evita accessi fuori limite
                    } */
			//vertici triangolo  normale
                    int id_p1 = index_grid[a][m];
                    int id_p2 = index_grid[a][m+1];
                    int id_p3 = index_grid[a+1][m];
					it = it+3;
                    Vector3d p1 = nuovi_vertici[id_p1];
                    Vector3d p2 = nuovi_vertici[id_p2];
                    Vector3d p3 = nuovi_vertici[id_p3];
			// prendo primo baricentro
                    Vector3d b_up = ((p1 + p2 + p3) / 3.0).normalized();
					it+=1;
                    int id_bar_up = TrovaVertice(b_up, nuovi_vertici);
                    if (id_bar_up == -1) {
                        id_bar_up = nuovi_vertici.size();
                        nuovi_vertici.push_back(b_up);
                    }
			//se a =0 prendo punti alla base della faccia di partenza
                    if (a == 0) {
                        Vector3d m12_up = ((p1 + p2) / 2.0).normalized();
                        int id_m12 = TrovaVertice(m12_up, nuovi_vertici);
                        if (id_m12 == -1) {
                            id_m12 = nuovi_vertici.size();
                            nuovi_vertici.push_back(m12_up);
                        }
						it+=1;
						nuove_facce.push_back({id_p1, id_bar_up, id_m12});
                        nuove_facce.push_back({id_p2, id_bar_up, id_m12});
						it_facce_lat +=1;
						it_facce_lat +=1;
                    }
		//in questi due if se m è all'inizio o alla fine della griglia prendo punt medi laterali
                    if (m == k - 2) {
                        Vector3d m23_up = ((p2 + p3) / 2.0).normalized();
						it+=1;
                        int id_m23 = TrovaVertice(m23_up, nuovi_vertici);
                        if (id_m23 == -1) {
                            id_m23 = nuovi_vertici.size();
                            nuovi_vertici.push_back(m23_up);
                        }

                        nuove_facce.push_back({id_p2, id_bar_up, id_m23});
                        nuove_facce.push_back({id_p3, id_bar_up, id_m23});
						it_facce_lat +=1;
						it_facce_lat +=1;
                    }

                    if (m == 0) {
                        Vector3d m31_up = ((p3 + p1) / 2.0).normalized();
                        int id_m31 = TrovaVertice(m31_up, nuovi_vertici);
                        if (id_m31 == -1) {
                            id_m31 = nuovi_vertici.size();
                            nuovi_vertici.push_back(m31_up);
                        }
						it+=1;

                        nuove_facce.push_back({id_p1, id_bar_up, id_m31});
                        nuove_facce.push_back({id_p3, id_bar_up, id_m31});
						it_facce_lat +=1;
						it_facce_lat +=1;
                    }

                    // triangolo inferiore
                    /* if (a < b - 1 && m+1 < index_grid[a].size() &&
                        m < index_grid[a+1].size() && 
                        m+1 < index_grid[a+1].size()) {

                        int id_d1 = index_grid[a][m+1];
                        int id_d2 = index_grid[a+1][m];
                        int id_d3 = index_grid[a+1][m+1];
						
                        Vector3d d1 = nuovi_vertici[id_d1];
                        Vector3d d2 = nuovi_vertici[id_d2];
                        Vector3d d3 = nuovi_vertici[id_d3]; */

			//uso static_cast che mi connverte da unsigned int a intero perchè m+1 lo è
                    if (a < b - 1 && m + 1 < static_cast<int>(index_grid[a + 1].size()))
					{
						it +=1;
						int id_d3 = index_grid[a+1][m+1];
						Vector3d d3 = nuovi_vertici[id_d3];
						Vector3d b_low = ((p2 +p3+d3) / 3.0).normalized();
						
						int id_bar_low = TrovaVertice(b_low, nuovi_vertici);
						if (id_bar_low == -1) 
						{
							id_bar_low = nuovi_vertici.size();
							nuovi_vertici.push_back(b_low);
						}
                    }
					if (a + 2 <= b &&m + 2 < (int)index_grid[a].size() &&m + 1 < (int)index_grid[a+1].size() &&m < (int)index_grid[a+2].size()) 
					{
						int id_p2_next = index_grid[a][m+2];
						int id_p3_next = index_grid[a+1][m+1];
						int id_p3_next_next = index_grid[a+2][m];

						Vector3d p2_next = nuovi_vertici[id_p2_next];
						Vector3d p3_next = nuovi_vertici[id_p3_next];
						Vector3d p3_next_next = nuovi_vertici[id_p3_next_next];

						Vector3d b_up_next = ((p2 + p2_next + p3_next) / 3.0).normalized();
						int id_bar_up_next = TrovaVertice(b_up_next, nuovi_vertici);
						if (id_bar_up_next == -1) {
							id_bar_up_next = nuovi_vertici.size();
							nuovi_vertici.push_back(b_up_next);
						}

						Vector3d b_up_next_next = ((p3 + p3_next + p3_next_next) / 3.0).normalized();
						int id_bar_up_next_next = TrovaVertice(b_up_next_next, nuovi_vertici);
						if (id_bar_up_next_next == -1) 
						{
							id_bar_up_next_next = nuovi_vertici.size();
							nuovi_vertici.push_back(b_up_next_next);
						}
						it_esa+=1;
						//// partizione
						int id_1=id_bar_up;
						int id_2=id_p2;
						int id_3=id_bar_up_next;
						int id_4=id_p3_next;
						int id_5=id_bar_up_next_next;
						int id_6=id_p3;
						
						//
                        // int id_bar_esagono;
						Vector3d b_low = ((p2 +p3+p3_next) / 3.0).normalized();
						
						int id_b_low = TrovaVertice(b_low, nuovi_vertici);
						if (id_b_low == -1) 
						{
							id_b_low = nuovi_vertici.size();
							nuovi_vertici.push_back(b_low);
						}
						 GeneraEsagono(
										nuovi_vertici, nuove_facce,
										id_1,id_2, id_3,
										id_4, id_5, id_6,id_b_low,it_es_tri); 
						
					}
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


/* Prende in input l'identificativo di un vertice `vertex_id`
e una struttura `PolyhedronMesh` che rappresenta un poliedro. Restituisce
un vettore di interi contenente gli ID delle facce adiacenti a quel vertice,
ordinati in modo tale che ogni faccia condivida un lato (due vertici) con la successiva.

Questo ordinamento è utile, per costruire il poliedro duale,
dove ogni vertice del duale si trova nel baricentro di una faccia e i vertici adiacenti
devono essere collegati in ordine corretto per formare una faccia planare.*/


std::vector<int> OrdinaFacceAttornoAlVertice(int vertex_id, const PolyhedronMesh& mesh) {
    std::vector<int> facce_condivise;

    // Trova tutte le facce che contengono il vertice `vertex_id`
    for (size_t fid = 0; fid < mesh.NumCell2Ds; fid++) {
        const auto& face = mesh.Cell2DsVertices[fid];
        if (std::find(face.begin(), face.end(), vertex_id) != face.end())
            facce_condivise.push_back(fid);
    }

    // Se ci sono 2 o meno facce, l'ordinamento non è necessario
    if (facce_condivise.size() <= 2)
        return facce_condivise;

    // Inizializza il vettore che conterrà le facce ordinate
    std::vector<int> ordered_faces;

    // Insieme per tenere traccia delle facce già usate (evita duplicati)
    std::set<int> usate;

    // Partiamo arbitrariamente dalla prima faccia trovata
    ordered_faces.push_back(facce_condivise[0]);
    usate.insert(facce_condivise[0]);

    // Flag per indicare se è stata trovata una nuova faccia adiacente
    bool trovato = true;

    // Continua finché non abbiamo ordinato tutte le facce o finché troviamo adiacenze
    while (ordered_faces.size() < facce_condivise.size() && trovato) {
        trovato = false;

        // Ultima faccia inserita nell'ordine
        int ultima = ordered_faces.back();
        const auto& f1 = mesh.Cell2DsVertices[ultima];

        // Cerchiamo una nuova faccia che sia adiacente all'ultima
        for (int fid : facce_condivise) {
            if (usate.count(fid)) continue;  // Salta se già usata

            const auto& f2 = mesh.Cell2DsVertices[fid];

            // Conta quanti vertici f1 e f2 hanno in comune
            int count_comuni = 0;
            for (int v : f2)
                if (std::find(f1.begin(), f1.end(), v) != f1.end())
                    count_comuni++;

            // Condizione di adiacenza: f2 contiene vertex_id e condivide 2 vertici con f1
            if (std::find(f2.begin(), f2.end(), vertex_id) != f2.end() && count_comuni == 2) {
                // Aggiunge la nuova faccia all'ordine
                ordered_faces.push_back(fid);
                usate.insert(fid);
                trovato = true;
                break;  // Passa alla prossima iterazione
            }
        }
    }

    // Restituisce la lista di facce ordinate attorno al vertice
    return ordered_faces;
}



/* Dato un poliedro triangolato `StartPolyhedron`, questa funzione costruisce il suo duale
`DualPolyhedron`, seguendo la regola della dualità:
- Ogni faccia del poliedro originale diventa un vertice del duale (posto nel baricentro).
- Ogni vertice del poliedro originale genera una faccia del duale, composta dai vertici
     duali delle facce che gli sono adiacenti (ordinati ciclicamente attorno al vertice).
- Gli spigoli del duale collegano i vertici duali di facce adiacenti.

 La funzione aggiorna le strutture Cell0Ds, Cell1Ds, Cell2Ds di `DualPolyhedron`.

 Ritorna true se la costruzione è andata a buon fine.*/

bool CostruisciDualMesh(const PolyhedronMesh& StartPolyhedron, PolyhedronMesh& DualPolyhedron)
{
    // Vettore dei nuovi vertici (baricentri delle facce originali)
    std::vector<Vector3d> nuovi_vertici;

    // Facce del duale: ciascuna è una lista di indici di baricentri (cioè vertici duali)
    std::vector<std::vector<int>> nuove_facce;

    // Mappa da ID faccia originale → ID vertice nel duale (baricentro)
    std::map<int, int> face_to_dual_vertex;

    // FASE 1: Calcola i baricentri delle facce originali → vertici del duale
    for (size_t fid = 0; fid < StartPolyhedron.NumCell2Ds; ++fid) {
        const auto& face = StartPolyhedron.Cell2DsVertices[fid];

        // Calcolo del baricentro come media delle coordinate dei vertici
        Vector3d bar = Vector3d::Zero();
        for (int idx : face)
            bar += StartPolyhedron.Cell0DsCoordinates.col(idx);
        bar /= face.size();

        // Proiezione del baricentro sulla sfera unitaria
        bar.normalize();

        // Controlla se il vertice è già presente (per evitare duplicati)
        int id = TrovaVertice(bar, nuovi_vertici);
        if (id == -1) {
            id = nuovi_vertici.size();
            nuovi_vertici.push_back(bar);
        }

        // Registra la corrispondenza faccia originale → vertice del duale
        face_to_dual_vertex[fid] = id;
    }

    // FASE 2: Ogni vertice del poliedro originale genera una faccia del duale
    for (size_t vid = 0; vid < StartPolyhedron.NumCell0Ds; ++vid) {
        // Ordina le facce attorno al vertice vid (in senso ciclico)
        std::vector<int> ordered_faces = OrdinaFacceAttornoAlVertice(vid, StartPolyhedron);

        // Costruisci la faccia del duale con i vertici corrispondenti ai baricentri ordinati
        std::vector<int> faccia_duale;
        for (int fid : ordered_faces)
            faccia_duale.push_back(face_to_dual_vertex[fid]);

        nuove_facce.push_back(faccia_duale);
    }

    // FASE 3: Costruzione dei vertici del duale
    DualPolyhedron.NumCell0Ds = nuovi_vertici.size();
    DualPolyhedron.Cell0DsId.resize(nuovi_vertici.size());
    DualPolyhedron.Cell0DsCoordinates = MatrixXd(3, nuovi_vertici.size());

    for (int i = 0; i < (int)nuovi_vertici.size(); ++i) {
        DualPolyhedron.Cell0DsId[i] = i;
        DualPolyhedron.Cell0DsCoordinates.col(i) = nuovi_vertici[i];
    }

    // FASE 4: Costruzione degli spigoli del duale
    // Uso un set per evitare duplicati, tenendo conto della direzione
    std::set<std::pair<int, int>> spigoli_set;

    for (const auto& faccia : nuove_facce) {
        for (size_t i = 0; i < faccia.size(); ++i) {
            int u = faccia[i];
            int v = faccia[(i + 1) % faccia.size()];
            // Salva l'arco (u,v) solo se (v,u) non è già stato aggiunto
            if (spigoli_set.count({v, u}) == 0)
                spigoli_set.insert({u, v});
        }
    }

    DualPolyhedron.NumCell1Ds = spigoli_set.size();
    DualPolyhedron.Cell1DsId.resize(spigoli_set.size());
    DualPolyhedron.Cell1DsExtrema = MatrixXi(2, spigoli_set.size());

    // Mappa da coppia di vertici → ID spigolo
    std::map<std::pair<int, int>, int> edge_map;
    int eid = 0;
    for (const auto& e : spigoli_set) {
        DualPolyhedron.Cell1DsId[eid] = eid;
        DualPolyhedron.Cell1DsExtrema(0, eid) = e.first;
        DualPolyhedron.Cell1DsExtrema(1, eid) = e.second;
        edge_map[e] = eid;
        eid++;
    }

    // FASE 5: Costruzione delle facce del duale
    DualPolyhedron.NumCell2Ds = nuove_facce.size();
    DualPolyhedron.Cell2DsId.resize(nuove_facce.size());
    DualPolyhedron.Cell2DsVertices.resize(nuove_facce.size());
    DualPolyhedron.Cell2DsEdges.resize(nuove_facce.size());

    for (int i = 0; i < (int)nuove_facce.size(); ++i) {
        DualPolyhedron.Cell2DsId[i] = i;

        // Vertici della faccia (convertiti a unsigned int)
        std::vector<unsigned int> vertici_unsigned(nuove_facce[i].begin(), nuove_facce[i].end());
        DualPolyhedron.Cell2DsVertices[i] = vertici_unsigned;

        // Trova gli spigoli corrispondenti
        std::vector<unsigned int> edges;
        for (size_t j = 0; j < nuove_facce[i].size(); ++j) {
            int u = nuove_facce[i][j];
            int v = nuove_facce[i][(j + 1) % nuove_facce[i].size()];
            if (edge_map.count({u, v})) {
                edges.push_back(edge_map[{u, v}]);
            } else if (edge_map.count({v, u})) {
                edges.push_back(edge_map[{v, u}]);
            } else {
                std::cerr << "Errore: spigolo non trovato tra " << u << " e " << v << std::endl;
            }
        }

        DualPolyhedron.Cell2DsEdges[i] = edges;
    }

    // Costruzione completata con successo
    return true;
}



bool EsportaMeshSuFile(const PolyhedronMesh& mesh, const string& outputDirectory)
{
    namespace fs = std::filesystem;
    fs::create_directories(outputDirectory);

    // === Cell0Ds.txt ===
    ofstream file0(outputDirectory + "/Cell0Ds.txt");
    file0 << "Id;X;Y;Z\n";
    for (size_t i = 0; i < mesh.NumCell0Ds; ++i)
        file0 << i << ";" << mesh.Cell0DsCoordinates(0, i) << ";"
              << mesh.Cell0DsCoordinates(1, i) << ";" << mesh.Cell0DsCoordinates(2, i) << "\n";
    file0.close();

    // === Cell1Ds.txt ===
    ofstream file1(outputDirectory + "/Cell1Ds.txt");
    file1 << "Id;Origine;Fine\n";
    for (size_t i = 0; i < mesh.NumCell1Ds; ++i)
        file1 << i << ";" << mesh.Cell1DsExtrema(0, i) << ";" << mesh.Cell1DsExtrema(1, i) << "\n";
    file1.close();

    // === Cell2Ds.txt ===
    ofstream file2(outputDirectory + "/Cell2Ds.txt");
    file2 << "Id;NumVertici;Vertici...;NumSpigoli;Spigoli...\n";
    for (size_t i = 0; i < mesh.NumCell2Ds; ++i) {
        file2 << i << ";" << mesh.Cell2DsVertices[i].size();
        for (int v : mesh.Cell2DsVertices[i]) file2 << ";" << v;
        file2 << ";" << mesh.Cell2DsEdges[i].size();
        for (int e : mesh.Cell2DsEdges[i]) file2 << ";" << e;
        file2 << "\n";
    }
    file2.close();

    // === Cell3Ds.txt ===
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
