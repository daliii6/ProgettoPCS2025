#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <list>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PolyhedronLibrary {

/*Struttura dati per rappresentare una mesh di poliedro.
 
 La mesh è definita da celle di dimensioni da 0D a 3D:
 - 0D: punti (nodi)
 - 1D: segmenti (spigoli)
 - 2D: facce (poligoni)
 - 3D: celle volumetriche
 */
struct PolyhedronMesh
{

    unsigned int NumCell0Ds; ///< Numero di celle 0D (vertici)
    unsigned int NumCell1Ds; ///< Numero di celle 1D (spigoli)
    unsigned int NumCell2Ds; ///< Numero di celle 2D (facce)
    unsigned int NumCell3Ds; ///< Numero di celle 3D (poliedri)

    // Identificativi univoci delle celle


    vector<unsigned int> Cell0DsId; ///< ID univoci delle celle 0D, dimensione: 1 x NumCell0Ds
    vector<unsigned int> Cell1DsId; ///< ID univoci delle celle 1D, dimensione: 1 x NumCell1Ds
    vector<unsigned int> Cell2DsId; ///< ID univoci delle celle 2D, dimensione: 1 x NumCell2Ds


    /*
    Matrice double (3 x NumCell0Ds): ogni colonna rappresenta le coordinate di un punto nello spazio 3D.
    */
    MatrixXd Cell0DsCoordinates;

    /*
    Matrice intera (2 x NumCell1Ds): ogni colonna contiene gli ID dei due nodi che definiscono lo spigolo (da -> a).
    */
    MatrixXi Cell1DsExtrema;

    /*
    Ogni elemento del vettore è un vettore contenente gli ID dei vertici che costituiscono la faccia.
    L'ordine è importante per la definizione dell'orientamento.
    */
    vector<vector<unsigned int>> Cell2DsVertices;

    /*
    Ogni elemento del vettore è un vettore contenente gli ID degli spigoli (celle 1D) che compongono la faccia.
    */
    vector<vector<unsigned int>> Cell2DsEdges;
};

}