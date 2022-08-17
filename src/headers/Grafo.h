/*
    CABEÇALHO
*/

#ifndef GRAFO_H
#define GRAFO_H 1

#include "Vertice.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <forward_list>
#include <list>

using namespace std;

class Grafo {

    //Atributes
private:
    int order;
    int number_edges;
    bool directed;
    bool weighted_edge;
    bool weighted_node;
    Vertice *first_node;
    Vertice *last_node;
    int cluster;
    string clusterType;
    vector<tuple<int, int>> clustersLimits;
    double clustersCapacity;
    float inferiorLimit;
    float upperLimit;
    float currentLimit;
    float maxBenefit;
    int fileType;

public:
    //Constructor
    Grafo(int order, bool directed, bool weighted_edge, bool weighted_node);
    Grafo(int order, int cluster, string clusterType, vector<tuple<int, int>> clustersLimits);
    Grafo(int order, int cluster, double clustersCapacity);
    Grafo(int inferiorLimit, int upperLimit);
    Grafo(int argc, const char **argv);

    //Destructor
    ~Grafo();

    //Getters
    int getOrder();

    int getNumberEdges();

    bool getDirected();

    bool getWeightedEdge();
    vector<tuple<int, int>> define_leitura();
    bool getWeightedNode();

    Vertice *getFirstNode();

    Vertice *getLastNode();

    //Other methods
    void insertNode(int id);
    void insertNodeAndWeight(int id, int weight);

    void insertEdge(int id, int target_id, float weight);

    void removeNode(int id);

    bool searchNode(int id);

    Vertice *getNode(int id);

    bool searchEdge(int id, int target_id);

    //methods phase1
    vector<Grafo*> guloso(bool random, double *result, float alfa);
    // vector<Graph*> gulosoRand(bool random, float *result, float alfa);
    void agmGuloso();
    void gulosoRandAdap();
    void gulosoReativo();

    float qualidadeSolucao(float resultadoObtido);
    Vertice* verticeValido(float min, float max);
    float distancia(int node1, int node2);
    void qualidade(float result);
    void imprimeCluster(vector<Grafo*> solucao, int option, float resultBeneficio);
    void imprimirVertices();
    void imprimirVertices2();

    //axiliar methods
    void printGraph();
    void printGraphDot(ofstream &file);
    void output(string output_file, vector<Grafo*> solution, float quality);
};

#endif // GRAFO_H