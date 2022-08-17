/*
    CABEÇALHO
*/

#ifndef VERTICE_H
#define VERTICE_H 1

#include "Aresta.h"

class Vertice {

    // Attributes
    private:
        Aresta* first_edge;
        Aresta* last_edge;
        int id;
        unsigned int in_degree;
        unsigned int out_degree;
        float weight;
        bool visited;
        Vertice* next_node;

    public:
        // Constructor
        Vertice(int id);
        // Destructor
        ~Vertice();
        // Getters
        Aresta* getFirstEdge();
        Aresta* getLastEdge();
        int getId();
        int getInDegree();
        int getOutDegree();
        float getWeight();
        Vertice* getNextNode();
        Aresta* getEdge(int id_alvo);
        bool getVisited();
        // Setters
        void setNextNode(Vertice* node);
        void setWeight(float weight);
        void setVisited(bool visited);
        // Other methods
        bool searchEdge(int target_id);
        void insertEdge(int target_id, float weight);
        void removeAllEdges();
        int removeEdge(int id, bool directed, Vertice* target_node);
        void incrementOutDegree();
        void decrementOutDegree();
        void incrementInDegree();
        void decrementInDegree();
        Aresta* hasEdgeBetween(int target_id);
        bool verifyEdge(int target_id);
        // Auxiliar methods
};

#endif // VERTICE_H