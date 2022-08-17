#include <iostream>
#include <fstream>

#include "Vertice.h"

using namespace std;

// Constructor
Vertice::Vertice(int id) {
    this->id = id;
    this->in_degree = 0;
    this->out_degree = 0;
    this->weight = 0;
    this->visited = false;
    this->first_edge = nullptr;
    this->last_edge = nullptr;
    this->next_node = nullptr;
};

// Destructor
Vertice::~Vertice() {
    Aresta* next_edge = this->first_edge;
    while(next_edge != nullptr) {
        Aresta* aux_edge = next_edge->getNextEdge();
        delete next_edge;
        next_edge = aux_edge;
    }
};

// Getters
Aresta* Vertice::getFirstEdge() {
    return this->first_edge;
}

Aresta* Vertice::getLastEdge() {
    return this->last_edge;
}

int Vertice::getId() {
    return this->id;
}

int Vertice::getInDegree() {
    return this->in_degree;
}

int Vertice::getOutDegree() {
    return this->out_degree;
}
bool Vertice::getVisited(){
    return visited;
}

float Vertice::getWeight() {
    return this->weight;
}

Vertice* Vertice::getNextNode() {
    return this->next_node;
}
Aresta* Vertice::getEdge(int id)
{
    for(Aresta *aresta_auxiliar = this->first_edge; aresta_auxiliar != nullptr; aresta_auxiliar = aresta_auxiliar->getNextEdge())
    {
        if(aresta_auxiliar->getTargetId() == id)
            return aresta_auxiliar;
    }
    return nullptr;
}

// Setters

void Vertice::setNextNode(Vertice* next_node) {
    this->next_node = next_node;
}

void Vertice::setWeight(float weight) {
    this->weight = weight;
}
void Vertice::setVisited(bool visitado) {
    this->visited = visitado;
}

// Other methods
void Vertice::insertEdge(int target_id, float weight) {
    // Verifies whether there are at least one edge in the node
    if(this->first_edge != nullptr) {
        // Allocating the new edge and keeping the integrity of the edge list
        Aresta* edge = new Aresta(target_id);
        edge->setWeight(weight);
        this->last_edge->setNextEdge(edge);
        this->last_edge = edge;
    } else {
         // Allocating the new edge and keeping the integrity of the edge list
        this->first_edge = new Aresta(target_id);
        this->first_edge->setWeight(weight);
        this->last_edge = this->first_edge;
    }
}

void Vertice::removeAllEdges() {
    // Verifies whether there are at least one edge in the node
    if(this->first_edge != nullptr) {
        Aresta* next = nullptr;
        Aresta* aux = this->first_edge;
        // Removing all edges of the node
        while(aux != nullptr) {
            next = aux->getNextEdge();
            delete aux;
        }
    }
    this->first_edge = this->last_edge = nullptr;
}

int Vertice::removeEdge(int id, bool directed, Vertice* target_node) {
    // Verifies whether the edge to remove is in the node
    if(this->searchEdge(id)) {
        Aresta* aux = this->first_edge;
        Aresta* previous = nullptr;
        // Searching for the edge to be removed
        while(aux->getTargetId() != id) {
            previous = aux;
            aux = aux->getNextEdge();
        }
        // Keeping the integrity of the edge list
        if(previous != nullptr)
            previous->setNextEdge(aux->getNextEdge());

        else
            this->first_edge = aux->getNextEdge();

        if(aux == this->last_edge)
            this->last_edge = previous;

        if(aux->getNextEdge() == this->last_edge)
            this->last_edge = aux->getNextEdge();

        delete aux;
        // Verifies whether the graph is directed
        if(directed)
            this->decrementOutDegree();

        else{

            this->decrementInDegree();
            target_node->decrementInDegree();

        }

        return 1;

    }

    return 0;

}

bool Vertice::searchEdge(int target_id) {
     // Verifies whether there are at least one edge in the node
    if(this->first_edge != nullptr) {
        // Searching for a specific edge of target id equal to target id
        for(Aresta* aux = this->first_edge; aux != nullptr; aux = aux->getNextEdge())
            if(aux->getTargetId() == target_id)
                return true;
    }
    return false;
}

void Vertice::incrementInDegree() {
    this->in_degree++;
}

void Vertice::incrementOutDegree() {
    this->out_degree++;
}

void Vertice::decrementInDegree() {
    this->in_degree--;
}

void Vertice::decrementOutDegree() {
    this->out_degree--;
}

Aresta* Vertice::hasEdgeBetween(int target_id) {
    for(Aresta *auxEdge = this->first_edge; auxEdge != nullptr; auxEdge = auxEdge->getNextEdge()) {
        if(auxEdge->getTargetId() == target_id)
            return auxEdge;
    }
    return nullptr;
}
bool Vertice::verifyEdge(int target_id) {
    for(Aresta *auxEdge = this->first_edge; auxEdge != nullptr; auxEdge = auxEdge->getNextEdge()) {
        if(auxEdge->getTargetId() == target_id)
            return true;
    }
    return false;
}