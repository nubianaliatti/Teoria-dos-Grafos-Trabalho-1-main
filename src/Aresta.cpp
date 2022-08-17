#include <iostream>
#include <fstream>

#include "Aresta.h"
#include "Vertice.h"

using namespace std;

Aresta::Aresta(int target_id)
{
    this->target_id = target_id;
    this->next_edge = nullptr;
    this->weight = 0;
}

// Destructor
Aresta::~Aresta()
{
    if (this->next_edge != nullptr)
    {
        delete this->next_edge;
        this->next_edge = nullptr;
    }
}

// Getters
int Aresta::getTargetId()
{
    return this->target_id;
}

int Aresta::getOrigem()
{
    return this->origem;
}

Aresta *Aresta::getNextEdge()
{
    return this->next_edge;
}

float Aresta::getWeight()
{
    return this->weight;
}

// Setters
void Aresta::setNextEdge(Aresta *edge)
{
    this->next_edge = edge;
}

void Aresta::setWeight(float weight)
{
    this->weight = weight;
}

void Aresta::setOrigem(int origem)
{
    this->origem = origem;
}